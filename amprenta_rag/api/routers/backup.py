"""Backup and disaster recovery API endpoints."""

import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, List, Optional
from uuid import UUID, uuid4

from fastapi import APIRouter, Depends, HTTPException, Query, status
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.backup.backup_engine import BackupEngine
from amprenta_rag.backup.project_export import export_project
from amprenta_rag.database.models import BackupRecord, ProjectExport
from amprenta_rag.models.auth import User

router = APIRouter(prefix="/backup", tags=["Backup"])


class TriggerBackupRequest(BaseModel):
    """Request to trigger a manual backup."""
    backup_type: str = "full"


class ProjectExportRequest(BaseModel):
    """Request to create a project export."""
    program_ids: Optional[List[UUID]] = None
    experiment_ids: Optional[List[UUID]] = None
    compound_ids: Optional[List[UUID]] = None
    include_related: bool = True


class BackupResponse(BaseModel):
    """Backup record response."""
    id: UUID
    backup_type: str
    status: str
    file_path: Optional[str]
    file_size_bytes: Optional[int]
    checksum_sha256: Optional[str]
    started_at: Optional[datetime]
    completed_at: Optional[datetime]
    error_message: Optional[str]
    created_at: datetime


class BackupHistoryResponse(BaseModel):
    """Backup history response."""
    items: List[BackupResponse]
    total: int
    page: int
    per_page: int


class ProjectExportResponse(BaseModel):
    """Project export creation response."""
    export_id: UUID
    message: str
    export_size_bytes: int
    entities_summary: str
    expires_at: datetime


@router.post("/database", status_code=status.HTTP_202_ACCEPTED)
async def trigger_manual_backup(
    request: TriggerBackupRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> Dict[str, str]:
    """Trigger a manual database backup (admin only).
    
    Args:
        request: Backup request parameters
        current_user: Current authenticated user
        db: Database session
        
    Returns:
        Dict with task information
        
    Raises:
        HTTPException: If user is not authorized
    """
    # TODO: Add admin role check when role system is implemented
    # For now, any authenticated user can trigger backups
    
    if request.backup_type not in ("full", "incremental"):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="backup_type must be 'full' or 'incremental'"
        )
    
    # Trigger backup task
    task = run_database_backup.delay(request.backup_type)
    
    return {
        "message": f"{request.backup_type.title()} backup initiated",
        "task_id": task.id,
        "backup_type": request.backup_type,
    }


@router.get("/history")
async def get_backup_history(
    page: int = Query(1, ge=1, description="Page number"),
    per_page: int = Query(50, ge=1, le=200, description="Items per page"),
    backup_type: Optional[str] = Query(None, description="Filter by backup type"),
    status_filter: Optional[str] = Query(None, alias="status", description="Filter by status"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> BackupHistoryResponse:
    """Get backup history with pagination and filtering.
    
    Args:
        page: Page number (1-based)
        per_page: Items per page
        backup_type: Filter by backup type
        status_filter: Filter by status
        current_user: Current authenticated user
        db: Database session
        
    Returns:
        Paginated backup history
    """
    # Build query
    query = db.query(BackupRecord)
    
    if backup_type:
        query = query.filter(BackupRecord.backup_type == backup_type)
    
    if status_filter:
        query = query.filter(BackupRecord.status == status_filter)
    
    # Get total count
    total = query.count()
    
    # Apply pagination
    offset = (page - 1) * per_page
    backups = query.order_by(BackupRecord.created_at.desc()).offset(offset).limit(per_page).all()
    
    # Convert to response format
    items = [
        BackupResponse(
            id=backup.id,
            backup_type=backup.backup_type,
            status=backup.status,
            file_path=backup.file_path,
            file_size_bytes=backup.file_size_bytes,
            checksum_sha256=backup.checksum_sha256,
            started_at=backup.started_at,
            completed_at=backup.completed_at,
            error_message=backup.error_message,
            created_at=backup.created_at,
        )
        for backup in backups
    ]
    
    return BackupHistoryResponse(
        items=items,
        total=total,
        page=page,
        per_page=per_page,
    )


@router.get("/{backup_id}")
async def get_backup_details(
    backup_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> BackupResponse:
    """Get backup details by ID.
    
    Args:
        backup_id: Backup UUID
        current_user: Current authenticated user
        db: Database session
        
    Returns:
        Backup details
        
    Raises:
        HTTPException: If backup not found
    """
    backup = db.query(BackupRecord).filter(BackupRecord.id == backup_id).first()
    
    if not backup:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Backup {backup_id} not found"
        )
    
    return BackupResponse(
        id=backup.id,
        backup_type=backup.backup_type,
        status=backup.status,
        file_path=backup.file_path,
        file_size_bytes=backup.file_size_bytes,
        checksum_sha256=backup.checksum_sha256,
        started_at=backup.started_at,
        completed_at=backup.completed_at,
        error_message=backup.error_message,
        created_at=backup.created_at,
    )


@router.get("/{backup_id}/download")
async def download_backup(
    backup_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> StreamingResponse:
    """Download backup file (streaming).
    
    Args:
        backup_id: Backup UUID
        current_user: Current authenticated user
        db: Database session
        
    Returns:
        Streaming response with backup file
        
    Raises:
        HTTPException: If backup not found or not completed
    """
    backup = db.query(BackupRecord).filter(BackupRecord.id == backup_id).first()
    
    if not backup:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Backup {backup_id} not found"
        )
    
    if backup.status != "completed":
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Backup is not completed (status: {backup.status})"
        )
    
    if not backup.file_path:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Backup file path not available"
        )
    
    # Download backup file to temporary location
    engine = BackupEngine()
    
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_path = Path(temp_file.name)
    
    try:
        # Download from storage
        engine.backup_client.download_file(backup.file_path, str(temp_path))
        
        # Create streaming response
        def file_generator():
            with open(temp_path, 'rb') as f:
                while chunk := f.read(8192):
                    yield chunk
        
        # Generate filename
        timestamp = backup.created_at.strftime("%Y%m%d_%H%M%S")
        filename = f"backup_{backup.backup_type}_{timestamp}.sql.gz"
        
        return StreamingResponse(
            file_generator(),
            media_type="application/gzip",
            headers={"Content-Disposition": f"attachment; filename={filename}"}
        )
        
    finally:
        # Clean up temp file
        if temp_path.exists():
            temp_path.unlink()


@router.post("/export", status_code=status.HTTP_201_CREATED)
async def create_project_export(
    request: ProjectExportRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> ProjectExportResponse:
    """Create a project export ZIP file.
    
    Args:
        request: Export request parameters
        current_user: Current authenticated user
        db: Database session
        
    Returns:
        Export information with download ID
        
    Raises:
        HTTPException: If no entities specified for export
    """
    # Validate request
    if not any([request.program_ids, request.experiment_ids, request.compound_ids]):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="At least one of program_ids, experiment_ids, or compound_ids must be specified"
        )
    
    try:
        # Generate export
        export_data = export_project(
            program_ids=request.program_ids,
            experiment_ids=request.experiment_ids,
            compound_ids=request.compound_ids,
            db=db,
            include_related=request.include_related,
        )
        
        export_size = len(export_data)
        entities_summary = f"programs: {len(request.program_ids or [])}, experiments: {len(request.experiment_ids or [])}, compounds: {len(request.compound_ids or [])}"
        
        # Store export using backup client
        engine = BackupEngine()
        export_id = uuid4()
        file_path = f"exports/{export_id}.zip"
        
        # Create temporary file and upload to storage
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_file.write(export_data)
            temp_path = Path(temp_file.name)
        
        try:
            # Upload to storage
            engine.backup_client.upload_file(str(temp_path), file_path)
            
            # Create ProjectExport record with 24-hour expiration
            expires_at = datetime.now(timezone.utc) + timedelta(hours=24)
            
            project_export = ProjectExport(
                id=export_id,
                file_path=file_path,
                file_size_bytes=export_size,
                entities_summary=entities_summary,
                created_by=current_user.id,
                expires_at=expires_at,
            )
            
            db.add(project_export)
            db.commit()
            
            return ProjectExportResponse(
                export_id=export_id,
                message="Project export created successfully",
                export_size_bytes=export_size,
                entities_summary=entities_summary,
                expires_at=expires_at,
            )
            
        finally:
            # Clean up temp file
            if temp_path.exists():
                temp_path.unlink()
        
    except Exception as e:
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Export failed: {str(e)}"
        )


@router.get("/export/{export_id}")
async def download_project_export(
    export_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> StreamingResponse:
    """Download project export ZIP file.
    
    Args:
        export_id: Export UUID
        current_user: Current authenticated user
        db: Database session
        
    Returns:
        Streaming response with export ZIP
        
    Raises:
        HTTPException: If export not found or expired
    """
    # Lookup ProjectExport by ID
    project_export = db.query(ProjectExport).filter(ProjectExport.id == export_id).first()
    
    if not project_export:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Export {export_id} not found"
        )
    
    # Check expiration (return 410 GONE if expired)
    if datetime.now(timezone.utc) > project_export.expires_at:
        # Clean up expired record
        db.delete(project_export)
        db.commit()
        
        raise HTTPException(
            status_code=status.HTTP_410_GONE,
            detail=f"Export {export_id} has expired and is no longer available"
        )
    
    # Download from storage and stream response
    engine = BackupEngine()
    
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_path = Path(temp_file.name)
    
    try:
        # Download from storage
        engine.backup_client.download_file(project_export.file_path, str(temp_path))
        
        # Create streaming response
        def file_generator():
            with open(temp_path, 'rb') as f:
                while chunk := f.read(8192):
                    yield chunk
        
        # Generate filename
        timestamp = project_export.created_at.strftime("%Y%m%d_%H%M%S")
        filename = f"project_export_{timestamp}.zip"
        
        # Delete record after successful download (one-time use)
        db.delete(project_export)
        db.commit()
        
        return StreamingResponse(
            file_generator(),
            media_type="application/zip",
            headers={"Content-Disposition": f"attachment; filename={filename}"}
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to download export: {str(e)}"
        )
        
    finally:
        # Clean up temp file
        if temp_path.exists():
            temp_path.unlink()
