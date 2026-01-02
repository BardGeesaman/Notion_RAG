"""Genomics API endpoints for VCF upload and ENA integration."""

from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, File, HTTPException, Query, UploadFile, status
from pydantic import BaseModel

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.database.models import User
from amprenta_rag.models.misc import AlignmentFile, GeneticVariant
from amprenta_rag.database.session import db_session
from amprenta_rag.ingestion.genomics.bam_parser import (
    parse_bam_header,
    fetch_reads,
    get_coverage,
    check_index_exists,
    ALIGNMENT_STORAGE_PATH,
)

router = APIRouter(prefix="/genomics", tags=["genomics"])

MAX_VCF_SIZE = 50_000_000  # 50MB


class VCFUploadResponse(BaseModel):
    """Response for VCF upload."""
    success: bool
    variants_count: int
    variants: list[dict] = []  # For preview mode
    imported_count: int = 0
    error: Optional[str] = None


class ENAStudyResult(BaseModel):
    """ENA study search result."""
    run_accession: str
    study_accession: Optional[str]
    title: Optional[str]
    organism: Optional[str]
    platform: Optional[str]
    fastq_count: int = 0


class ENASearchResponse(BaseModel):
    """Response for ENA search."""
    success: bool
    query: str
    results: List[ENAStudyResult]
    total: int
    error: Optional[str] = None


class ENAIngestRequest(BaseModel):
    """Request for ENA ingestion."""
    study_ids: List[str]


class ENAIngestResponse(BaseModel):
    """Response for ENA ingestion."""
    success: bool
    job_id: Optional[str] = None
    message: str


# --- Alignment Schemas ---

class AlignmentUploadResponse(BaseModel):
    """Response for alignment upload."""
    success: bool
    alignment_id: Optional[UUID] = None
    filename: str
    file_format: str
    has_index: bool
    stats: Optional[dict] = None
    error: Optional[str] = None


class AlignmentListItem(BaseModel):
    """Alignment file list item."""
    id: UUID
    filename: str
    file_format: str
    has_index: bool
    reference_genome: Optional[str]
    total_reads: Optional[int]
    created_at: datetime


class AlignmentDetail(BaseModel):
    """Full alignment file details."""
    id: UUID
    filename: str
    file_format: str
    has_index: bool
    reference_genome: Optional[str]
    num_references: Optional[int]
    read_groups: Optional[list]
    total_reads: Optional[int]
    mapped_reads: Optional[int]
    unmapped_reads: Optional[int]
    duplicate_rate: Optional[float]
    mean_coverage: Optional[float]
    experiment_id: Optional[UUID]
    created_at: datetime


class ReadResponse(BaseModel):
    """Response for reads endpoint."""
    success: bool
    region: str
    reads: list[dict]
    total_in_region: int
    page: int
    page_size: int
    error: Optional[str] = None


class CoverageResponse(BaseModel):
    """Response for coverage endpoint."""
    success: bool
    region: str
    bin_size: int
    coverage: list[dict]
    error: Optional[str] = None


@router.post("/variants/upload", response_model=VCFUploadResponse)
async def upload_vcf(
    file: UploadFile = File(...),
    experiment_id: Optional[UUID] = Query(None),
    preview_only: bool = Query(False, description="Parse without saving"),
    current_user: User = Depends(get_current_user),
) -> VCFUploadResponse:
    """
    Upload and parse VCF file.
    
    - preview_only=True: Parse and return variants without saving
    - preview_only=False: Parse and create GeneticVariant records
    """
    # Validate file type
    if not file.filename or not file.filename.endswith((".vcf", ".vcf.gz")):
        raise HTTPException(
            status_code=status.HTTP_415_UNSUPPORTED_MEDIA_TYPE,
            detail="File must be .vcf or .vcf.gz",
        )
    
    # Read and validate size
    content = await file.read()
    if len(content) > MAX_VCF_SIZE:
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail=f"File too large. Maximum: {MAX_VCF_SIZE // 1_000_000}MB",
        )
    
    try:
        from amprenta_rag.ingestion.genomics.vcf_parser import parse_vcf_bytes
        
        variants = parse_vcf_bytes(content)
        
        if preview_only:
            # Return first 100 for preview
            return VCFUploadResponse(
                success=True,
                variants_count=len(variants),
                variants=variants[:100],
            )
        
        # Import to database
        imported = 0
        with db_session() as db:
            for v in variants:
                variant = GeneticVariant(
                    gene=v.get("gene") or "Unknown",
                    variant=v.get("variant", ""),
                    organism="",  # User can update later
                    chromosome=v.get("chromosome"),
                    position=v.get("position"),
                    reference_allele=v.get("reference_allele"),
                    alternate_allele=v.get("alternate_allele"),
                    quality=v.get("quality"),
                    vcf_filter=v.get("vcf_filter"),
                    experiment_id=experiment_id,
                )
                db.add(variant)
                imported += 1
            
            db.commit()
        
        return VCFUploadResponse(
            success=True,
            variants_count=len(variants),
            imported_count=imported,
        )
    
    except Exception as e:
        return VCFUploadResponse(
            success=False,
            variants_count=0,
            error=str(e),
        )


@router.get("/ena/search", response_model=ENASearchResponse)
def search_ena_studies(
    q: str = Query(..., description="Search keyword"),
    organism: Optional[str] = Query(None, description="Filter by organism"),
    library_strategy: Optional[str] = Query(None, description="e.g., RNA-Seq, WGS"),
    max_results: int = Query(20, ge=1, le=100),
    current_user: User = Depends(get_current_user),
) -> ENASearchResponse:
    """Search ENA for genomics studies."""
    try:
        from amprenta_rag.ingestion.repositories.ena import ENARepository
        
        ena = ENARepository()
        
        # Build filters
        filters = {}
        if organism:
            filters["organism"] = organism
        if library_strategy:
            filters["library_strategy"] = library_strategy
        
        # Search
        study_ids = ena.search_studies(
            keywords=[q],
            filters=filters if filters else None,
            max_results=max_results,
        )
        
        # Fetch metadata for results
        results = []
        for study_id in study_ids[:max_results]:
            metadata = ena.fetch_study_metadata(study_id)
            if metadata:
                results.append(ENAStudyResult(
                    run_accession=study_id,
                    study_accession=metadata.raw_metadata.get("study_accession"),
                    title=metadata.title,
                    organism=metadata.organism[0] if metadata.organism else None,
                    platform=metadata.platform,
                    fastq_count=len(metadata.raw_metadata.get("fastq_ftp", "").split(";")) if metadata.raw_metadata.get("fastq_ftp") else 0,
                ))
        
        return ENASearchResponse(
            success=True,
            query=q,
            results=results,
            total=len(results),
        )
    
    except Exception as e:
        return ENASearchResponse(
            success=False,
            query=q,
            results=[],
            total=0,
            error=str(e),
        )


@router.post("/ena/ingest", response_model=ENAIngestResponse)
def ingest_ena_studies(
    request: ENAIngestRequest,
    current_user: User = Depends(get_current_user),
) -> ENAIngestResponse:
    """Trigger ENA study ingestion (creates Celery job)."""
    if not request.study_ids:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="No study IDs provided",
        )
    
    # For now, return a placeholder - Celery task integration can be added later
    # The existing genomics Celery tasks can be extended
    return ENAIngestResponse(
        success=True,
        job_id=None,  # Would be Celery task ID
        message=f"Queued {len(request.study_ids)} studies for ingestion (Celery integration pending)",
    )


# --- Alignment Endpoints ---

MAX_ALIGNMENT_SIZE = 500_000_000  # 500MB


@router.post("/alignments/upload", response_model=AlignmentUploadResponse)
async def upload_alignment(
    file: UploadFile = File(...),
    index_file: UploadFile = File(None),
    experiment_id: Optional[UUID] = Query(None),
    current_user: User = Depends(get_current_user),
) -> AlignmentUploadResponse:
    """
    Upload BAM/CRAM alignment file with optional index.
    
    - file: BAM (.bam) or CRAM (.cram) file
    - index_file: Optional .bai or .crai index file
    """
    import uuid as uuid_module
    
    # Validate file type
    if not file.filename:
        raise HTTPException(status_code=400, detail="Filename required")
    
    filename_lower = file.filename.lower()
    if filename_lower.endswith(".bam"):
        file_format = "BAM"
    elif filename_lower.endswith(".cram"):
        file_format = "CRAM"
    else:
        raise HTTPException(
            status_code=status.HTTP_415_UNSUPPORTED_MEDIA_TYPE,
            detail="File must be .bam or .cram",
        )
    
    # Read and validate size
    content = await file.read()
    if len(content) > MAX_ALIGNMENT_SIZE:
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail=f"File too large. Maximum: {MAX_ALIGNMENT_SIZE // 1_000_000}MB",
        )
    
    # Validate magic bytes
    if file_format == "BAM" and not content[:4] == b'\x1f\x8b\x08\x04':  # gzip magic for BAM
        raise HTTPException(status_code=400, detail="Invalid BAM file format")
    elif file_format == "CRAM" and not content[:4] == b'CRAM':
        raise HTTPException(status_code=400, detail="Invalid CRAM file format")
    
    try:
        # Create storage directory
        alignment_id = uuid_module.uuid4()
        storage_dir = Path(ALIGNMENT_STORAGE_PATH) / str(alignment_id)
        storage_dir.mkdir(parents=True, exist_ok=True)
        
        # Save main file
        file_path = storage_dir / file.filename
        with open(file_path, "wb") as f:
            f.write(content)
        
        # Save index file if provided
        index_path = None
        has_index = False
        if index_file and index_file.filename:
            index_content = await index_file.read()
            index_path = storage_dir / index_file.filename
            with open(index_path, "wb") as f:
                f.write(index_content)
            has_index = True
        else:
            # Check if index exists alongside file
            has_index, existing_index = check_index_exists(file_path)
            if existing_index:
                index_path = existing_index
        
        # Parse header and stats
        header_info = parse_bam_header(file_path)
        
        # Save to database
        with db_session() as db:
            alignment = AlignmentFile(
                id=alignment_id,
                filename=file.filename,
                file_path=str(file_path),
                index_file_path=str(index_path) if index_path else None,
                file_format=file_format,
                has_index=has_index,
                reference_genome=header_info.get("reference_genome"),
                num_references=header_info.get("num_references"),
                read_groups=header_info.get("read_groups"),
                experiment_id=experiment_id,
                created_by_id=current_user.id,
            )
            db.add(alignment)
            db.commit()
        
        return AlignmentUploadResponse(
            success=True,
            alignment_id=alignment_id,
            filename=file.filename,
            file_format=file_format,
            has_index=has_index,
            stats=header_info,
        )
    
    except Exception as e:
        # Cleanup on failure
        if 'storage_dir' in locals() and storage_dir.exists():
            shutil.rmtree(storage_dir, ignore_errors=True)
        return AlignmentUploadResponse(
            success=False,
            filename=file.filename or "",
            file_format=file_format if 'file_format' in locals() else "UNKNOWN",
            has_index=False,
            error=str(e),
        )


@router.get("/alignments", response_model=list[AlignmentListItem])
def list_alignments(
    experiment_id: Optional[UUID] = Query(None),
    has_index: Optional[bool] = Query(None),
    limit: int = Query(100, ge=1, le=500),
    current_user: User = Depends(get_current_user),
) -> list[AlignmentListItem]:
    """List alignment files with optional filters."""
    with db_session() as db:
        query = db.query(AlignmentFile)
        
        if experiment_id:
            query = query.filter(AlignmentFile.experiment_id == experiment_id)
        if has_index is not None:
            query = query.filter(AlignmentFile.has_index == has_index)
        
        alignments = query.order_by(AlignmentFile.created_at.desc()).limit(limit).all()
        
        return [
            AlignmentListItem(
                id=a.id,
                filename=a.filename,
                file_format=a.file_format,
                has_index=a.has_index,
                reference_genome=a.reference_genome,
                total_reads=a.total_reads,
                created_at=a.created_at,
            )
            for a in alignments
        ]


@router.get("/alignments/{alignment_id}", response_model=AlignmentDetail)
def get_alignment(
    alignment_id: UUID,
    current_user: User = Depends(get_current_user),
) -> AlignmentDetail:
    """Get alignment file details."""
    with db_session() as db:
        alignment = db.query(AlignmentFile).filter(
            AlignmentFile.id == alignment_id
        ).first()
        
        if not alignment:
            raise HTTPException(status_code=404, detail="Alignment not found")
        
        return AlignmentDetail(
            id=alignment.id,
            filename=alignment.filename,
            file_format=alignment.file_format,
            has_index=alignment.has_index,
            reference_genome=alignment.reference_genome,
            num_references=alignment.num_references,
            read_groups=alignment.read_groups,
            total_reads=alignment.total_reads,
            mapped_reads=alignment.mapped_reads,
            unmapped_reads=alignment.unmapped_reads,
            duplicate_rate=alignment.duplicate_rate,
            mean_coverage=alignment.mean_coverage,
            experiment_id=alignment.experiment_id,
            created_at=alignment.created_at,
        )


@router.get("/alignments/{alignment_id}/reads", response_model=ReadResponse)
def get_alignment_reads(
    alignment_id: UUID,
    region: str = Query(..., description="Region in chr:start-end format"),
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=1000),
    current_user: User = Depends(get_current_user),
) -> ReadResponse:
    """
    Fetch reads from a genomic region.
    
    Requires index file (.bai/.crai) for the alignment.
    """
    with db_session() as db:
        alignment = db.query(AlignmentFile).filter(
            AlignmentFile.id == alignment_id
        ).first()
        
        if not alignment:
            raise HTTPException(status_code=404, detail="Alignment not found")
        
        if not alignment.has_index:
            raise HTTPException(
                status_code=400,
                detail="Index file required for region queries. Upload .bai/.crai file.",
            )
        
        file_path = Path(alignment.file_path)
        if not file_path.exists():
            raise HTTPException(status_code=404, detail="Alignment file not found on disk")
    
    try:
        # Validate region format
        import re
        if not re.match(r'^[\w]+:\d+-\d+$', region):
            raise HTTPException(
                status_code=400,
                detail="Invalid region format. Use chr:start-end (e.g., chr1:1000-2000)",
            )
        
        offset = (page - 1) * page_size
        reads = fetch_reads(file_path, region, offset=offset, limit=page_size)
        
        return ReadResponse(
            success=True,
            region=region,
            reads=[r.model_dump() for r in reads],
            total_in_region=len(reads),  # Note: actual total requires separate count
            page=page,
            page_size=page_size,
        )
    
    except ValueError as e:
        return ReadResponse(
            success=False,
            region=region,
            reads=[],
            total_in_region=0,
            page=page,
            page_size=page_size,
            error=str(e),
        )


@router.get("/alignments/{alignment_id}/coverage", response_model=CoverageResponse)
def get_alignment_coverage(
    alignment_id: UUID,
    region: str = Query(..., description="Region in chr:start-end format"),
    bin_size: int = Query(100, ge=10, le=10000),
    current_user: User = Depends(get_current_user),
) -> CoverageResponse:
    """Get coverage histogram for a genomic region."""
    with db_session() as db:
        alignment = db.query(AlignmentFile).filter(
            AlignmentFile.id == alignment_id
        ).first()
        
        if not alignment:
            raise HTTPException(status_code=404, detail="Alignment not found")
        
        if not alignment.has_index:
            raise HTTPException(
                status_code=400,
                detail="Index file required for coverage queries.",
            )
        
        file_path = Path(alignment.file_path)
        if not file_path.exists():
            raise HTTPException(status_code=404, detail="Alignment file not found on disk")
    
    try:
        import re
        if not re.match(r'^[\w]+:\d+-\d+$', region):
            raise HTTPException(
                status_code=400,
                detail="Invalid region format. Use chr:start-end",
            )
        
        coverage_data = get_coverage(file_path, region, bin_size=bin_size)
        
        return CoverageResponse(
            success=True,
            region=region,
            bin_size=bin_size,
            coverage=coverage_data,
        )
    
    except ValueError as e:
        return CoverageResponse(
            success=False,
            region=region,
            bin_size=bin_size,
            coverage=[],
            error=str(e),
        )
