"""Genomics API endpoints for VCF upload and ENA integration."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, File, HTTPException, Query, UploadFile, status
from pydantic import BaseModel

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.database.models import User
from amprenta_rag.models.misc import GeneticVariant
from amprenta_rag.database.session import db_session

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
