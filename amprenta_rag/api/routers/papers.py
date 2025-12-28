"""
API router for scientific papers ingestion and search.

Provides endpoints for searching PubMed/PMC, ingesting papers with
deduplication, and retrieving paper content.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel, ConfigDict, Field
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.database.models import Literature
from amprenta_rag.ingestion.papers import PubMedRepository
from amprenta_rag.ingestion.papers.embedding import embed_paper_sections
from amprenta_rag.ingestion.papers.jats_parser import PaperSection, parse_jats_xml
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

router = APIRouter()


# Pydantic schemas
class PaperSearchRequest(BaseModel):
    """Request for paper search."""

    query: str = Field(..., min_length=1, description="Search query")
    source: str = Field(default="pubmed", description="Paper source (pubmed, biorxiv)")
    limit: int = Field(default=20, ge=1, le=100, description="Maximum results")
    offset: int = Field(default=0, ge=0, description="Pagination offset")


class PaperSearchResult(BaseModel):
    """Single search result."""

    paper_id: str
    title: str
    abstract: Optional[str] = None
    authors: List[str] = Field(default_factory=list)
    journal: Optional[str] = None
    year: Optional[int] = None
    doi: Optional[str] = None
    pmid: Optional[str] = None
    url: Optional[str] = None

    model_config = ConfigDict(from_attributes=True)


class PaperSearchResponse(BaseModel):
    """Response for paper search."""

    results: List[PaperSearchResult]
    total: int
    offset: int
    limit: int


class PaperIngestRequest(BaseModel):
    """Request for paper ingestion."""

    pmid: Optional[str] = Field(None, description="PubMed ID")
    doi: Optional[str] = Field(None, description="DOI")
    fetch_fulltext: bool = Field(default=False, description="Attempt to fetch full text")


class PaperIngestResponse(BaseModel):
    """Response for paper ingestion."""

    literature_id: UUID
    pmid: Optional[str] = None
    doi: Optional[str] = None
    title: str
    already_exists: bool = Field(description="True if paper was already in database")
    chunks_created: int = Field(default=0, description="Number of RAG chunks created")


class PaperSectionResponse(BaseModel):
    """Paper section for response."""

    title: str
    content: str
    order: int

    model_config = ConfigDict(from_attributes=True)


class PaperDetailResponse(BaseModel):
    """Detailed paper response with sections."""

    id: UUID
    pmid: Optional[str] = None
    pmc_id: Optional[str] = None
    doi: Optional[str] = None
    title: str
    abstract: Optional[str] = None
    journal: Optional[str] = None
    year: Optional[int] = None
    mesh_terms: List[str] = Field(default_factory=list)
    sections: List[PaperSectionResponse] = Field(default_factory=list)

    model_config = ConfigDict(from_attributes=True)


@router.post("/search", response_model=PaperSearchResponse)
async def search_papers(
    request: PaperSearchRequest,
    db: Session = Depends(get_database_session),
) -> PaperSearchResponse:
    """
    Search for scientific papers.

    Currently supports PubMed search. Returns metadata for matching papers.
    """
    if request.source.lower() != "pubmed":
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported source: {request.source}. Only 'pubmed' is currently supported.",
        )

    try:
        repo = PubMedRepository()
        pmids = repo.search_papers(request.query, max_results=request.limit + request.offset)

        # Apply pagination
        paginated_pmids = pmids[request.offset : request.offset + request.limit]

        # Fetch metadata for each result
        results = []
        for pmid in paginated_pmids:
            metadata = repo.fetch_metadata(pmid)
            if metadata:
                results.append(
                    PaperSearchResult(
                        paper_id=metadata.paper_id,
                        title=metadata.title,
                        abstract=metadata.abstract,
                        authors=metadata.authors,
                        journal=metadata.journal,
                        year=metadata.year,
                        doi=metadata.doi,
                        pmid=metadata.pmid,
                        url=metadata.url,
                    )
                )

        return PaperSearchResponse(
            results=results,
            total=len(pmids),
            offset=request.offset,
            limit=request.limit,
        )

    except Exception as e:
        logger.error("[PAPERS_API] Search failed: %r", e)
        raise HTTPException(status_code=500, detail=f"Search failed: {str(e)}")


@router.post("/ingest", response_model=PaperIngestResponse, status_code=201)
async def ingest_paper(
    request: PaperIngestRequest,
    db: Session = Depends(get_database_session),
) -> PaperIngestResponse:
    """
    Ingest a paper by PMID or DOI.

    Checks for duplicates using DOI/PMID (P0-2 deduplication).
    Optionally fetches and embeds full text if available.
    """
    if not request.pmid and not request.doi:
        raise HTTPException(
            status_code=400, detail="Either pmid or doi must be provided"
        )

    # Check for existing paper (P0-2 deduplication)
    existing = None
    if request.pmid:
        existing = db.query(Literature).filter(Literature.pmid == request.pmid).first()
    if not existing and request.doi:
        existing = db.query(Literature).filter(Literature.doi == request.doi).first()

    if existing:
        logger.info(
            "[PAPERS_API] Paper already exists: PMID=%s, DOI=%s",
            request.pmid,
            request.doi,
        )
        return PaperIngestResponse(
            literature_id=existing.id,
            pmid=existing.pmid,
            doi=existing.doi,
            title=existing.title or "Unknown",
            already_exists=True,
            chunks_created=0,
        )

    # Fetch metadata from PubMed
    try:
        repo = PubMedRepository()

        # Search by PMID if provided, otherwise try DOI
        if request.pmid:
            metadata = repo.fetch_metadata(request.pmid)
        elif request.doi:
            # Search by DOI to get PMID, then fetch
            pmids = repo.search_papers(f'"{request.doi}"[DOI]', max_results=1)
            if not pmids:
                raise HTTPException(status_code=404, detail="Paper not found")
            metadata = repo.fetch_metadata(pmids[0])
        else:
            raise HTTPException(status_code=400, detail="PMID or DOI required")

        if not metadata:
            raise HTTPException(status_code=404, detail="Paper not found")

        # Create Literature record
        literature = Literature(
            title=metadata.title,
            abstract=metadata.abstract,
            doi=metadata.doi,
            journal=metadata.journal,
            year=metadata.year,
            pmid=metadata.pmid,
            pmc_id=metadata.pmc_id,
            mesh_terms=metadata.mesh_terms,
            source="pubmed",
            url=metadata.url,
            source_type="scientific_paper",
            full_text_available=False,  # Will update if full text fetched
        )
        db.add(literature)
        db.flush()  # Get literature.id

        logger.info(
            "[PAPERS_API] Created Literature record: PMID=%s, ID=%s",
            metadata.pmid,
            literature.id,
        )

        # TODO: Fetch and embed full text if requested and available
        # For now, just commit metadata
        chunks_created = 0
        if request.fetch_fulltext:
            logger.info(
                "[PAPERS_API] Full text fetching not yet implemented for PMID=%s",
                metadata.pmid,
            )

        db.commit()

        return PaperIngestResponse(
            literature_id=literature.id,
            pmid=literature.pmid,
            doi=literature.doi,
            title=literature.title,
            already_exists=False,
            chunks_created=chunks_created,
        )

    except HTTPException:
        raise
    except Exception as e:
        db.rollback()
        logger.error("[PAPERS_API] Ingestion failed: %r", e)
        raise HTTPException(status_code=500, detail=f"Ingestion failed: {str(e)}")


@router.get("/{paper_id}", response_model=PaperDetailResponse)
async def get_paper(
    paper_id: UUID,
    db: Session = Depends(get_database_session),
) -> PaperDetailResponse:
    """
    Get paper details by ID.

    Returns paper metadata and sections if available.
    """
    literature = db.query(Literature).filter(Literature.id == paper_id).first()

    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")

    # TODO: Retrieve sections from RAG chunks when full text is available
    sections = []

    return PaperDetailResponse(
        id=literature.id,
        pmid=literature.pmid,
        pmc_id=literature.pmc_id,
        doi=literature.doi,
        title=literature.title or "Unknown",
        abstract=literature.abstract,
        journal=literature.journal,
        year=literature.year,
        mesh_terms=literature.mesh_terms or [],
        sections=sections,
    )

