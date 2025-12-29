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
from amprenta_rag.ingestion.papers.semantic_scholar import SemanticScholarRepository
from amprenta_rag.ingestion.papers.openalex import OpenAlexRepository
from amprenta_rag.models.content import PaperCitation
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


class CitationResponse(BaseModel):
    """Citation information response."""

    paper_id: Optional[UUID] = None
    doi: Optional[str] = None
    title: str
    is_influential: bool = False
    citation_context: Optional[str] = None

    model_config = ConfigDict(from_attributes=True)


class EnrichPaperRequest(BaseModel):
    """Request for paper enrichment from external sources."""

    fetch_citations: bool = Field(default=True, description="Fetch and store citations")
    fetch_references: bool = Field(default=True, description="Fetch and store references")


class EnrichPaperResponse(BaseModel):
    """Response for paper enrichment."""

    enriched: bool
    citations_added: int = 0
    references_added: int = 0
    semantic_scholar_id: Optional[str] = None
    openalex_id: Optional[str] = None
    citation_count: Optional[int] = None


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


@router.get("/{paper_id}/citations", response_model=List[CitationResponse])
async def get_paper_citations(
    paper_id: UUID,
    limit: int = Query(50, ge=1, le=200),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_database_session),
) -> List[CitationResponse]:
    """
    Get papers that cite this paper.

    Returns citations stored in the PaperCitation table.
    """
    # Verify paper exists
    literature = db.query(Literature).filter(Literature.id == paper_id).first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    # Get citations
    citations = (
        db.query(PaperCitation)
        .filter(PaperCitation.cited_paper_id == paper_id)
        .order_by(PaperCitation.is_influential.desc())
        .offset(offset)
        .limit(limit)
        .all()
    )
    
    results = []
    for citation in citations:
        # Try to get citing paper details from Literature table
        citing_paper = citation.citing_paper if citation.citing_paper else None
        
        results.append(
            CitationResponse(
                paper_id=citing_paper.id if citing_paper else None,
                doi=citing_paper.doi if citing_paper else None,
                title=citing_paper.title if citing_paper else "Unknown",
                is_influential=citation.is_influential,
                citation_context=citation.citation_context,
            )
        )
    
    return results


@router.get("/{paper_id}/references", response_model=List[CitationResponse])
async def get_paper_references(
    paper_id: UUID,
    limit: int = Query(50, ge=1, le=200),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_database_session),
) -> List[CitationResponse]:
    """
    Get papers that this paper cites (references).

    Returns references stored in the PaperCitation table.
    """
    # Verify paper exists
    literature = db.query(Literature).filter(Literature.id == paper_id).first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    # Get references
    references = (
        db.query(PaperCitation)
        .filter(PaperCitation.citing_paper_id == paper_id)
        .order_by(PaperCitation.is_influential.desc())
        .offset(offset)
        .limit(limit)
        .all()
    )
    
    results = []
    for ref in references:
        # Try to get cited paper details from Literature table or use stored metadata
        cited_paper = ref.cited_paper if ref.cited_paper else None
        
        results.append(
            CitationResponse(
                paper_id=cited_paper.id if cited_paper else None,
                doi=ref.cited_doi or (cited_paper.doi if cited_paper else None),
                title=ref.cited_title or (cited_paper.title if cited_paper else "Unknown"),
                is_influential=ref.is_influential,
                citation_context=None,  # References don't have context
            )
        )
    
    return results


@router.post("/{paper_id}/enrich", response_model=EnrichPaperResponse)
async def enrich_paper(
    paper_id: UUID,
    request: EnrichPaperRequest,
    db: Session = Depends(get_database_session),
) -> EnrichPaperResponse:
    """
    Enrich paper with Semantic Scholar and OpenAlex metadata.

    Fetches additional metadata and optionally builds citation graph.
    """
    # Get paper
    literature = db.query(Literature).filter(Literature.id == paper_id).first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    citations_added = 0
    references_added = 0
    
    # Try Semantic Scholar enrichment
    s2_repo = SemanticScholarRepository()
    s2_paper_id = None
    
    try:
        # Search S2 by DOI or PMID
        if literature.doi:
            s2_paper_id = f"DOI:{literature.doi}"
        elif literature.pmid:
            s2_paper_id = f"PMID:{literature.pmid}"
        
        if s2_paper_id:
            s2_metadata = s2_repo.fetch_metadata(s2_paper_id)
            if s2_metadata:
                literature.semantic_scholar_id = s2_metadata.paper_id
                # Update citation counts if available in S2 metadata
                logger.info("[PAPERS_API] Enriched with Semantic Scholar: %s", s2_metadata.paper_id)
        
        # Fetch citations if requested
        if request.fetch_citations and s2_paper_id:
            try:
                citations = s2_repo.get_citations(s2_paper_id, limit=100)
                for cite in citations:
                    citing_paper_s2_id = cite.get("paperId")
                    if citing_paper_s2_id:
                        # Check if we already have this citation
                        existing = db.query(PaperCitation).filter(
                            PaperCitation.citing_paper_id == paper_id,
                        ).first()
                        
                        if not existing:
                            citation = PaperCitation(
                                citing_paper_id=paper_id,  # This will need to be matched later
                                cited_paper_id=paper_id,
                                cited_title=cite.get("title", ""),
                                is_influential=cite.get("isInfluential", False),
                            )
                            db.add(citation)
                            citations_added += 1
            except Exception as e:
                logger.warning("[PAPERS_API] Failed to fetch S2 citations: %r", e)
        
        # Fetch references if requested
        if request.fetch_references and s2_paper_id:
            try:
                references = s2_repo.get_references(s2_paper_id, limit=100)
                for ref in references:
                    cited_paper_s2_id = ref.get("paperId")
                    if cited_paper_s2_id:
                        # Check if we already have this reference
                        existing = db.query(PaperCitation).filter(
                            PaperCitation.citing_paper_id == paper_id,
                        ).first()
                        
                        if not existing:
                            reference = PaperCitation(
                                citing_paper_id=paper_id,
                                cited_paper_id=None,  # Don't know if it's in our system
                                cited_title=ref.get("title", ""),
                                is_influential=ref.get("isInfluential", False),
                            )
                            db.add(reference)
                            references_added += 1
            except Exception as e:
                logger.warning("[PAPERS_API] Failed to fetch S2 references: %r", e)
    
    except Exception as e:
        logger.warning("[PAPERS_API] Semantic Scholar enrichment failed: %r", e)
    
    # Try OpenAlex enrichment
    try:
        oa_repo = OpenAlexRepository()
        if literature.doi:
            oa_work_id = f"https://doi.org/{literature.doi}"
            oa_metadata = oa_repo.get_work(oa_work_id)
            if oa_metadata:
                literature.openalex_id = oa_metadata.paper_id
                logger.info("[PAPERS_API] Enriched with OpenAlex: %s", oa_metadata.paper_id)
    except Exception as e:
        logger.warning("[PAPERS_API] OpenAlex enrichment failed: %r", e)
    
    db.commit()
    
    return EnrichPaperResponse(
        enriched=True,
        citations_added=citations_added,
        references_added=references_added,
        semantic_scholar_id=literature.semantic_scholar_id,
        openalex_id=literature.openalex_id,
        citation_count=literature.citation_count,
    )

