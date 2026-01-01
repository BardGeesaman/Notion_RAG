"""
API router for scientific papers ingestion and search.

Provides endpoints for searching PubMed/PMC, ingesting papers with
deduplication, and retrieving paper content.
"""

from __future__ import annotations

import asyncio
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, File, HTTPException, Query, UploadFile
from pydantic import BaseModel, ConfigDict, Field
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select

from amprenta_rag.api.async_dependencies import get_async_database_session
from amprenta_rag.database.models import Literature
from amprenta_rag.ingestion.papers import PubMedRepository
from amprenta_rag.ingestion.papers.semantic_scholar import SemanticScholarRepository
from amprenta_rag.ingestion.papers.openalex import OpenAlexRepository
from amprenta_rag.models.content import PaperCitation, PublicationExtraction, SupplementaryFile
from amprenta_rag.extraction.supplementary_parser import parse_supplementary_file
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

router = APIRouter()


# Sync helper functions for external API calls
def _sync_search_papers(query: str, source: str, limit: int, offset: int):
    """Sync helper for paper search via external APIs."""
    if source.lower() == "pubmed":
        repo = PubMedRepository()
        pmids = repo.search_papers(query, max_results=limit + offset)
        paginated_pmids = pmids[offset : offset + limit]
        
        results = []
        for pmid in paginated_pmids:
            metadata = repo.fetch_metadata(pmid)
            if metadata:
                results.append({
                    "paper_id": metadata.paper_id,
                    "title": metadata.title,
                    "abstract": metadata.abstract,
                    "authors": metadata.authors,
                    "journal": metadata.journal,
                    "year": metadata.year,
                    "doi": metadata.doi,
                    "pmid": metadata.pmid,
                    "url": metadata.url,
                })
        
        return {"results": results, "total": len(pmids)}
    else:
        raise ValueError(f"Unsupported source: {source}")


def _sync_fetch_paper_metadata(pmid: str = None, doi: str = None):
    """Sync helper for fetching paper metadata from PubMed."""
    repo = PubMedRepository()
    
    if pmid:
        return repo.fetch_metadata(pmid)
    elif doi:
        # Search by DOI to get PMID, then fetch
        pmids = repo.search_papers(f'"{doi}"[DOI]', max_results=1)
        if not pmids:
            return None
        return repo.fetch_metadata(pmids[0])
    else:
        return None


def _sync_enrich_with_semantic_scholar(paper_id: str, fetch_citations: bool, fetch_references: bool):
    """Sync helper for Semantic Scholar enrichment."""
    s2_repo = SemanticScholarRepository()
    
    result = {
        "metadata": None,
        "citations": [],
        "references": []
    }
    
    try:
        # Fetch metadata
        metadata = s2_repo.fetch_metadata(paper_id)
        result["metadata"] = metadata
        
        # Fetch citations if requested
        if fetch_citations:
            citations = s2_repo.get_citations(paper_id, limit=100)
            result["citations"] = citations
        
        # Fetch references if requested  
        if fetch_references:
            references = s2_repo.get_references(paper_id, limit=100)
            result["references"] = references
            
    except Exception as e:
        logger.warning("Semantic Scholar API call failed: %r", e)
        
    return result


def _sync_enrich_with_openalex(work_id: str):
    """Sync helper for OpenAlex enrichment."""
    try:
        oa_repo = OpenAlexRepository()
        return oa_repo.get_work(work_id)
    except Exception as e:
        logger.warning("OpenAlex API call failed: %r", e)
        return None


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


class ExperimentResponse(BaseModel):
    """Extracted experiment response."""

    id: UUID
    experiment_type: str
    cell_line: Optional[str] = None
    treatment: Optional[str] = None
    concentration: Optional[str] = None
    timepoint: Optional[str] = None
    replicate_count: Optional[int] = None
    measured_entities: List[str] = Field(default_factory=list)
    key_findings: Optional[str] = None
    extraction_confidence: Optional[int] = None

    model_config = ConfigDict(from_attributes=True)


class ExtractExperimentsResponse(BaseModel):
    """Response for experiment extraction from PDF."""

    experiments_extracted: int
    extraction_confidence: float


class SupplementaryUploadResponse(BaseModel):
    """Response for supplementary file upload."""

    file_id: UUID
    filename: str
    schema_type: str
    row_count: int
    tables_extracted: int


class LinkDatasetRequest(BaseModel):
    """Request to link supplementary file to dataset."""

    supplementary_file_id: UUID
    dataset_id: UUID


class LinkDatasetResponse(BaseModel):
    """Response for dataset linking."""

    linked: bool
    supplementary_file_id: UUID
    dataset_id: UUID


@router.post("/search", response_model=PaperSearchResponse)
async def search_papers(
    request: PaperSearchRequest,
    db: AsyncSession = Depends(get_async_database_session),
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
        # Search papers using async thread pool
        search_result = await asyncio.to_thread(
            _sync_search_papers,
            request.query,
            request.source,
            request.limit,
            request.offset
        )
        
        # Convert to response objects
        results = [
            PaperSearchResult(**result_data) 
            for result_data in search_result["results"]
        ]

        return PaperSearchResponse(
            results=results,
            total=search_result["total"],
            offset=request.offset,
            limit=request.limit,
        )

    except Exception as e:
        logger.error("[PAPERS_API] Search failed: %r", e)
        raise HTTPException(status_code=500, detail=f"Search failed: {str(e)}")


@router.post("/ingest", response_model=PaperIngestResponse, status_code=201)
async def ingest_paper(
    request: PaperIngestRequest,
    db: AsyncSession = Depends(get_async_database_session),
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
        result = await db.execute(
            select(Literature).filter(Literature.pmid == request.pmid)
        )
        existing = result.scalars().first()
    if not existing and request.doi:
        result = await db.execute(
            select(Literature).filter(Literature.doi == request.doi)
        )
        existing = result.scalars().first()

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

    # Fetch metadata from PubMed using async thread pool
    try:
        metadata = await asyncio.to_thread(
            _sync_fetch_paper_metadata,
            request.pmid,
            request.doi
        )

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
        await db.flush()  # Get literature.id

        logger.info(
            "[PAPERS_API] Created Literature record: PMID=%s, ID=%s",
            metadata.pmid,
            literature.id,
        )

        # NOTE: Full text embedding tracked in ROADMAP
        # For now, just commit metadata
        chunks_created = 0
        if request.fetch_fulltext:
            logger.info(
                "[PAPERS_API] Full text fetching not yet implemented for PMID=%s",
                metadata.pmid,
            )

        await db.commit()

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
    db: AsyncSession = Depends(get_async_database_session),
) -> PaperDetailResponse:
    """
    Get paper details by ID.

    Returns paper metadata and sections if available.
    """
    result = await db.execute(
        select(Literature).filter(Literature.id == paper_id)
    )
    literature = result.scalars().first()

    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")

    # NOTE: RAG chunk section retrieval tracked in ROADMAP
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
    db: AsyncSession = Depends(get_async_database_session),
) -> List[CitationResponse]:
    """
    Get papers that cite this paper.

    Returns citations stored in the PaperCitation table.
    """
    # Verify paper exists
    result = await db.execute(
        select(Literature).filter(Literature.id == paper_id)
    )
    literature = result.scalars().first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    # Get citations
    citations_result = await db.execute(
        select(PaperCitation)
        .filter(PaperCitation.cited_paper_id == paper_id)
        .order_by(PaperCitation.is_influential.desc())
        .offset(offset)
        .limit(limit)
    )
    citations = citations_result.scalars().all()
    
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
    db: AsyncSession = Depends(get_async_database_session),
) -> List[CitationResponse]:
    """
    Get papers that this paper cites (references).

    Returns references stored in the PaperCitation table.
    """
    # Verify paper exists
    result = await db.execute(
        select(Literature).filter(Literature.id == paper_id)
    )
    literature = result.scalars().first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    # Get references
    references_result = await db.execute(
        select(PaperCitation)
        .filter(PaperCitation.citing_paper_id == paper_id)
        .order_by(PaperCitation.is_influential.desc())
        .offset(offset)
        .limit(limit)
    )
    references = references_result.scalars().all()
    
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
    db: AsyncSession = Depends(get_async_database_session),
) -> EnrichPaperResponse:
    """
    Enrich paper with Semantic Scholar and OpenAlex metadata.

    Fetches additional metadata and optionally builds citation graph.
    """
    # Get paper
    result = await db.execute(
        select(Literature).filter(Literature.id == paper_id)
    )
    literature = result.scalars().first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    citations_added = 0
    references_added = 0
    
    # Try Semantic Scholar enrichment using async thread pool
    s2_paper_id = None
    
    try:
        # Search S2 by DOI or PMID
        if literature.doi:
            s2_paper_id = f"DOI:{literature.doi}"
        elif literature.pmid:
            s2_paper_id = f"PMID:{literature.pmid}"
        
        if s2_paper_id:
            s2_result = await asyncio.to_thread(
                _sync_enrich_with_semantic_scholar,
                s2_paper_id,
                request.fetch_citations,
                request.fetch_references
            )
            
            # Process metadata
            if s2_result["metadata"]:
                literature.semantic_scholar_id = s2_result["metadata"].paper_id
                logger.info("[PAPERS_API] Enriched with Semantic Scholar: %s", s2_result["metadata"].paper_id)
            
            # Process citations (papers that CITE our paper)
            for cite in s2_result["citations"]:
                citing_paper_s2_id = cite.get("paperId")
                if citing_paper_s2_id:
                    cite_title = cite.get("title", "")
                    # Check if we already have this citation (avoid duplicates by title)
                    existing_result = await db.execute(
                        select(PaperCitation).filter(
                            PaperCitation.cited_paper_id == paper_id,
                            PaperCitation.cited_title == cite_title,
                        )
                    )
                    existing = existing_result.scalars().first()
                    
                    if not existing and cite_title:
                        # External paper cites our paper
                        external_ids = cite.get("externalIds", {}) or {}
                        citation = PaperCitation(
                            citing_paper_id=None,  # External paper not in system yet
                            cited_paper_id=paper_id,  # Our paper being cited
                            cited_doi=external_ids.get("DOI"),
                            cited_title=cite_title,
                            citation_context="",  # S2 citations API doesn't include context
                            is_influential=cite.get("isInfluential", False),
                        )
                        db.add(citation)
                        citations_added += 1
            
            # Process references (papers that our paper CITES)
            for ref in s2_result["references"]:
                cited_paper_s2_id = ref.get("paperId")
                if cited_paper_s2_id:
                    ref_title = ref.get("title", "")
                    # Check if we already have this reference (avoid duplicates by title)
                    existing_result = await db.execute(
                        select(PaperCitation).filter(
                            PaperCitation.citing_paper_id == paper_id,
                            PaperCitation.cited_title == ref_title,
                        )
                    )
                    existing = existing_result.scalars().first()
                    
                    if not existing and ref_title:
                        # Our paper cites an external paper
                        external_ids = ref.get("externalIds", {}) or {}
                        reference = PaperCitation(
                            citing_paper_id=paper_id,  # Our paper doing the citing
                            cited_paper_id=None,  # External paper not in system yet
                            cited_doi=external_ids.get("DOI"),
                            cited_title=ref_title,
                            citation_context=None,  # References don't have context
                            is_influential=ref.get("isInfluential", False),
                        )
                        db.add(reference)
                        references_added += 1
    
    except Exception as e:
        logger.warning("[PAPERS_API] Semantic Scholar enrichment failed: %r", e)
    
    # Try OpenAlex enrichment using async thread pool
    try:
        if literature.doi:
            oa_work_id = f"https://doi.org/{literature.doi}"
            oa_metadata = await asyncio.to_thread(
                _sync_enrich_with_openalex,
                oa_work_id
            )
            if oa_metadata:
                literature.openalex_id = oa_metadata.paper_id
                logger.info("[PAPERS_API] Enriched with OpenAlex: %s", oa_metadata.paper_id)
    except Exception as e:
        logger.warning("[PAPERS_API] OpenAlex enrichment failed: %r", e)
    
    await db.commit()
    
    return EnrichPaperResponse(
        enriched=True,
        citations_added=citations_added,
        references_added=references_added,
        semantic_scholar_id=literature.semantic_scholar_id,
        openalex_id=literature.openalex_id,
        citation_count=literature.citation_count,
    )


@router.post("/{paper_id}/extract", response_model=ExtractExperimentsResponse)
async def extract_experiments_from_paper(
    paper_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
) -> ExtractExperimentsResponse:
    """
    Extract structured experiment data from paper PDF.

    NOTE: Requires PDF to be available. Currently placeholder - full implementation
    requires PDF storage integration.
    """
    result = await db.execute(
        select(Literature).filter(Literature.id == paper_id)
    )
    literature = result.scalars().first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    # Placeholder: Would need PDF bytes from storage
    # For now, return empty extraction
    logger.warning("[PAPERS_API] PDF extraction not yet implemented - requires PDF storage")
    
    return ExtractExperimentsResponse(
        experiments_extracted=0,
        extraction_confidence=0.0,
    )


@router.post("/{paper_id}/supplementary", response_model=SupplementaryUploadResponse)
async def upload_supplementary_file(
    paper_id: UUID,
    file: UploadFile = File(...),
    db: AsyncSession = Depends(get_async_database_session),
) -> SupplementaryUploadResponse:
    """
    Upload and parse supplementary file for a paper.

    Accepts Excel or CSV files, auto-detects schema, and stores parsed data.
    """
    import tempfile
    from pathlib import Path
    
    # Verify paper exists
    result = await db.execute(
        select(Literature).filter(Literature.id == paper_id)
    )
    literature = result.scalars().first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    # Save uploaded file temporarily
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=Path(file.filename or "file").suffix) as tmp_file:
            content = await file.read()
            tmp_file.write(content)
            tmp_path = tmp_file.name
        
        # Parse file
        result = parse_supplementary_file(tmp_path)
        
        # Extract first table's schema info
        schema_type = "generic"
        row_count = 0
        columns = []
        data = {}
        
        if result.tables:
            first_table = result.tables[0]
            schema_type = first_table.detected_schema.schema_type
            row_count = first_table.row_count
            columns = first_table.columns
            data = {"tables": [t.dict() for t in result.tables]}
        
        # Store in database
        supp_file = SupplementaryFile(
            literature_id=paper_id,
            filename=file.filename or "unknown",
            file_type=result.file_type,
            schema_type=schema_type,
            row_count=row_count,
            columns=columns,
            data=data,
        )
        db.add(supp_file)
        await db.commit()
        db.refresh(supp_file)
        
        logger.info(
            "[PAPERS_API] Uploaded supplementary file %s for paper %s",
            file.filename,
            paper_id,
        )
        
        # Cleanup temp file
        Path(tmp_path).unlink()
        
        return SupplementaryUploadResponse(
            file_id=supp_file.id,
            filename=supp_file.filename,
            schema_type=schema_type,
            row_count=row_count,
            tables_extracted=len(result.tables),
        )
    
    except Exception as e:
        logger.error("[PAPERS_API] Supplementary upload failed: %r", e)
        raise HTTPException(status_code=500, detail=f"Upload failed: {str(e)}")


@router.get("/{paper_id}/experiments", response_model=List[ExperimentResponse])
async def get_paper_experiments(
    paper_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
) -> List[ExperimentResponse]:
    """
    Get extracted experiments for a paper.

    Returns all PublicationExtraction records for the specified paper.
    """
    result = await db.execute(
        select(Literature).filter(Literature.id == paper_id)
    )
    literature = result.scalars().first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    experiments_result = await db.execute(
        select(PublicationExtraction).filter(PublicationExtraction.literature_id == paper_id)
    )
    experiments = experiments_result.scalars().all()
    
    return [ExperimentResponse.model_validate(exp) for exp in experiments]


@router.post("/{paper_id}/link-dataset", response_model=LinkDatasetResponse)
async def link_supplementary_to_dataset(
    paper_id: UUID,
    request: LinkDatasetRequest,
    db: AsyncSession = Depends(get_async_database_session),
) -> LinkDatasetResponse:
    """
    Link a supplementary file to an existing Dataset.

    Updates the SupplementaryFile record to reference a Dataset ID.
    """
    # Verify paper exists
    result = await db.execute(
        select(Literature).filter(Literature.id == paper_id)
    )
    literature = result.scalars().first()
    if not literature:
        raise HTTPException(status_code=404, detail="Paper not found")
    
    # Get supplementary file
    supp_result = await db.execute(
        select(SupplementaryFile).filter(
            SupplementaryFile.id == request.supplementary_file_id,
            SupplementaryFile.literature_id == paper_id,
        )
    )
    supp_file = supp_result.scalars().first()
    
    if not supp_file:
        raise HTTPException(status_code=404, detail="Supplementary file not found for this paper")
    
    # Update linked dataset
    supp_file.linked_dataset_id = request.dataset_id
    await db.commit()
    
    logger.info(
        "[PAPERS_API] Linked supplementary file %s to dataset %s",
        request.supplementary_file_id,
        request.dataset_id,
    )
    
    return LinkDatasetResponse(
        linked=True,
        supplementary_file_id=request.supplementary_file_id,
        dataset_id=request.dataset_id,
    )

