"""Data Catalog API endpoints."""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, status
from sqlalchemy.orm import Session

from amprenta_rag.api import schemas
from amprenta_rag.api.dependencies import get_db, get_current_user
from amprenta_rag.database.models import User
from amprenta_rag.services import catalog as service

router = APIRouter(prefix="/catalog", tags=["Data Catalog"])


# ============================================================================
# CATALOG BROWSING
# ============================================================================

@router.get("/entries", response_model=List[schemas.CatalogEntryResponse])
def list_catalog_entries(
    category: Optional[str] = Query(None, description="Filter by category"),
    search: Optional[str] = Query(None, description="Search entity names"),
    db: Session = Depends(get_db),
):
    """List all catalog entries with optional filtering."""
    entries = service.get_catalog_entries(db, category=category, search=search)
    return [schemas.CatalogEntryResponse.model_validate(entry) for entry in entries]


@router.get("/entries/{entity_type}", response_model=schemas.CatalogEntryDetailResponse)
def get_catalog_entry(
    entity_type: str,
    db: Session = Depends(get_db),
):
    """Get single catalog entry with all columns."""
    entry = service.get_catalog_entry(db, entity_type)
    if not entry:
        raise HTTPException(status_code=404, detail=f"Entity type '{entity_type}' not found")
    return schemas.CatalogEntryDetailResponse.model_validate(entry)


@router.put("/entries/{entity_type}", response_model=schemas.CatalogEntryResponse)
def update_catalog_entry(
    entity_type: str,
    data: schemas.CatalogEntryUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Update catalog entry metadata (admin only)."""
    entry = service.update_catalog_entry(
        db, entity_type, 
        description=data.description, 
        category=data.category
    )
    if not entry:
        raise HTTPException(status_code=404, detail=f"Entity type '{entity_type}' not found")
    return schemas.CatalogEntryResponse.model_validate(entry)


@router.post("/refresh", response_model=dict)
def refresh_catalog(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Trigger catalog refresh (discover new models). Admin only."""
    count = service.refresh_catalog(db)
    lineage_count = service.detect_lineage_from_fks(db)
    return {"entries_updated": count, "lineage_edges_created": lineage_count}


# ============================================================================
# COLUMN SEARCH
# ============================================================================

@router.get("/columns/search", response_model=List[schemas.ColumnSearchResult])
def search_columns(
    q: str = Query(..., min_length=2, description="Search query"),
    limit: int = Query(50, le=200),
    db: Session = Depends(get_db),
):
    """Search columns by name or description across all entities."""
    columns = service.search_columns(db, query=q, limit=limit)
    results = []
    for col in columns:
        results.append(schemas.ColumnSearchResult(
            entity_type=col.catalog_entry.entity_type,
            column_name=col.column_name,
            display_name=col.display_name,
            data_type=col.data_type,
            description=col.description,
        ))
    return results


@router.get("/columns/{entity_type}/{column_name}", response_model=schemas.ColumnMetadataResponse)
def get_column_metadata(
    entity_type: str,
    column_name: str,
    db: Session = Depends(get_db),
):
    """Get detailed column metadata."""
    col = service.get_column_metadata(db, entity_type, column_name)
    if not col:
        raise HTTPException(status_code=404, detail="Column not found")
    return schemas.ColumnMetadataResponse.model_validate(col)


@router.put("/columns/{column_id}", response_model=schemas.ColumnMetadataResponse)
def update_column_metadata(
    column_id: UUID,
    data: schemas.ColumnMetadataUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Update column metadata (display name, description, glossary link)."""
    col = service.update_column_metadata(
        db, column_id,
        display_name=data.display_name,
        description=data.description,
        glossary_term_id=data.glossary_term_id
    )
    if not col:
        raise HTTPException(status_code=404, detail="Column not found")
    return schemas.ColumnMetadataResponse.model_validate(col)


# ============================================================================
# GLOSSARY
# ============================================================================

@router.get("/glossary", response_model=List[schemas.GlossaryTermResponse])
def list_glossary_terms(
    category: Optional[str] = Query(None),
    search: Optional[str] = Query(None),
    db: Session = Depends(get_db),
):
    """List glossary terms with optional filtering."""
    terms = service.get_glossary_terms(db, category=category, search=search)
    return [schemas.GlossaryTermResponse.model_validate(term) for term in terms]


@router.post("/glossary", response_model=schemas.GlossaryTermResponse, status_code=201)
def create_glossary_term(
    data: schemas.GlossaryTermCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Create new glossary term."""
    try:
        term = service.create_glossary_term(
            db,
            term=data.term,
            definition=data.definition,
            category=data.category,
            synonyms=data.synonyms or [],
            created_by_id=current_user.id
        )
        return schemas.GlossaryTermResponse.model_validate(term)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/glossary/{term_id}", response_model=schemas.GlossaryTermResponse)
def get_glossary_term(
    term_id: UUID,
    db: Session = Depends(get_db),
):
    """Get single glossary term."""
    term = service.get_glossary_term(db, term_id)
    if not term:
        raise HTTPException(status_code=404, detail="Term not found")
    return schemas.GlossaryTermResponse.model_validate(term)


@router.put("/glossary/{term_id}", response_model=schemas.GlossaryTermResponse)
def update_glossary_term(
    term_id: UUID,
    data: schemas.GlossaryTermUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Update glossary term."""
    term = service.update_glossary_term(
        db, term_id,
        definition=data.definition,
        category=data.category,
        synonyms=data.synonyms
    )
    if not term:
        raise HTTPException(status_code=404, detail="Term not found")
    return schemas.GlossaryTermResponse.model_validate(term)


@router.delete("/glossary/{term_id}", status_code=204)
def delete_glossary_term(
    term_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Delete glossary term (unlinks from columns first)."""
    success = service.delete_glossary_term(db, term_id)
    if not success:
        raise HTTPException(status_code=404, detail="Term not found")


# ============================================================================
# LINEAGE
# ============================================================================

@router.get("/lineage/{entity_type}/{entity_id}", response_model=schemas.LineageGraphResponse)
def get_lineage_graph(
    entity_type: str,
    entity_id: UUID,
    depth: int = Query(3, ge=1, le=5),
    direction: str = Query("both", pattern="^(upstream|downstream|both)$"),
    db: Session = Depends(get_db),
):
    """
    Get lineage graph for an entity.
    
    Returns Cytoscape.js compatible JSON with nodes and edges.
    Max 500 nodes returned for performance.
    """
    graph = service.get_lineage_graph(
        db, entity_type, entity_id, 
        depth=depth, 
        direction=direction
    )
    return schemas.LineageGraphResponse.model_validate(graph)


@router.post("/lineage", response_model=dict, status_code=201)
def add_lineage_edge(
    data: schemas.LineageEdgeCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Add manual lineage edge between entities."""
    try:
        edge = service.add_lineage_edge(
            db,
            source_type=data.source_type,
            source_id=data.source_id,
            target_type=data.target_type,
            target_id=data.target_id,
            relationship_type=data.relationship_type,
            transformation=data.transformation
        )
        return {"id": str(edge.id), "created": True}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))