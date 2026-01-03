"""Data Catalog service layer."""

from __future__ import annotations

import logging
from datetime import datetime
from typing import List, Optional
from uuid import UUID

from sqlalchemy import and_, or_, text, inspect, func
from sqlalchemy.orm import Session

from amprenta_rag.database.models import CatalogEntry, ColumnMetadata, GlossaryTerm, DataLineageEdge

logger = logging.getLogger(__name__)


# ============================================================================
# CATALOG ENTRY OPERATIONS
# ============================================================================

def get_catalog_entries(
    db: Session,
    category: Optional[str] = None,
    search: Optional[str] = None,
    limit: int = 100,
) -> List[CatalogEntry]:
    """
    List catalog entries with optional filters.
    
    Args:
        db: Database session
        category: Filter by category (Core, Chemistry, Omics, Admin)
        search: Search in display_name and description
        limit: Maximum results to return
    
    Returns:
        List of catalog entries
    """
    query = db.query(CatalogEntry).filter(CatalogEntry.is_active.is_(True))
    
    if category:
        query = query.filter(CatalogEntry.category == category)
    
    if search:
        search_pattern = f"%{search}%"
        query = query.filter(
            or_(
                CatalogEntry.display_name.ilike(search_pattern),
                CatalogEntry.description.ilike(search_pattern),
                CatalogEntry.entity_type.ilike(search_pattern),
            )
        )
    
    return query.order_by(CatalogEntry.display_name).limit(limit).all()


def get_catalog_entry(db: Session, entity_type: str) -> Optional[CatalogEntry]:
    """
    Get catalog entry by entity type with columns.
    
    Args:
        db: Database session
        entity_type: Entity type identifier
    
    Returns:
        Catalog entry with columns loaded, or None if not found
    """
    return db.query(CatalogEntry).filter(
        CatalogEntry.entity_type == entity_type,
        CatalogEntry.is_active.is_(True)
    ).first()


def update_catalog_entry(
    db: Session,
    entity_type: str,
    description: Optional[str] = None,
    category: Optional[str] = None,
    row_count: Optional[int] = None,
) -> Optional[CatalogEntry]:
    """
    Update catalog entry metadata.
    
    Args:
        db: Database session
        entity_type: Entity type to update
        description: New description
        category: New category
        row_count: Updated row count
    
    Returns:
        Updated catalog entry or None if not found
    """
    entry = get_catalog_entry(db, entity_type)
    if not entry:
        return None
    
    if description is not None:
        entry.description = description
    if category is not None:
        entry.category = category
    if row_count is not None:
        entry.row_count = row_count
        entry.last_refreshed = datetime.utcnow()
    
    db.commit()
    db.refresh(entry)
    logger.info(f"Updated catalog entry: {entity_type}")
    return entry


# ============================================================================
# COLUMN OPERATIONS
# ============================================================================

def search_columns(
    db: Session,
    query: str,
    limit: int = 50,
) -> List[ColumnMetadata]:
    """
    Search columns by name or description.
    
    Args:
        db: Database session
        query: Search term
        limit: Maximum results
    
    Returns:
        List of matching columns with catalog entry loaded
    """
    search_pattern = f"%{query}%"
    
    return db.query(ColumnMetadata).join(CatalogEntry).filter(
        and_(
            CatalogEntry.is_active.is_(True),
            or_(
                ColumnMetadata.column_name.ilike(search_pattern),
                ColumnMetadata.display_name.ilike(search_pattern),
                ColumnMetadata.description.ilike(search_pattern),
            )
        )
    ).limit(limit).all()


# ============================================================================
# GLOSSARY OPERATIONS
# ============================================================================

def create_glossary_term(
    db: Session,
    term: str,
    definition: str,
    category: Optional[str] = None,
    synonyms: Optional[List[str]] = None,
    related_terms: Optional[List[str]] = None,
    source: Optional[str] = None,
    created_by_id: Optional[UUID] = None,
) -> GlossaryTerm:
    """
    Create a new glossary term.
    
    Args:
        db: Database session
        term: Term name (must be unique)
        definition: Term definition
        category: Optional category
        synonyms: List of synonyms
        related_terms: List of related term names
        source: Source reference
        created_by_id: User who created the term
    
    Returns:
        Created glossary term
    
    Raises:
        ValueError: If term already exists
    """
    # Check if term already exists
    existing = db.query(GlossaryTerm).filter(GlossaryTerm.term == term).first()
    if existing:
        raise ValueError(f"Glossary term '{term}' already exists")
    
    glossary_term = GlossaryTerm(
        term=term,
        definition=definition,
        category=category,
        synonyms=synonyms,
        related_terms=related_terms,
        source=source,
        created_by_id=created_by_id,
    )
    
    db.add(glossary_term)
    db.commit()
    db.refresh(glossary_term)
    logger.info(f"Created glossary term: {term}")
    return glossary_term


def get_glossary_terms(
    db: Session,
    category: Optional[str] = None,
    search: Optional[str] = None,
    limit: int = 100,
) -> List[GlossaryTerm]:
    """
    List glossary terms with optional filters.
    
    Args:
        db: Database session
        category: Filter by category
        search: Search in term, definition, synonyms
        limit: Maximum results
    
    Returns:
        List of glossary terms
    """
    query = db.query(GlossaryTerm)
    
    if category:
        query = query.filter(GlossaryTerm.category == category)
    
    if search:
        search_pattern = f"%{search}%"
        query = query.filter(
            or_(
                GlossaryTerm.term.ilike(search_pattern),
                GlossaryTerm.definition.ilike(search_pattern),
                GlossaryTerm.synonyms.op("@>")(f'["{search}"]'),  # JSONB contains
            )
        )
    
    return query.order_by(GlossaryTerm.term).limit(limit).all()


def update_glossary_term(
    db: Session,
    term_id: UUID,
    term: Optional[str] = None,
    definition: Optional[str] = None,
    category: Optional[str] = None,
    synonyms: Optional[List[str]] = None,
    related_terms: Optional[List[str]] = None,
    source: Optional[str] = None,
) -> Optional[GlossaryTerm]:
    """
    Update glossary term.
    
    Args:
        db: Database session
        term_id: Term ID to update
        term: New term name
        definition: New definition
        category: New category
        synonyms: New synonyms list
        related_terms: New related terms list
        source: New source
    
    Returns:
        Updated term or None if not found
    """
    glossary_term = db.query(GlossaryTerm).filter(GlossaryTerm.id == term_id).first()
    if not glossary_term:
        return None
    
    if term is not None:
        # Check uniqueness if changing term name
        if term != glossary_term.term:
            existing = db.query(GlossaryTerm).filter(GlossaryTerm.term == term).first()
            if existing:
                raise ValueError(f"Term '{term}' already exists")
        glossary_term.term = term
    
    if definition is not None:
        glossary_term.definition = definition
    if category is not None:
        glossary_term.category = category
    if synonyms is not None:
        glossary_term.synonyms = synonyms
    if related_terms is not None:
        glossary_term.related_terms = related_terms
    if source is not None:
        glossary_term.source = source
    
    db.commit()
    db.refresh(glossary_term)
    logger.info(f"Updated glossary term: {glossary_term.term}")
    return glossary_term


def delete_glossary_term(db: Session, term_id: UUID) -> bool:
    """
    Delete glossary term (unlinks columns first).
    
    Args:
        db: Database session
        term_id: Term ID to delete
    
    Returns:
        True if deleted, False if not found
    """
    glossary_term = db.query(GlossaryTerm).filter(GlossaryTerm.id == term_id).first()
    if not glossary_term:
        return False
    
    # Unlink columns first
    db.query(ColumnMetadata).filter(
        ColumnMetadata.glossary_term_id == term_id
    ).update({"glossary_term_id": None})
    
    # Delete term
    db.delete(glossary_term)
    db.commit()
    logger.info(f"Deleted glossary term: {glossary_term.term}")
    return True


# ============================================================================
# DATA LINEAGE OPERATIONS
# ============================================================================

def add_lineage_edge(
    db: Session,
    source_type: str,
    source_id: UUID,
    target_type: str,
    target_id: UUID,
    relationship_type: str,
    transformation: Optional[str] = None,
    edge_metadata: Optional[dict] = None,
) -> DataLineageEdge:
    """
    Add data lineage edge.
    
    Args:
        db: Database session
        source_type: Source entity type
        source_id: Source entity ID
        target_type: Target entity type
        target_id: Target entity ID
        relationship_type: Type of relationship (derived_from, input_to, generated_by)
        transformation: Optional transformation description
        edge_metadata: Optional metadata dict
    
    Returns:
        Created lineage edge
    
    Raises:
        ValueError: If edge already exists (due to unique constraint)
    """
    # Check if edge already exists
    existing = db.query(DataLineageEdge).filter(
        DataLineageEdge.source_type == source_type,
        DataLineageEdge.source_id == source_id,
        DataLineageEdge.target_type == target_type,
        DataLineageEdge.target_id == target_id,
    ).first()
    
    if existing:
        raise ValueError(f"Lineage edge already exists: {source_type}:{source_id} -> {target_type}:{target_id}")
    
    edge = DataLineageEdge(
        source_type=source_type,
        source_id=source_id,
        target_type=target_type,
        target_id=target_id,
        relationship_type=relationship_type,
        transformation=transformation,
        edge_metadata=edge_metadata,
    )
    
    db.add(edge)
    db.commit()
    db.refresh(edge)
    logger.info(f"Added lineage edge: {source_type}:{source_id} -> {target_type}:{target_id}")
    return edge


def get_entity_lineage(
    db: Session,
    entity_type: str,
    entity_id: UUID,
    max_nodes: int = 500,
) -> dict:
    """
    Get immediate parents and children for an entity.
    
    P1 FIX: Limited to 500 nodes max for query safety.
    Query timeout will be implemented in Batch 2.
    
    Args:
        db: Database session
        entity_type: Entity type
        entity_id: Entity ID
        max_nodes: Maximum nodes to return (safety limit)
    
    Returns:
        Dict with 'parents' and 'children' lists
    """
    # Get parents (entities that flow into this one)
    parents = db.query(DataLineageEdge).filter(
        DataLineageEdge.target_type == entity_type,
        DataLineageEdge.target_id == entity_id,
    ).limit(max_nodes // 2).all()
    
    # Get children (entities this one flows into)
    children = db.query(DataLineageEdge).filter(
        DataLineageEdge.source_type == entity_type,
        DataLineageEdge.source_id == entity_id,
    ).limit(max_nodes // 2).all()
    
    return {
        "parents": [
            {
                "type": edge.source_type,
                "id": str(edge.source_id),
                "relationship": edge.relationship_type,
                "transformation": edge.transformation,
            }
            for edge in parents
        ],
        "children": [
            {
                "type": edge.target_type,
                "id": str(edge.target_id),
                "relationship": edge.relationship_type,
                "transformation": edge.transformation,
            }
            for edge in children
        ],
        "total_parents": len(parents),
        "total_children": len(children),
    }


# ============================================================================
# CATALOG REFRESH OPERATIONS
# ============================================================================

def refresh_catalog_stats(db: Session, entity_type: str) -> bool:
    """
    Refresh row count statistics for a catalog entry.
    
    Args:
        db: Database session
        entity_type: Entity type to refresh
    
    Returns:
        True if refreshed, False if entry not found
    """
    entry = get_catalog_entry(db, entity_type)
    if not entry:
        return False
    
    try:
        # Get row count from actual table
        result = db.execute(text(f"SELECT COUNT(*) FROM {entry.table_name}"))
        row_count = result.scalar()
        
        entry.row_count = row_count
        entry.last_refreshed = datetime.utcnow()
        db.commit()
        
        logger.info(f"Refreshed {entity_type}: {row_count} rows")
        return True
    except Exception as e:
        logger.error(f"Failed to refresh {entity_type}: {e}")
        return False


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def get_catalog_stats(db: Session) -> dict:
    """
    Get overall catalog statistics.
    
    Returns:
        Dict with entry counts by category and total rows
    """
    entries = db.query(CatalogEntry).filter(CatalogEntry.is_active.is_(True)).all()
    
    stats = {
        "total_entries": len(entries),
        "total_rows": sum(e.row_count or 0 for e in entries),
        "categories": {},
        "last_refresh": None,
    }
    
    # Group by category
    for entry in entries:
        category = entry.category or "Unknown"
        if category not in stats["categories"]:
            stats["categories"][category] = {"count": 0, "rows": 0}
        
        stats["categories"][category]["count"] += 1
        stats["categories"][category]["rows"] += entry.row_count or 0
        
        # Track most recent refresh
        if entry.last_refreshed:
            if not stats["last_refresh"] or entry.last_refreshed > stats["last_refresh"]:
                stats["last_refresh"] = entry.last_refreshed
    
    return stats


def search_glossary(db: Session, query: str, limit: int = 20) -> List[dict]:
    """
    Search glossary terms and return simplified results.
    
    Args:
        db: Database session
        query: Search query
        limit: Maximum results
    
    Returns:
        List of term dicts with matched context
    """
    search_pattern = f"%{query}%"
    
    terms = db.query(GlossaryTerm).filter(
        or_(
            GlossaryTerm.term.ilike(search_pattern),
            GlossaryTerm.definition.ilike(search_pattern),
        )
    ).limit(limit).all()
    
    results = []
    for term in terms:
        # Highlight matched context
        match_type = "term" if query.lower() in term.term.lower() else "definition"
        
        results.append({
            "id": str(term.id),
            "term": term.term,
            "definition": term.definition[:200] + "..." if len(term.definition) > 200 else term.definition,
            "category": term.category,
            "match_type": match_type,
            "synonyms": term.synonyms or [],
        })
    
    return results


# ============================================================================
# AUTO-DISCOVERY ENGINE
# ============================================================================

# Entity category mapping
ENTITY_CATEGORIES = {
    # Core
    "Dataset": "Core",
    "Experiment": "Core",
    "Program": "Core",
    "Feature": "Core",
    "Signature": "Core",
    # Chemistry
    "Compound": "Chemistry",
    "Sample": "Chemistry",
    "HTSCampaign": "Chemistry",
    "HTSPlate": "Chemistry",
    "HTSWell": "Chemistry",
    "HTSResult": "Chemistry",
    "CompoundPlate": "Chemistry",
    "CompoundRequest": "Chemistry",
    # Omics
    "SingleCellDataset": "Omics",
    "LINCSSignature": "Omics",
    "Variant": "Omics",
    # Admin
    "User": "Admin",
    "Company": "Admin",
    "Team": "Admin",
    "AuditLog": "Admin",
    # Default
    "_default": "Other",
}


def refresh_catalog(db: Session) -> int:
    """
    Auto-discover tables and columns from SQLAlchemy models.
    
    Returns count of catalog entries created/updated.
    
    Uses SQLAlchemy's inspect() to introspect:
    - Table names
    - Column names, types, nullable, PKs, FKs
    - Row counts (cached)
    """
    from amprenta_rag.database.base import Base
    
    count = 0
    mapper_registry = Base.registry.mappers
    
    for mapper in mapper_registry:
        model_class = mapper.class_
        table = mapper.local_table
        
        if table is None:
            continue
            
        entity_type = model_class.__name__
        table_name = table.name
        
        # Skip internal tables
        if table_name.startswith('_') or table_name in ('alembic_version',):
            continue
        
        # Get or create catalog entry
        entry = db.query(CatalogEntry).filter(CatalogEntry.entity_type == entity_type).first()
        
        if not entry:
            entry = CatalogEntry(
                entity_type=entity_type,
                table_name=table_name,
                display_name=_humanize(entity_type),
                category=ENTITY_CATEGORIES.get(entity_type, ENTITY_CATEGORIES["_default"]),
            )
            db.add(entry)
            db.flush()
        
        # Update row count
        try:
            entry.row_count = db.query(model_class).count()
        except Exception:
            entry.row_count = None
        
        entry.last_refreshed = func.now()
        
        # Discover columns
        inspector = inspect(model_class)
        existing_columns = {c.column_name for c in entry.columns}
        
        for column in inspector.columns:
            if column.name in existing_columns:
                continue
                
            col_meta = ColumnMetadata(
                catalog_entry_id=entry.id,
                column_name=column.name,
                display_name=_humanize(column.name),
                data_type=_extract_column_type(column.type),
                is_nullable=column.nullable or False,
                is_primary_key=column.primary_key,
                is_foreign_key=len(column.foreign_keys) > 0,
                foreign_key_target=_get_fk_target(column) if column.foreign_keys else None,
            )
            db.add(col_meta)
        
        count += 1
    
    db.commit()
    return count


def _humanize(name: str) -> str:
    """Convert CamelCase or snake_case to human readable."""
    import re
    # CamelCase to spaces
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1 \2', name)
    # snake_case to spaces
    s2 = re.sub('([a-z0-9])([A-Z])', r'\1 \2', s1)
    return s2.replace('_', ' ').title()


def _extract_column_type(sa_type) -> str:
    """Convert SQLAlchemy type to catalog type string."""
    type_name = type(sa_type).__name__.lower()
    
    type_map = {
        'uuid': 'uuid',
        'string': 'string',
        'varchar': 'string',
        'text': 'text',
        'integer': 'integer',
        'biginteger': 'integer',
        'smallinteger': 'integer',
        'float': 'float',
        'numeric': 'decimal',
        'boolean': 'boolean',
        'datetime': 'datetime',
        'date': 'date',
        'time': 'time',
        'json': 'json',
        'jsonb': 'json',
        'array': 'array',
        'enum': 'enum',
    }
    
    return type_map.get(type_name, 'unknown')


def _get_fk_target(column) -> Optional[str]:
    """Get foreign key target as 'table.column' string."""
    if not column.foreign_keys:
        return None
    fk = list(column.foreign_keys)[0]
    return str(fk.target_fullname)


def get_sample_values(db: Session, table_name: str, column_name: str, limit: int = 5) -> List[str]:
    """
    Get sample values for a column.
    
    Safety: Uses parameterized query and LIMIT at DB level.
    """
    # Safety: validate table/column names (alphanumeric + underscore only)
    import re
    if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$', table_name):
        return []
    if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$', column_name):
        return []
    
    try:
        query = text(f"""
            SELECT DISTINCT "{column_name}"::text 
            FROM "{table_name}" 
            WHERE "{column_name}" IS NOT NULL 
            LIMIT :limit
        """)
        result = db.execute(query, {"limit": limit})
        return [str(row[0]) for row in result if row[0] is not None]
    except Exception:
        return []


def detect_lineage_from_fks(db: Session) -> int:
    """
    Auto-detect lineage edges from foreign key relationships.
    
    Creates 'references' relationship type edges.
    Returns count of edges created.
    """
    count = 0
    
    for entry in db.query(CatalogEntry).all():
        for col in entry.columns:
            if col.is_foreign_key and col.foreign_key_target:
                # Parse target: "table.column" or "schema.table.column"
                parts = col.foreign_key_target.split('.')
                target_table = parts[-2] if len(parts) >= 2 else parts[0]
                
                # Find target catalog entry
                target_entry = db.query(CatalogEntry).filter(
                    CatalogEntry.table_name == target_table
                ).first()
                
                if target_entry:
                    # Check if edge exists
                    existing = db.query(DataLineageEdge).filter(
                        DataLineageEdge.source_type == entry.entity_type,
                        DataLineageEdge.target_type == target_entry.entity_type,
                        DataLineageEdge.relationship_type == "references"
                    ).first()
                    
                    if not existing:
                        # Create edge (source references target)
                        edge = DataLineageEdge(
                            source_type=entry.entity_type,
                            source_id=entry.id,  # Use catalog entry ID for schema-level lineage
                            target_type=target_entry.entity_type,
                            target_id=target_entry.id,
                            relationship_type="references",
                            transformation=f"FK: {col.column_name}",
                        )
                        db.add(edge)
                        count += 1
    
    db.commit()
    return count
