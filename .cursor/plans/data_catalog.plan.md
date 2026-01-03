# Data Catalog System

## Overview

A foundational data governance feature providing entity discovery, column-level metadata, data lineage visualization, and business glossary integration.

**Priority:** P1  
**Effort:** 1 session (MVP)  
**Dependencies:** None

## Goals

1. Help scientists discover what data exists in the platform
2. Understand column meanings, types, and relationships
3. Visualize data lineage (source → transformations → outputs)
4. Bridge business terminology with technical field names

---

## Batch 1: Database Models + Service Foundation

### Models (`amprenta_rag/models/catalog.py`)

```python
class CatalogEntry(Base):
    """Represents a registered entity type (table) in the catalog."""
    __tablename__ = "catalog_entries"
    
    id: Mapped[UUID] = mapped_column(primary_key=True, default=uuid4)
    entity_type: Mapped[str] = mapped_column(String(100), unique=True)  # e.g., "Dataset", "Compound"
    table_name: Mapped[str] = mapped_column(String(100))  # e.g., "datasets", "compounds"
    display_name: Mapped[str] = mapped_column(String(200))
    description: Mapped[Optional[str]] = mapped_column(Text)
    category: Mapped[str] = mapped_column(String(50))  # "Core", "Chemistry", "Omics", "Admin"
    row_count: Mapped[Optional[int]] = mapped_column(Integer)  # Cached count
    last_refreshed: Mapped[Optional[datetime]] = mapped_column(DateTime)
    is_active: Mapped[bool] = mapped_column(Boolean, default=True)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    
    # Relationships
    columns: Mapped[List["ColumnMetadata"]] = relationship(back_populates="catalog_entry")


class ColumnMetadata(Base):
    """Column-level metadata for catalog entries."""
    __tablename__ = "column_metadata"
    
    id: Mapped[UUID] = mapped_column(primary_key=True, default=uuid4)
    catalog_entry_id: Mapped[UUID] = mapped_column(ForeignKey("catalog_entries.id"))
    column_name: Mapped[str] = mapped_column(String(100))
    display_name: Mapped[Optional[str]] = mapped_column(String(200))
    data_type: Mapped[str] = mapped_column(String(50))  # "string", "integer", "uuid", "datetime", "json"
    description: Mapped[Optional[str]] = mapped_column(Text)
    is_nullable: Mapped[bool] = mapped_column(Boolean, default=True)
    is_primary_key: Mapped[bool] = mapped_column(Boolean, default=False)
    is_foreign_key: Mapped[bool] = mapped_column(Boolean, default=False)
    foreign_key_target: Mapped[Optional[str]] = mapped_column(String(200))  # "table.column"
    example_values: Mapped[Optional[str]] = mapped_column(Text)  # JSON array of examples
    glossary_term_id: Mapped[Optional[UUID]] = mapped_column(ForeignKey("glossary_terms.id"))
    
    # Relationships
    catalog_entry: Mapped["CatalogEntry"] = relationship(back_populates="columns")
    glossary_term: Mapped[Optional["GlossaryTerm"]] = relationship()


class GlossaryTerm(Base):
    """Business glossary terms with definitions."""
    __tablename__ = "glossary_terms"
    
    id: Mapped[UUID] = mapped_column(primary_key=True, default=uuid4)
    term: Mapped[str] = mapped_column(String(200), unique=True)  # e.g., "IC50"
    definition: Mapped[str] = mapped_column(Text)
    category: Mapped[Optional[str]] = mapped_column(String(100))  # "Chemistry", "Biology", "Statistics"
    synonyms: Mapped[Optional[str]] = mapped_column(Text)  # JSON array
    related_terms: Mapped[Optional[str]] = mapped_column(Text)  # JSON array of term IDs
    source: Mapped[Optional[str]] = mapped_column(String(200))  # "Internal", "PubChem", "UniProt"
    created_by_id: Mapped[Optional[UUID]] = mapped_column(ForeignKey("users.id"))
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    updated_at: Mapped[datetime] = mapped_column(DateTime, default=func.now(), onupdate=func.now())


class DataLineageEdge(Base):
    """Represents data flow between entities."""
    __tablename__ = "data_lineage_edges"
    
    id: Mapped[UUID] = mapped_column(primary_key=True, default=uuid4)
    source_type: Mapped[str] = mapped_column(String(100))  # "Dataset", "Experiment"
    source_id: Mapped[UUID] = mapped_column(PostgresUUID)
    target_type: Mapped[str] = mapped_column(String(100))
    target_id: Mapped[UUID] = mapped_column(PostgresUUID)
    relationship_type: Mapped[str] = mapped_column(String(50))  # "derived_from", "input_to", "generated_by"
    transformation: Mapped[Optional[str]] = mapped_column(String(200))  # "normalization", "batch_correction"
    metadata: Mapped[Optional[dict]] = mapped_column(JSONB)  # Additional context
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    
    __table_args__ = (
        Index("ix_lineage_source", "source_type", "source_id"),
        Index("ix_lineage_target", "target_type", "target_id"),
    )
```

### Service Layer (`amprenta_rag/services/catalog.py`)

```python
# Core functions:
def refresh_catalog() -> int:
    """Auto-discover tables and columns from SQLAlchemy models. Returns count of entries updated."""

def get_catalog_entries(category: Optional[str] = None, search: Optional[str] = None) -> List[CatalogEntry]:
    """List catalog entries with optional filtering."""

def get_catalog_entry(entity_type: str) -> Optional[CatalogEntry]:
    """Get single catalog entry with columns."""

def update_catalog_entry(entity_type: str, description: str, category: str) -> CatalogEntry:
    """Update catalog entry metadata."""

def search_columns(query: str, limit: int = 50) -> List[ColumnMetadata]:
    """Search columns by name or description across all entities."""

def get_column_metadata(entity_type: str, column_name: str) -> Optional[ColumnMetadata]:
    """Get detailed column metadata."""

def update_column_metadata(column_id: UUID, display_name: str, description: str, glossary_term_id: Optional[UUID]) -> ColumnMetadata:
    """Update column metadata and glossary link."""

# Glossary functions:
def create_glossary_term(term: str, definition: str, category: Optional[str], synonyms: List[str]) -> GlossaryTerm:
    """Create new glossary term."""

def get_glossary_terms(category: Optional[str] = None, search: Optional[str] = None) -> List[GlossaryTerm]:
    """List glossary terms with filtering."""

def update_glossary_term(term_id: UUID, definition: str, category: str, synonyms: List[str]) -> GlossaryTerm:
    """Update glossary term."""

def delete_glossary_term(term_id: UUID) -> bool:
    """Delete glossary term (unlinks from columns first)."""

# Lineage functions:
def add_lineage_edge(source_type: str, source_id: UUID, target_type: str, target_id: UUID, relationship_type: str, transformation: Optional[str] = None) -> DataLineageEdge:
    """Record data lineage relationship."""

def get_lineage_graph(entity_type: str, entity_id: UUID, depth: int = 3, direction: str = "both") -> dict:
    """Get lineage graph for an entity (upstream/downstream/both). Returns Cytoscape.js JSON."""

def get_entity_lineage(entity_type: str, entity_id: UUID) -> dict:
    """Get immediate lineage (parents + children) for an entity."""
```

### Alembic Migration

- Create `catalog_entries`, `column_metadata`, `glossary_terms`, `data_lineage_edges` tables
- Add indexes for search performance

### Tests (Batch 1)

- 5 model tests (relationships, constraints)
- 8 service tests (CRUD, search, lineage graph)

---

## Batch 2: Auto-Discovery Engine

### Schema Introspection (`amprenta_rag/services/catalog.py`)

```python
def _discover_models() -> List[dict]:
    """Introspect SQLAlchemy models to build catalog entries."""
    # Uses SQLAlchemy's inspect() to get:
    # - Table names
    # - Column names, types, nullable, PKs, FKs
    # - Relationships
    
def _extract_column_type(sa_type) -> str:
    """Convert SQLAlchemy type to catalog type string."""
    
def _get_sample_values(table_name: str, column_name: str, limit: int = 5) -> List[str]:
    """Get sample values for a column (for discovery UI)."""

def _count_rows(table_name: str) -> int:
    """Get row count for table (cached)."""

# Category mapping
ENTITY_CATEGORIES = {
    "Dataset": "Core",
    "Experiment": "Core", 
    "Compound": "Chemistry",
    "Sample": "Chemistry",
    "Signature": "Omics",
    "Feature": "Omics",
    "User": "Admin",
    "Program": "Admin",
    # ... etc
}
```

### Lineage Auto-Detection

```python
def _detect_lineage_from_fks() -> List[dict]:
    """Auto-detect lineage edges from foreign key relationships."""

def _detect_lineage_from_audit() -> List[dict]:
    """Detect lineage from audit_logs (derived_from events)."""
```

### Tests (Batch 2)

- 4 introspection tests
- 3 lineage detection tests

---

## Batch 3: API Endpoints

### Router (`amprenta_rag/api/routers/catalog.py`)

```python
# Catalog Browsing
GET  /api/v1/catalog/entries                    # List all catalog entries
GET  /api/v1/catalog/entries/{entity_type}      # Get single entry with columns
PUT  /api/v1/catalog/entries/{entity_type}      # Update entry metadata
POST /api/v1/catalog/refresh                    # Trigger catalog refresh

# Column Search
GET  /api/v1/catalog/columns/search             # Search columns across entities
GET  /api/v1/catalog/columns/{entity_type}/{column_name}  # Get column details
PUT  /api/v1/catalog/columns/{column_id}        # Update column metadata

# Glossary
GET  /api/v1/catalog/glossary                   # List glossary terms
POST /api/v1/catalog/glossary                   # Create term
GET  /api/v1/catalog/glossary/{term_id}         # Get term
PUT  /api/v1/catalog/glossary/{term_id}         # Update term
DELETE /api/v1/catalog/glossary/{term_id}       # Delete term

# Lineage
GET  /api/v1/catalog/lineage/{entity_type}/{entity_id}  # Get lineage graph
POST /api/v1/catalog/lineage                    # Add lineage edge (manual)
```

### Schemas (`amprenta_rag/api/schemas.py`)

```python
class CatalogEntryResponse(BaseModel):
    id: UUID
    entity_type: str
    table_name: str
    display_name: str
    description: Optional[str]
    category: str
    row_count: Optional[int]
    column_count: int
    last_refreshed: Optional[datetime]

class ColumnMetadataResponse(BaseModel):
    id: UUID
    column_name: str
    display_name: Optional[str]
    data_type: str
    description: Optional[str]
    is_nullable: bool
    is_primary_key: bool
    is_foreign_key: bool
    foreign_key_target: Optional[str]
    example_values: Optional[List[str]]
    glossary_term: Optional[GlossaryTermBrief]

class GlossaryTermCreate(BaseModel):
    term: str
    definition: str
    category: Optional[str]
    synonyms: Optional[List[str]]

class GlossaryTermResponse(BaseModel):
    id: UUID
    term: str
    definition: str
    category: Optional[str]
    synonyms: Optional[List[str]]
    related_terms: Optional[List[str]]
    source: Optional[str]
    created_at: datetime

class LineageGraphResponse(BaseModel):
    nodes: List[dict]  # Cytoscape.js nodes
    edges: List[dict]  # Cytoscape.js edges
    center_entity: dict

class LineageEdgeCreate(BaseModel):
    source_type: str
    source_id: UUID
    target_type: str
    target_id: UUID
    relationship_type: str
    transformation: Optional[str]
```

### Tests (Batch 3)

- 10 API tests covering all endpoints

---

## Batch 4: Dashboard UI

### Page (`scripts/dashboard/pages/data_catalog.py`)

**4-Tab Layout:**

#### Tab 1: Browse Entities
- Card grid of entity types grouped by category (Core, Chemistry, Omics, Admin)
- Each card shows: entity name, description, row count, column count
- Click to expand column details
- Search/filter by category

#### Tab 2: Column Search
- Global search box for column names/descriptions
- Results table: Entity | Column | Type | Description | Glossary Term
- Click column to see full details + sample values
- "Where is this column used?" - shows all tables with matching column names

#### Tab 3: Data Lineage
- Entity type + ID selector (or deep-link from other pages)
- Cytoscape.js graph visualization
- Upstream (where data came from) in blue
- Downstream (where data flows to) in green
- Center entity highlighted
- Click nodes to navigate
- Depth slider (1-5 hops)
- Export lineage as PNG/JSON

#### Tab 4: Business Glossary
- Alphabetical term list with search
- Term cards: term, definition, category, synonyms
- Add/Edit/Delete terms (admin only)
- "Link to Column" button to associate with columns
- Import/Export glossary (CSV)

### Deep-Link Support

```python
# URL parameters for cross-page linking
?entity_type=Dataset&entity_id=xxx  # Opens lineage for specific entity
?search=patient_id                  # Opens column search with query
?term=IC50                          # Opens glossary to specific term
```

### Integration Points

- Add "📖 View in Catalog" button to entity detail pages
- Add "🌳 View Lineage" button to Dataset/Experiment pages
- Show glossary tooltip on hover for linked columns

---

## Batch 5: Tests

### Service Tests (`amprenta_rag/tests/services/test_catalog_service.py`)

```python
# Catalog Entry Tests (5)
- test_refresh_catalog_discovers_models
- test_get_catalog_entries_filters_by_category
- test_search_catalog_entries
- test_update_catalog_entry
- test_get_catalog_entry_with_columns

# Column Tests (4)
- test_search_columns_by_name
- test_search_columns_by_description
- test_update_column_metadata
- test_link_column_to_glossary

# Glossary Tests (5)
- test_create_glossary_term
- test_get_glossary_terms_with_search
- test_update_glossary_term
- test_delete_glossary_term_unlinks_columns
- test_glossary_synonyms_searchable

# Lineage Tests (4)
- test_add_lineage_edge
- test_get_lineage_graph_upstream
- test_get_lineage_graph_downstream
- test_get_lineage_graph_bidirectional
```

### API Tests (`amprenta_rag/tests/api/test_catalog_api.py`)

- 10 API endpoint tests

---

## Implementation Order

| Batch | Content | Estimated |
|-------|---------|-----------|
| 1 | Models + Service foundation | 30 min |
| 2 | Auto-discovery engine | 20 min |
| 3 | API endpoints (12 endpoints) | 25 min |
| 4 | Dashboard UI (4 tabs) | 35 min |
| 5 | Tests (28 tests) | 25 min |

**Total: ~2.5 hours**

---

## Success Criteria

1. ✅ Scientists can browse all entity types with schemas
2. ✅ Column search works across the platform
3. ✅ Lineage graph renders for any entity
4. ✅ Glossary terms can be created and linked to columns
5. ✅ Auto-refresh discovers new models
6. ✅ 28+ tests passing

---

## Future Enhancements (P2/P3)

- Data quality scores per column (completeness, uniqueness)
- Column-level usage analytics (which columns are queried most)
- Schema change detection and alerts
- Automatic glossary suggestions from LLM
- Data classification (PII, PHI, sensitive)
- Column-level access controls

