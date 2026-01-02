# UI Feature Gaps - Comprehensive Fix (v2)

**Updated based on Reviewer feedback: Added E2E tests, file upload security**Address all identified UI gaps with complete test coverage.

## Scope

| Gap | Solution | Batch ||-----|----------|-------|| Data Lifecycle E2E tests | Add 6 E2E tests for existing dashboard | 0 || Lifecycle stats endpoint | Add `GET /lifecycle/stats` | 1 || Lifecycle audit endpoint | Add `GET /lifecycle/audit` | 1 || OCR/Scraper/Normalizer UI | New "AI Extraction Tools" page | 2 || Conflict Resolver UI | Add "Conflicts" tab to Mapping Refresh | 3 || KEGG Refresh visibility | Add "KEGG Cache" tab to Mapping Refresh | 3 || New page E2E tests | 5 E2E tests for new UI components | 4 |---

## Batch 0: Data Lifecycle E2E Tests (ADDED per Reviewer)

**New File**: `amprenta_rag/tests/e2e/test_data_lifecycle_e2e.py`6 E2E tests covering the existing Data Lifecycle dashboard:

```python
def test_lifecycle_page_loads():
    """Page renders with 4 tabs."""

def test_quarantine_check_impact():
    """Check Impact button returns impact data."""

def test_quarantine_restore_entity():
    """Restore action updates entity status."""

def test_bulk_preview_shows_impact():
    """Preview button shows aggregated impact."""

def test_bulk_execute_requires_preview():
    """Execute button disabled until preview done."""

def test_audit_tab_shows_history():
    """Audit tab displays lifecycle changes."""
```

---

## Batch 1: Data Lifecycle Dashboard Completion

**API Endpoints** - Add to `amprenta_rag/api/routers/lifecycle.py`:

```python
@router.get("/stats")
def get_lifecycle_stats() -> dict:
    """Get entity counts grouped by lifecycle_status."""
    # Query each entity type, group by lifecycle_status
    # Return: {"dataset": {"active": 100, "quarantined": 5, ...}, ...}

@router.get("/audit")
def get_lifecycle_audit(
    entity_type: Optional[str] = None,
    status_change: Optional[str] = None,
    days: int = 7,
    skip: int = 0,
    limit: int = 50,
) -> dict:
    """Query AuditLog for lifecycle_status_change actions."""
```

**Dashboard Updates** - Update `scripts/dashboard/pages/data_lifecycle.py`:

- Overview tab: Replace placeholder with real stats from `/lifecycle/stats`
- Audit tab: Replace placeholder with paginated results from `/lifecycle/audit`

---

## Batch 2: AI Extraction Tools Dashboard (NEW)

**New Page**: `scripts/dashboard/pages/ai_extraction.py`3-tab interface:

```javascript
Tab 1: OCR
├── File upload (PDF/image, max 10MB)
├── Language selector (eng, fra, deu, etc.)
├── "Extract Text" button
└── Results display with word count

Tab 2: Web Scraper  
├── URL input field
├── "Scrape" button
├── Content preview (title, author, word count)
└── Full text expandable section

Tab 3: Entity Normalizer
├── Entity type selector (gene, compound, disease)
├── Text input for entity name
├── "Lookup" button
└── Normalized result display (UniProt ID, PubChem CID, etc.)
```

**New API Endpoints** - Create `amprenta_rag/api/routers/extraction.py`:

```python
MAX_FILE_SIZE = 10_000_000  # 10MB
ALLOWED_CONTENT_TYPES = ["application/pdf", "image/png", "image/jpeg", "image/tiff"]

@router.post("/ocr")
async def extract_ocr(file: UploadFile) -> dict:
    """Extract text from uploaded file using OCR.
    
    P1 Security:
    - File size limit: 10MB
    - Content type validation: PDF, PNG, JPEG, TIFF only
    """
    if file.size and file.size > MAX_FILE_SIZE:
        raise HTTPException(413, "File too large (max 10MB)")
    if file.content_type not in ALLOWED_CONTENT_TYPES:
        raise HTTPException(415, f"Unsupported file type: {file.content_type}")
    # ... OCR logic

@router.post("/scrape")
def scrape_url(url: str) -> dict:
    """Scrape content from URL.
    
    P1 Security:
    - Inherits rate limiting from WebScraper (0.5s per domain)
    """

@router.post("/normalize")
def normalize_entity(entity_type: str, name: str) -> dict:
    """Normalize entity name to standard identifiers."""
```

---

## Batch 3: Sync Management (Extend Mapping Refresh)

**Update**: `scripts/dashboard/pages/mapping_refresh.py`Add 2 new tabs (total 6 tabs):

```javascript
Existing tabs: Status, Statistics, Lookup, Jobs
New tabs:
├── Tab 5: Conflicts
│   ├── Pending conflicts table
│   ├── Conflict detail view (local vs external values)
│   ├── Resolution actions (prefer_external, prefer_local, ignore)
│   └── Uses existing GET /sync/conflicts and POST /sync/conflicts/{id}/resolve
│
└── Tab 6: KEGG Cache
    ├── Expiring mappings count (next 7 days)
    ├── Last refresh timestamp
    ├── Manual refresh trigger button
    └── Refresh history log
```

**New API Endpoint** - Add to `amprenta_rag/api/routers/mappings.py`:

```python
@router.get("/kegg/status")
def get_kegg_cache_status() -> dict:
    """Get KEGG cache status: expiring count, last refresh, etc."""
```

---

## Batch 4: Tests (EXPANDED per Reviewer)

### API Tests (8 tests)

- `test_lifecycle_stats_endpoint`
- `test_lifecycle_audit_endpoint`
- `test_ocr_endpoint_success`
- `test_ocr_endpoint_file_too_large`
- `test_ocr_endpoint_invalid_type`
- `test_scrape_endpoint`
- `test_normalize_endpoint`
- `test_kegg_status_endpoint`

### Import Tests (3 tests)

- Verify all new dashboard pages import without error

### E2E Tests for New Pages (5 tests) - ADDED per Reviewer

**AI Extraction Tools** (`test_ai_extraction_e2e.py`):

```python
def test_extraction_page_loads()
def test_ocr_tab_file_upload()
def test_normalizer_lookup()
```

**Mapping Refresh Extensions** (`test_mapping_refresh_extended_e2e.py`):

```python
def test_conflicts_tab_shows_pending()
def test_kegg_cache_tab_shows_status()
```

---

## Test Summary (Updated)

| Category | Count ||----------|-------|| Batch 0: Data Lifecycle E2E | 6 || Batch 4: API tests | 8 || Batch 4: Import tests | 3 || Batch 4: New page E2E | 5 || Batch 4: Lifecycle API extensions | 3 || **Total** | **25** |---

## File Summary

| File | Action ||------|--------|| `amprenta_rag/tests/e2e/test_data_lifecycle_e2e.py` | NEW - 6 E2E tests || `amprenta_rag/api/routers/lifecycle.py` | Add 2 endpoints || `amprenta_rag/api/routers/extraction.py` | NEW - 3 endpoints || `amprenta_rag/api/routers/mappings.py` | Add 1 endpoint || `amprenta_rag/api/main.py` | Register extraction router || `scripts/dashboard/pages/data_lifecycle.py` | Update 2 tabs || `scripts/dashboard/pages/ai_extraction.py` | NEW - 3 tabs || `scripts/dashboard/pages/mapping_refresh.py` | Add 2 tabs || `scripts/dashboard/core/config.py` | Register AI Extraction page || `amprenta_rag/tests/api/test_extraction_api.py` | NEW - 5 tests || `amprenta_rag/tests/api/test_lifecycle_api.py` | Add 2 tests || `amprenta_rag/tests/e2e/test_ai_extraction_e2e.py` | NEW - 3 tests || `amprenta_rag/tests/e2e/test_mapping_refresh_extended_e2e.py` | NEW - 2 tests |---

## P2 Observations (Deferred)

1. **Bulk restore API**: Consider `POST /lifecycle/bulk/restore` for symmetry
2. **SSRF prevention**: Block internal IPs in scraper