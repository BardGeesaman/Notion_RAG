# BAM/CRAM Alignment Viewing

Extend the genomics module with binary alignment file (BAM/CRAM) support for viewing sequence alignments, coverage statistics, and read-level detail.

**Status**: Approved with P1 fixes incorporated

---

## Architecture

```
Upload Flow:
  Dashboard Upload → POST /alignments/upload (file + index) → bam_parser.py → AlignmentFile Model → DB

View Flow:
  Browse Tab → GET /alignments
  Region View → GET /alignments/{id}/reads?region=chr:start-end&page=1
  Coverage → GET /alignments/{id}/coverage?region=chr:start-end
```

---

## Scope

| Component | File | Complexity |
|-----------|------|------------|
| BAM/CRAM Parser | `amprenta_rag/ingestion/genomics/bam_parser.py` | Medium |
| AlignmentFile Model | `amprenta_rag/models/misc.py` | Low |
| Migration | `alembic/versions/xxxx_add_alignment_files.py` | Low |
| API Endpoints (5) | `amprenta_rag/api/routers/genomics.py` | Medium |
| Dashboard UI | `scripts/dashboard/pages/alignments.py` | Medium |
| Tests (12) | `tests/` | Medium |

---

## Batch 1: Parser + Model + Migration

### 1. BAM Parser (`amprenta_rag/ingestion/genomics/bam_parser.py`)

```python
def parse_bam_header(bam_path: Path) -> dict:
    """Extract header info: reference sequences, read groups, programs."""
    
def get_alignment_stats(bam_path: Path) -> AlignmentStats:
    """Calculate: total reads, mapped, unmapped, duplicate rate, coverage."""
    
def fetch_reads(bam_path: Path, region: str, offset: int, limit: int) -> list[dict]:
    """Fetch reads in region (chr:start-end format) with pagination."""
    
def get_coverage(bam_path: Path, region: str, bin_size: int = 100) -> list[int]:
    """Calculate coverage histogram for a region."""
    
def check_index_exists(bam_path: Path) -> bool:
    """Check if .bai/.crai index file exists."""
```

### 2. AlignmentFile Model (extend `models/misc.py`)

```python
class AlignmentFile(Base):
    __tablename__ = "alignment_files"
    
    id = Column(UUID, primary_key=True)
    filename = Column(String(500), nullable=False)
    file_path = Column(String(1000), nullable=False)
    index_file_path = Column(String(1000), nullable=True)  # P1 FIX
    file_format = Column(String(10), nullable=False)  # "BAM" or "CRAM"
    has_index = Column(Boolean, default=False)  # P1 FIX
    
    # Header metadata
    reference_genome = Column(String(100))
    num_references = Column(Integer)
    read_groups = Column(JSON)
    
    # Stats
    total_reads = Column(BigInteger)
    mapped_reads = Column(BigInteger)
    unmapped_reads = Column(BigInteger)
    duplicate_rate = Column(Float)
    mean_coverage = Column(Float)
    
    # Relationships
    experiment_id = Column(UUID, ForeignKey("experiments.id"), nullable=True)
    created_by_id = Column(UUID, ForeignKey("users.id"))
    created_at = Column(DateTime, server_default=func.now())
```

### 3. Alembic Migration

- Add `alignment_files` table
- Indexes on `experiment_id`, `created_by_id`, `has_index`

### 4. Storage Configuration (P1 FIX)

Environment variable: `ALIGNMENT_STORAGE_PATH`
- Default: `./data/alignments/`
- Production: S3 path or mounted volume

File organization:
```
{ALIGNMENT_STORAGE_PATH}/{alignment_id}/
  ├── {filename}         # BAM/CRAM file
  └── {filename}.bai     # Index file (or .crai)
```

---

## Batch 2: API Endpoints

Extend `amprenta_rag/api/routers/genomics.py`:

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/alignments/upload` | POST | Upload BAM/CRAM + index file (P1 FIX) |
| `/alignments` | GET | List alignment files with filters |
| `/alignments/{id}` | GET | Get alignment file details |
| `/alignments/{id}/reads` | GET | Fetch reads in region with pagination (P1 FIX) |
| `/alignments/{id}/coverage` | GET | Get coverage histogram |

### Upload Endpoint (P1 FIX: Accept index file)

```python
@router.post("/alignments/upload")
async def upload_alignment(
    file: UploadFile = File(...),              # BAM/CRAM
    index_file: UploadFile = File(None),       # .bai/.crai (optional)
    experiment_id: Optional[UUID] = Query(None),
    current_user: User = Depends(get_current_user),
) -> AlignmentUploadResponse:
```

### Reads Endpoint (P1 FIX: Pagination)

```python
@router.get("/alignments/{id}/reads")
def fetch_reads(
    id: UUID,
    region: str = Query(..., description="chr:start-end"),
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=1000),
    current_user: User = Depends(get_current_user),
) -> ReadResponse:
```

### Security/Validation

- Magic byte validation: BAM (`BAM\x01`), CRAM (`CRAM`)
- File size limit: 500MB (configurable via `MAX_ALIGNMENT_SIZE`)
- Region format validation: `^chr[\dXYMT]+:\d+-\d+$`
- Read count limits: max 5000 per response
- Index required for region queries (return 400 if missing)

---

## Batch 3: Dashboard UI

**New Page**: `scripts/dashboard/pages/alignments.py`

### 4-Tab Layout (P1 FIX: All tabs must have E2E coverage)

| Tab | Content |
|-----|---------|
| **Browse** | List uploaded alignments, filter by experiment, download links |
| **Upload** | File uploaders (BAM + index), experiment selector, progress |
| **View** | Region input, coverage chart (Plotly), read pileup table |
| **Stats** | Summary metrics, read group info, reference sequences |

### Register in Config

Add to `scripts/dashboard/core/config.py`:
- `PAGE_REGISTRY["Alignments"]`
- Add to "Discovery" group

---

## Batch 4: Tests (12 total)

### Parser Tests (4)

```python
# amprenta_rag/tests/ingestion/test_bam_parser.py
def test_parse_header_extracts_references()
def test_get_alignment_stats_calculates_correctly()
def test_fetch_reads_returns_paginated_results()
def test_get_coverage_returns_histogram()
```

### API Tests (4)

```python
# amprenta_rag/tests/api/test_genomics_api.py
def test_alignment_upload_validates_format()
def test_alignment_list_returns_filtered()
def test_alignment_reads_requires_index()
def test_alignment_coverage_returns_data()
```

### E2E Tests (4) - P1 FIX: Cover all 4 tabs

```python
# amprenta_rag/tests/e2e/test_alignments_e2e.py
def test_alignment_page_loads()           # All 4 tabs visible
def test_alignment_upload_form()          # Upload tab works
def test_alignment_region_view()          # View tab with region input
def test_alignment_stats_display()        # Stats tab shows metrics
```

---

## Verification Checkpoints

```bash
# After Batch 1
ruff check amprenta_rag/ingestion/genomics/bam_parser.py amprenta_rag/models/misc.py
alembic upgrade head

# After Batch 2
python -c "from amprenta_rag.api.routers.genomics import router; print(f'Routes: {len(router.routes)}')"

# After Batch 3
python -c "from scripts.dashboard.core.config import PAGE_REGISTRY; print('Alignments' in PAGE_REGISTRY)"

# After Batch 4
pytest amprenta_rag/tests/ingestion/test_bam_parser.py -v
pytest amprenta_rag/tests/api/test_genomics_api.py -v -k alignment
pytest amprenta_rag/tests/e2e/test_alignments_e2e.py -v
```

---

## P2 Observations (Future Work)

- Auto-index generation: Celery task to run `pysam.index()` if no index provided
- CRAM reference genome: May need reference FASTA handling
- IGV.js integration for genome browser visualization
- Streaming JSON response for large regions
- Soft delete with file cleanup from storage
- Per-request file handle opening for concurrent access

