# Code Architecture & Patterns

**Last Updated:** 2025-01-01  
**Status:** Production guide based on comprehensive codebase audit

## Overview

This guide documents the established architectural patterns across the Amprenta multi-omics platform. It serves as the authoritative reference for maintaining consistency when adding new features or refactoring existing code.

The codebase follows a layered architecture with clear separation of concerns:
- **Services Layer**: Business logic and data operations
- **API Layer**: REST endpoints and request handling  
- **Database Layer**: Data models and persistence
- **Sync/Adapters**: External system integrations
- **ML Layer**: Machine learning and prediction services
- **Dashboard Layer**: Streamlit-based user interface
- **Jobs Layer**: Background task processing

Each layer has established patterns that should be followed for consistency and maintainability.

---

## Services Layer (`amprenta_rag/services/`)

**Pattern**: Functional (93%) with limited class-based for complex state management  
**Files**: 14 functional services, 1 class-based service

### Session Management
**Standard**: Context manager pattern with `db_session()`

```python
# ✅ PREFERRED: Context manager pattern (from activity.py)
def log_activity(
    event_type: ActivityEventType,
    target_type: str,
    target_id: UUID,
    target_name: str,
    actor_id: Optional[UUID] = None,
) -> Optional[ActivityEvent]:
    try:
        with db_session() as db:
            event = ActivityEvent(
                event_type=event_type,
                target_type=target_type,
                target_id=target_id,
                target_name=target_name,
                actor_id=actor_id,
            )
            db.add(event)
            db.commit()
            db.expunge(event)
            return event
    except Exception as e:
        logger.error(f"Failed to log activity: {e}")
        return None
```

**Alternative**: Dependency injection for API-called services

```python
# ✅ ACCEPTABLE: Dependency injection (from ip_service.py)
def create_disclosure(
    title: str,
    description: str,
    inventors: List[dict],
    user_id: UUID,
    company_id: UUID,
    db: Session,
) -> InventionDisclosure:
    disclosure = InventionDisclosure(
        title=title,
        description=description,
        created_by_id=user_id,
        company_id=company_id,
    )
    db.add(disclosure)
    db.commit()
    return disclosure
```

### Error Handling
**Standard**: Try/except with logging, return `Optional[Model]` or `None`

```python
# ✅ STANDARD: Error handling pattern (from entity_reviews.py)
def create_review(
    entity_type: str,
    entity_id: UUID,
    submitted_by_id: UUID,
    db: Optional[Session] = None
) -> Optional[EntityReview]:
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            # Business logic here
            review = EntityReview(
                entity_type=entity_type,
                entity_id=entity_id,
                reviewer_id=submitted_by_id,
                status="draft"
            )
            session.add(review)
            session.commit()
            session.expunge(review)
            return review
            
    except Exception as e:
        logger.error(f"Failed to create review: {e}")
        return None
```

### Class-based Services
**When to use**: Complex state management or multiple related operations

```python
# ✅ CLASS-BASED: For complex state (from notebook_review.py)
class NotebookReviewService:
    def __init__(self, db_session_factory):
        self.db_session = db_session_factory
        self._cache = {}
    
    def capture_snapshot(self, review_id: UUID, notebook_path: str) -> NotebookSnapshot:
        # Complex state management logic
        pass
```

---

## API Layer (`amprenta_rag/api/routers/`)

**Pattern**: 100% functional endpoints with dependency injection  
**Files**: ~25 routers with consistent patterns

### Dependency Injection
**Standard**: Use `Depends()` for all external dependencies

```python
# ✅ STANDARD: Dependency injection (from inline_annotations.py)
@router.post("", response_model=InlineAnnotationResponse, status_code=status.HTTP_201_CREATED)
def create_inline_annotation(
    annotation: InlineAnnotationCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> InlineAnnotationResponse:
    try:
        new_annotation = create_annotation(
            entity_type=annotation.entity_type,
            entity_id=annotation.entity_id,
            content=annotation.content,
            user_id=current_user.id,
            db=db,
        )
        return InlineAnnotationResponse.model_validate(new_annotation)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
```

### Schema Usage
**Standard**: Always use Pydantic request/response models

```python
# ✅ STANDARD: Pydantic schemas (from programs.py)
@router.post("/", response_model=Program, status_code=201)
async def create_program(
    program: ProgramCreate,  # Request schema
    db: Session = Depends(get_database_session),
) -> Program:  # Response schema
    return program_service.create_program(db, program)
```

### Error Handling
**Standard**: HTTPException with standardized status codes

```python
# ✅ STANDARD: Error handling (from sharing.py)
@router.post("/entities/{entity_type}/{entity_id}/share")
def share_entity_endpoint(
    entity_type: str,
    entity_id: UUID,
    share_request: EntityShareCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityShareResponse:
    # Validate entity type
    valid_types = ["dataset", "experiment", "compound", "signature"]
    if entity_type not in valid_types:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type. Must be one of: {', '.join(valid_types)}"
        )
    
    # Check authorization
    if not can_share_entity(str(current_user.id), entity_type, str(entity_id), db):
        raise HTTPException(status_code=403, detail="Not authorized to share this entity")
    
    try:
        share = share_entity(entity_type=entity_type, entity_id=entity_id, ...)
        return EntityShareResponse.model_validate(share)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to create share: {e}")
```

---

## Database Layer (`amprenta_rag/database/`)

**Pattern**: Repository pattern with domain-specific CRUD modules  
**Files**: Base session manager, ~50 models, modular CRUD

### Model Definitions
**Standard**: UUID primary keys with relationship() and back_populates

```python
# ✅ STANDARD: Model definition (from models.py)
class Program(Base):
    __tablename__ = "programs"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    description = Column(Text)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc))
    
    # Relationships with back_populates
    experiments: Mapped[List["Experiment"]] = relationship(
        "Experiment", 
        secondary=program_experiment_assoc,
        back_populates="programs"
    )
    
    datasets: Mapped[List["Dataset"]] = relationship(
        "Dataset",
        secondary=program_dataset_assoc,
        back_populates="programs"
    )
```

### Session Management
**Standard**: Context manager with automatic rollback

```python
# ✅ STANDARD: Session context manager (from session.py)
@contextmanager
def db_session() -> Generator[Session, None, None]:
    """Context manager for database sessions with automatic cleanup."""
    db_gen = None
    db = None
    
    try:
        db_gen = get_db()
        db = next(db_gen)
        
        yield db
        db.commit()
        
    except Exception as e:
        if db:
            try:
                db.rollback()
            except Exception:
                pass
        raise e
        
    finally:
        if db:
            try:
                db.close()
            except Exception:
                pass
```

### CRUD Operations
**Standard**: Repository pattern in domain-specific modules

```python
# ✅ STANDARD: Repository pattern (from crud_datasets.py)
def create_dataset(
    db: Session,
    name: str,
    omics_type: str,
    description: Optional[str] = None,
    file_paths: Optional[List[str]] = None,
) -> Dataset:
    """Create a new dataset with validation."""
    dataset = Dataset(
        id=uuid4(),
        name=name,
        omics_type=omics_type,
        description=description,
        file_paths=file_paths or [],
    )
    
    db.add(dataset)
    db.commit()
    db.refresh(dataset)
    
    logger.info(f"Created dataset: {name}")
    return dataset
```

---

## Sync/Adapters (`amprenta_rag/sync/adapters/`)

**Pattern**: BaseSyncAdapter ABC with 3 abstract methods  
**Files**: 1 base interface, 5 adapter implementations

### Base Interface
**Standard**: Abstract base class with consistent method signatures

```python
# ✅ STANDARD: Base adapter (from base.py)
class BaseSyncAdapter(ABC):
    """Base adapter for synchronizing records from an external data source."""
    
    source: str  # "chembl", "pubchem", etc.
    
    @abstractmethod
    async def fetch_records(self, since: datetime | None) -> AsyncIterator[dict]:
        """Fetch records from external source, optionally since last sync."""
    
    @abstractmethod
    def compute_checksum(self, record: dict) -> str:
        """Compute MD5 hash for change detection."""
    
    @abstractmethod
    def map_to_entity(self, record: dict, db_session) -> tuple[str, UUID | None]:
        """Map external record to local entity. Returns (entity_type, entity_id)."""
```

### Implementation Pattern
**Standard**: Inherit from BaseSyncAdapter with async fetch_records

```python
# ✅ STANDARD: Adapter implementation (from uniprot_mapping.py)
class UniProtMappingAdapter(BaseSyncAdapter):
    source = "uniprot"
    
    async def fetch_records(self, since: datetime | None) -> AsyncIterator[dict]:
        """Download latest UniProt ID mappings."""
        url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
        
        async with httpx.AsyncClient() as client:
            async with client.stream("GET", url) as response:
                response.raise_for_status()
                
                async for line in self._parse_mapping_file(response):
                    yield {
                        "uniprot_id": line[0],
                        "db_name": line[1], 
                        "db_id": line[2],
                    }
    
    def compute_checksum(self, record: dict) -> str:
        """Compute MD5 hash for change detection."""
        content = f"{record['uniprot_id']}|{record['db_name']}|{record['db_id']}"
        return hashlib.md5(content.encode()).hexdigest()
    
    def map_to_entity(self, record: dict, db_session) -> tuple[str, UUID | None]:
        """Map to ID mapping entity."""
        return ("id_mapping", None)
```

---

## ML/Analysis Patterns (`amprenta_rag/ml/`)

**Pattern**: Class-based predictors with registry pattern  
**Files**: 2 main predictor classes, ensemble patterns

### Predictor Classes
**Standard**: Class-based with `__init__`, `predict()`, registry integration

```python
# ✅ STANDARD: Predictor pattern (from admet/predictor.py)
class ADMETPredictor:
    """ADMET prediction service using registered models."""
    
    def __init__(self):
        self.registry = get_registry()
    
    def predict(
        self,
        smiles_list: List[str],
        endpoints: Optional[List[str]] = None,
        include_shap: bool = False,
    ) -> List[Dict[str, Any]]:
        """Predict ADMET properties for compounds."""
        endpoints = endpoints or list(ADMET_MODELS.keys())
        results = []
        
        for smiles in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    results.append({"smiles": smiles, "error": "Invalid SMILES"})
                    continue
                
                pred = {"smiles": smiles}
                
                for endpoint in endpoints:
                    model_name = ADMET_MODELS.get(endpoint)
                    if not model_name:
                        continue
                        
                    try:
                        ml_model = self.registry.get_active_model(model_name)
                        # Prediction logic here
                        pred[endpoint] = prediction_value
                    except Exception as e:
                        pred[endpoint] = {"error": str(e)}
                
                results.append(pred)
                
            except Exception as e:
                results.append({"smiles": smiles, "error": str(e)})
        
        return results
```

### Ensemble Pattern
**Standard**: Ensemble + calibration + applicability

```python
# ✅ STANDARD: Ensemble pattern (from qsar/trainer.py)
def train_target_model(
    target: str,
    n_models: int = 5,
    calibration_method: str = "isotonic",
) -> Dict[str, Any]:
    """Train ensemble model with calibration and applicability."""
    
    # Train ensemble
    ensemble = BootstrapEnsemble(n_models=n_models).fit(X_train, y_train)
    
    # Calibrate on calibration split
    cal_mean, _cal_std = ensemble.predict_proba(X_cal)
    calibrator = CalibrationWrapper(method=calibration_method).fit(cal_mean, y_cal)
    
    # Applicability checker
    app = ApplicabilityChecker(threshold=0.3).fit(X_train)
    
    # Create artifact
    artifact = {
        "ensemble": ensemble.to_artifact(),
        "calibrator": calibrator,
        "applicability": {"threshold": app.threshold, "centroid": app.training_centroid},
        "metadata": {"target": target, "train_size": X_train.shape[0]}
    }
    
    return artifact
```

---

## Dashboard Patterns (`scripts/dashboard/`)

**Pattern**: `render_*_page()` naming with tab-based layouts  
**Files**: ~30 pages with consistent structure

### Page Structure
**Standard**: `render_*_page()` naming convention

```python
# ✅ STANDARD: Page structure (from chemistry/main.py)
def render_chemistry_page() -> None:
    """Render the Chemistry page with modular tab handlers."""
    st.header("⚗️ Chemistry")
    
    tab_compounds, tab_reg, tab_struct = st.tabs([
        "Compounds",
        "Register Compound", 
        "Structure Search"
    ])
    
    render_compounds_tab(
        tab_compounds,
        db_session=db_session,
        compound_model=Compound,
    )
    
    render_structure_tab(
        tab_struct,
        db_session=db_session,
        substructure_search=substructure_search,
    )
```

### Import Pattern
**Standard**: Direct module imports to avoid Streamlit caching issues

```python
# ✅ STANDARD: Import pattern (from IMPORT_PATTERN.md)
# Import the models module directly
import amprenta_rag.database.models as db_models

# Access models via module attribute
Compound = db_models.Compound
HTSCampaign = db_models.HTSCampaign

# ❌ AVOID: Direct imports (can fail with Streamlit caching)
# from amprenta_rag.database.models import Compound, HTSCampaign
```

### API Helpers
**Standard**: httpx client with error handling

```python
# ✅ STANDARD: API helper pattern (from docking_runs.py)
API_BASE = os.environ.get("API_URL", "http://localhost:8000")

def _api_get(endpoint: str, params: Optional[Dict] = None) -> Optional[Dict]:
    """Generic API GET with error handling."""
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{API_BASE}{endpoint}", params=params)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {e}")
        return None

# Usage in page
def render_docking_runs_page() -> None:
    structures = _api_get("/api/structures", {"limit": 200})
    if not structures:
        st.info("No structures available")
        return
```

---

## Background Jobs (`amprenta_rag/jobs/`)

**Pattern**: Celery tasks with consistent retry and queue patterns  
**Files**: 11 task modules with standardized structure

### Task Definition
**Standard**: `@celery_app.task(bind=True)` with retry configuration

```python
# ✅ STANDARD: Task pattern (from tasks/mapping_refresh.py)
@celery_app.task(bind=True, max_retries=3, default_retry_delay=300, queue='scheduled')
def refresh_uniprot_mappings_task(self) -> Dict[str, any]:
    """Weekly UniProt ID mapping refresh (bulk download)."""
    try:
        logger.info("Starting UniProt mappings refresh task")
        
        from amprenta_rag.services.id_mapping_service import refresh_uniprot_mappings
        count = asyncio.run(refresh_uniprot_mappings())
        
        logger.info(f"UniProt refresh completed: {count} mappings updated")
        return {"status": "success", "count": count}
        
    except Exception as exc:
        logger.exception("UniProt refresh failed")
        
        # Retry with exponential backoff
        if self.request.retries < self.max_retries:
            countdown = 300 * (2 ** self.request.retries)  # 5, 10, 20 minutes
            raise self.retry(exc=exc, countdown=countdown)
        
        return {"status": "failed", "error": str(exc)}
```

### Error Handling with Database Updates
**Standard**: Update job status in database, graceful retry

```python
# ✅ STANDARD: Error handling (from tasks/extraction.py)
@celery_app.task(bind=True, max_retries=3, default_retry_delay=60, queue='high')
def process_extraction_job(self, job_id: str) -> dict:
    """Process batch document extraction."""
    job_uuid = UUID(job_id)
    
    try:
        svc = ExtractionBatchService(db_session)
        job = svc.process_job(job_uuid)
        return {"status": job.status, "job_id": job_id}
    
    except Exception as exc:
        # Update job status to failed if possible
        try:
            with db_session() as db:
                job = db.query(ExtractionJob).filter(ExtractionJob.id == job_uuid).first()
                if job is not None:
                    job.status = "failed"
                    job.completed_at = datetime.now(timezone.utc)
                    db.commit()
        except Exception:
            pass  # Don't let DB errors prevent retry logic
        
        # Retry with exponential backoff
        if self.request.retries >= self.max_retries:
            return {"status": "failed", "error": str(exc), "job_id": job_id}
        
        self.retry(exc=exc, countdown=60 * (2 ** self.request.retries))
```

### Queue Configuration
**Standard**: Task-specific queue routing

```python
# ✅ STANDARD: Queue configuration (from config.py)
task_routes = {
    'amprenta_rag.jobs.tasks.genomics.*': {'queue': 'default'},
    'amprenta_rag.jobs.tasks.docking.*': {'queue': 'high'},
    'amprenta_rag.jobs.tasks.extraction.*': {'queue': 'high'},
    'amprenta_rag.jobs.tasks.sync.*': {'queue': 'low'},
    'amprenta_rag.jobs.tasks.single_cell.*': {'queue': 'default'},
}

task_queues = (
    Queue('default', routing_key='default'),
    Queue('high', routing_key='high'),
    Queue('low', routing_key='low'),
    Queue('scheduled', routing_key='scheduled'),
)
```

---

## Testing Patterns

### Unit Tests
**Standard**: Mock external dependencies, test business logic

```python
# ✅ STANDARD: Unit test pattern
@patch('amprenta_rag.services.activity.db_session')
def test_log_activity_success(mock_db_session):
    """Test activity logging with mocked database."""
    mock_db = MagicMock()
    mock_db_session.return_value.__enter__.return_value = mock_db
    
    result = log_activity(
        event_type=ActivityEventType.COMPOUND_ADDED,
        target_type="compound",
        target_id=uuid4(),
        target_name="Test Compound"
    )
    
    assert result is not None
    mock_db.add.assert_called_once()
    mock_db.commit.assert_called_once()
```

### Integration Tests
**Standard**: Real database with transaction rollback

```python
# ✅ STANDARD: Integration test pattern
@pytest.mark.integration
def test_create_program_integration(integration_client, db_session):
    """Test program creation with real database."""
    program_data = {
        "name": f"Test Program {uuid4().hex[:8]}",
        "description": "Integration test program"
    }
    
    response = integration_client.post("/api/v1/programs", json=program_data)
    
    assert response.status_code == 201
    data = response.json()
    
    # Verify in real database
    program_id = UUID(data["id"])
    db_program = db_session.query(Program).filter(Program.id == program_id).first()
    assert db_program.name == program_data["name"]
```

### E2E Tests
**Standard**: Playwright with semantic selectors

```python
# ✅ STANDARD: E2E test pattern
def test_create_annotation_workflow(page, streamlit_server):
    """Test full annotation creation workflow."""
    page.goto(f"{streamlit_server}/Lab_Notebook")
    
    # Use semantic selectors
    page.get_by_role("button", name="View Annotations").click()
    page.get_by_text("Create New Annotation").click()
    
    # Fill form
    page.fill('textarea[aria-label="Annotation Content"]', "Test annotation")
    page.get_by_role("button", name="Create Annotation").click()
    
    # Verify result
    expect(page.get_by_text("Annotation created successfully")).to_be_visible()
```

---

## Anti-patterns to Avoid

### ❌ Mixed Session Management
```python
# DON'T: Mix session patterns in same layer
def bad_service_function():
    with db_session() as db:  # Context manager
        # ... some logic
        other_function(db)    # Passing session around
        
def other_function(db: Session):  # Dependency injection
    # This creates confusion about session ownership
```

### ❌ Hardcoded URLs
```python
# DON'T: Hardcode API URLs
def bad_api_call():
    response = requests.get("http://localhost:8000/api/data")  # Hardcoded!

# DO: Use configuration
API_BASE = os.environ.get("API_URL", "http://localhost:8000")
response = requests.get(f"{API_BASE}/api/data")
```

### ❌ Skipping Tests
```python
# DON'T: Skip tests (forbidden in this codebase)
@pytest.mark.skip("TODO: fix later")  # NEVER do this
def test_important_feature():
    pass

# DO: Fix the test or delete it
def test_important_feature():
    # Proper test implementation
    assert actual_behavior() == expected_result
```

### ❌ Implicit Dependencies
```python
# DON'T: Use global state or implicit dependencies
current_user = None  # Global state

@router.get("/data")
def get_data():
    if current_user is None:  # Implicit dependency
        raise HTTPException(401)

# DO: Use explicit dependency injection
@router.get("/data")
def get_data(current_user: User = Depends(get_current_user)):
    # Explicit dependency
```

---

## Summary

This architecture guide documents the established patterns across the Amprenta platform. Following these patterns ensures:

- **Consistency** across the codebase
- **Maintainability** for future development
- **Testability** with clear separation of concerns
- **Reliability** through proven error handling patterns

When implementing new features, refer to the relevant section and follow the documented patterns. When in doubt, examine existing code in the referenced files for concrete examples.

For questions about architectural decisions or pattern exceptions, refer to the development team or create an architecture decision record (ADR).
