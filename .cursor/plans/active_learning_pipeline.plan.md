# Active Learning Pipeline

## Overview
Human-in-the-loop ML model refinement using uncertainty sampling. Scientists label high-uncertainty compounds to iteratively improve QSAR/ADMET models.

**Timeline:** 2-3 weeks
**Dependencies:** Existing BootstrapEnsemble, MLModelRegistry, MLModel table

---

## Architecture

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│  Model          │───▶│  Uncertainty     │───▶│  Label Queue    │
│  (Ensemble)     │    │  Sampling        │    │  (DB table)     │
└─────────────────┘    └──────────────────┘    └─────────────────┘
                                                       │
                                                       ▼
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│  Retrain        │◀───│  Labeling UI     │◀───│  Scientist      │
│  Pipeline       │    │  (Dashboard)     │    │  Review         │
└─────────────────┘    └──────────────────┘    └─────────────────┘
```

---

## Data Model

### LabelQueue Table
```python
class LabelQueueItem(Base):
    __tablename__ = "label_queue"
    
    id = Column(UUID, primary_key=True)
    
    # What to label
    compound_id = Column(UUID, ForeignKey("compounds.id"), nullable=False)
    model_id = Column(UUID, ForeignKey("ml_models.id"), nullable=False)
    
    # Uncertainty info
    prediction = Column(Float)  # Model's current prediction
    uncertainty = Column(Float)  # Std deviation or entropy
    applicability_distance = Column(Float)  # Distance from training domain
    
    # Selection metadata
    selection_strategy = Column(String(50))  # uncertainty, margin, entropy, random
    selection_batch = Column(Integer)  # Which iteration added this
    priority_score = Column(Float)  # Combined ranking score
    
    # Labeling
    status = Column(String(50), default="pending")  # pending, in_progress, labeled, skipped
    assigned_to = Column(UUID, ForeignKey("users.id"), nullable=True)
    assigned_at = Column(DateTime, nullable=True)
    
    # Label result
    label = Column(Float, nullable=True)  # The human-provided label
    label_source = Column(String(100), nullable=True)  # "experimental", "literature", "expert_estimate"
    label_confidence = Column(String(20), nullable=True)  # high, medium, low
    labeled_at = Column(DateTime, nullable=True)
    labeled_by = Column(UUID, ForeignKey("users.id"), nullable=True)
    notes = Column(Text, nullable=True)
    
    created_at = Column(DateTime)
    updated_at = Column(DateTime)
    
    # P1 FIX: Unique constraint to prevent duplicate selections
    __table_args__ = (
        UniqueConstraint("compound_id", "model_id", "selection_batch", name="uq_label_queue_compound_model_batch"),
        Index("ix_label_queue_status_model", "status", "model_id"),
    )
```

### ActiveLearningCycle Table
```python
class ActiveLearningCycle(Base):
    __tablename__ = "active_learning_cycles"
    
    id = Column(UUID, primary_key=True)
    model_id = Column(UUID, ForeignKey("ml_models.id"), nullable=False)
    
    # Cycle info
    cycle_number = Column(Integer, nullable=False)
    selection_strategy = Column(String(50))
    batch_size = Column(Integer)
    
    # Status
    status = Column(String(50))  # selecting, labeling, training, complete
    
    # Metrics before/after
    metrics_before = Column(JSON)  # {"auc": 0.85, "ece": 0.08}
    metrics_after = Column(JSON)   # {"auc": 0.88, "ece": 0.06}
    
    # Stats
    items_selected = Column(Integer)
    items_labeled = Column(Integer)
    items_skipped = Column(Integer)
    
    started_at = Column(DateTime)
    completed_at = Column(DateTime, nullable=True)
```

---

## Service Layer

### `amprenta_rag/ml/active_learning.py`

```python
class ActiveLearningService:
    """Active learning for iterative model improvement."""
    
    def __init__(self, model_id: UUID):
        self.model_id = model_id
        self.model = load_model(model_id)
    
    def select_samples(
        self,
        pool: List[str],  # SMILES to consider
        strategy: str = "uncertainty",
        batch_size: int = 50,
        exclude_labeled: bool = True
    ) -> List[dict]:
        """Select most informative samples for labeling.
        
        Strategies:
        - uncertainty: Highest prediction variance
        - margin: Smallest margin between classes
        - entropy: Highest prediction entropy
        - random: Random baseline
        - hybrid: Weighted combination
        """
        
    def add_to_queue(
        self,
        compound_ids: List[UUID],
        selection_batch: int
    ) -> int:
        """Add selected compounds to labeling queue."""
        
    def get_pending_items(
        self,
        user_id: Optional[UUID] = None,
        limit: int = 20
    ) -> List[LabelQueueItem]:
        """Get items pending labeling."""
        
    def submit_label(
        self,
        item_id: UUID,
        label: float,
        source: str,
        confidence: str,
        notes: Optional[str],
        user_id: UUID
    ) -> LabelQueueItem:
        """Submit a label for a queue item."""
        
    def trigger_retrain(
        self,
        cycle_id: UUID,
        min_new_labels: int = 10
    ) -> Optional[UUID]:
        """Retrain model with newly labeled data.
        
        P1 FIX - Validation Handling:
        - Uses SAME held-out validation set from original training
        - New labels added to training set only
        - Before/after metrics computed on identical validation set
        
        P1 FIX - Queue Item Handling After Retrain:
        - Items with status="labeled" → archived (linked to old model version)
        - Items with status="pending" → re-scored with new model uncertainty
        - New model version created, old model marked as "superseded"
        """
        
    def get_cycle_stats(self, model_id: UUID) -> dict:
        """Get active learning cycle statistics."""
```

### Uncertainty Sampling Strategies

```python
def uncertainty_sampling(predictions: np.ndarray, uncertainties: np.ndarray) -> np.ndarray:
    """Rank by highest uncertainty (std deviation)."""
    return np.argsort(-uncertainties)

def margin_sampling(predictions: np.ndarray) -> np.ndarray:
    """Rank by smallest margin from decision boundary (0.5)."""
    margins = np.abs(predictions - 0.5)
    return np.argsort(margins)

def entropy_sampling(predictions: np.ndarray) -> np.ndarray:
    """Rank by highest prediction entropy."""
    p = np.clip(predictions, 1e-10, 1-1e-10)
    entropy = -p * np.log(p) - (1-p) * np.log(1-p)
    return np.argsort(-entropy)

def hybrid_sampling(
    predictions: np.ndarray,
    uncertainties: np.ndarray,
    applicability: np.ndarray,
    weights: dict = {"uncertainty": 0.5, "margin": 0.3, "applicability": 0.2}
) -> np.ndarray:
    """Combined ranking with configurable weights."""
```

---

## API Endpoints

### `amprenta_rag/api/routers/active_learning.py`

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | /active-learning/select | Select samples for labeling |
| GET | /active-learning/queue | Get labeling queue |
| GET | /active-learning/queue/{item_id} | Get queue item details |
| POST | /active-learning/queue/{item_id}/label | Submit label |
| POST | /active-learning/queue/{item_id}/skip | Skip item |
| POST | /active-learning/queue/{item_id}/assign | Assign to user |
| GET | /active-learning/cycles | List AL cycles for model |
| POST | /active-learning/cycles/{cycle_id}/retrain | Trigger retraining |
| GET | /active-learning/stats/{model_id} | Get AL statistics |

---

## Dashboard

### `scripts/dashboard/pages/active_learning.py`

**Tab 1: Label Queue**
- Table of pending items with: Compound, Structure (2D), Prediction, Uncertainty, Priority
- Click to open labeling modal
- Filters: status, assigned_to, model
- Bulk actions: assign, skip

**Tab 2: Labeling Interface**
- Selected compound details (structure, properties, similar compounds)
- Model's prediction with uncertainty visualization
- Label input form:
  - Value (numeric or class)
  - Source dropdown (experimental, literature, expert)
  - Confidence (high/medium/low)
  - Notes field
- Navigation: Next, Skip, Save & Next

**Tab 3: Cycle History**
- Timeline of AL cycles per model
- Before/after metrics comparison
- Items labeled per cycle
- Model improvement visualization (learning curve)

**Tab 4: Sample Selection**
- Select model from dropdown
- Choose pool (all unlabeled, specific dataset, compound list)
- Configure strategy and batch size
- Preview selected samples
- "Add to Queue" button

---

## Implementation Batches

### Batch 1: Data Model + Service (Day 1-3)
- [ ] Create LabelQueueItem model
- [ ] Create ActiveLearningCycle model
- [ ] Generate Alembic migration
- [ ] Implement ActiveLearningService core methods
- [ ] Implement sampling strategies

### Batch 2: API Endpoints (Day 4-5)
- [ ] Create active_learning router
- [ ] Implement all 9 endpoints
- [ ] Add Pydantic schemas
- [ ] Register router

### Batch 3: Dashboard (Day 6-8)
- [ ] Tab 1: Label Queue table
- [ ] Tab 2: Labeling Interface with compound viz
- [ ] Tab 3: Cycle History with charts
- [ ] Tab 4: Sample Selection wizard

### Batch 4: Integration + Tests (Day 9-10)
- [ ] Connect to existing BootstrapEnsemble
- [ ] Implement retrain trigger with Celery
- [ ] Service tests (12)
- [ ] API tests (12)
- [ ] E2E tests (6)
- [ ] Documentation

---

## Success Criteria

- [ ] Scientists can view and label high-uncertainty compounds
- [ ] System tracks labeling cycles with before/after metrics
- [ ] Model retraining works with new labels
- [ ] 4 sampling strategies implemented
- [ ] 30+ tests passing
- [ ] Dashboard fully functional

---

## TODOs

- [ ] Batch 1: Data Model + Service
- [ ] Batch 2: API Endpoints
- [ ] Batch 3: Dashboard (4 tabs)
- [ ] Batch 4: Integration + Tests


