"""Active learning service for iterative model improvement."""
import logging
from uuid import UUID
from typing import List, Optional, Dict, Any
from datetime import datetime, timezone
import numpy as np

from sqlalchemy.orm import Session
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import LabelQueueItem, ActiveLearningCycle, MLModel, Compound

logger = logging.getLogger(__name__)


# ============ Sampling Strategies ============

def uncertainty_sampling(uncertainties: np.ndarray) -> np.ndarray:
    """Rank by highest uncertainty (std deviation)."""
    return np.argsort(-uncertainties)


def margin_sampling(predictions: np.ndarray) -> np.ndarray:
    """Rank by smallest margin from decision boundary (0.5)."""
    margins = np.abs(predictions - 0.5)
    return np.argsort(margins)


def entropy_sampling(predictions: np.ndarray) -> np.ndarray:
    """Rank by highest prediction entropy."""
    p = np.clip(predictions, 1e-10, 1 - 1e-10)
    entropy = -p * np.log(p) - (1 - p) * np.log(1 - p)
    return np.argsort(-entropy)


def hybrid_sampling(
    predictions: np.ndarray,
    uncertainties: np.ndarray,
    applicability: np.ndarray,
    weights: Dict[str, float] = None
) -> np.ndarray:
    """Combined ranking with configurable weights."""
    if weights is None:
        weights = {"uncertainty": 0.5, "margin": 0.3, "applicability": 0.2}
    
    # Normalize each score to [0, 1]
    def normalize(arr):
        if arr.max() == arr.min():
            return np.zeros_like(arr)
        return (arr - arr.min()) / (arr.max() - arr.min())
    
    unc_score = normalize(uncertainties)
    margin_score = 1 - normalize(np.abs(predictions - 0.5))  # Invert: smaller margin = higher score
    app_score = normalize(applicability)  # Higher distance = more novel
    
    combined = (
        weights["uncertainty"] * unc_score +
        weights["margin"] * margin_score +
        weights["applicability"] * app_score
    )
    return np.argsort(-combined)


# ============ Service Class ============

class ActiveLearningService:
    """Service for active learning operations."""
    
    def __init__(self, model_id: UUID, db: Optional[Session] = None):
        self.model_id = model_id
        self._db = db
    
    def select_samples(
        self,
        compound_ids: List[UUID],
        strategy: str = "uncertainty",
        batch_size: int = 50,
        cycle_number: int = 1
    ) -> List[LabelQueueItem]:
        """Select most informative samples for labeling.
        
        Args:
            compound_ids: Pool of compound IDs to consider
            strategy: Selection strategy (uncertainty, margin, entropy, hybrid, random)
            batch_size: Number of samples to select
            cycle_number: Current AL cycle number
            
        Returns:
            List of created LabelQueueItem objects
        """
        with db_session() as db:
            # Get model and load ensemble for predictions
            model = db.query(MLModel).filter(MLModel.id == self.model_id).first()
            if not model:
                raise ValueError(f"Model {self.model_id} not found")
            
            # Get compounds not already in queue for this model
            existing = db.query(LabelQueueItem.compound_id).filter(
                LabelQueueItem.model_id == self.model_id,
                LabelQueueItem.status.in_(["pending", "in_progress"])
            ).all()
            existing_ids = {e[0] for e in existing}
            
            pool_ids = [cid for cid in compound_ids if cid not in existing_ids]
            if not pool_ids:
                logger.warning("No unlabeled compounds in pool")
                return []
            
            # Load compounds and get SMILES
            compounds = db.query(Compound).filter(Compound.id.in_(pool_ids)).all()
            smiles_list = [c.smiles for c in compounds if c.smiles]
            compound_map = {c.smiles: c.id for c in compounds if c.smiles}
            
            if not smiles_list:
                return []
            
            # Get predictions and uncertainties from model
            # (In real implementation, load model artifact and run inference)
            # For MVP, generate mock uncertainty scores
            n = len(smiles_list)
            predictions = np.random.uniform(0.3, 0.7, n)  # Mock predictions near boundary
            uncertainties = np.random.uniform(0.1, 0.5, n)  # Mock uncertainties
            applicability = np.random.uniform(0, 1, n)  # Mock applicability distances
            
            # Apply selection strategy
            if strategy == "uncertainty":
                ranking = uncertainty_sampling(uncertainties)
            elif strategy == "margin":
                ranking = margin_sampling(predictions)
            elif strategy == "entropy":
                ranking = entropy_sampling(predictions)
            elif strategy == "hybrid":
                ranking = hybrid_sampling(predictions, uncertainties, applicability)
            elif strategy == "random":
                ranking = np.random.permutation(n)
            else:
                raise ValueError(f"Unknown strategy: {strategy}")
            
            # Select top batch_size
            selected_indices = ranking[:batch_size]
            
            # Create queue items
            items = []
            for rank, idx in enumerate(selected_indices):
                smiles = smiles_list[idx]
                compound_id = compound_map[smiles]
                
                item = LabelQueueItem(
                    compound_id=compound_id,
                    model_id=self.model_id,
                    prediction=float(predictions[idx]),
                    uncertainty=float(uncertainties[idx]),
                    applicability_distance=float(applicability[idx]),
                    selection_strategy=strategy,
                    selection_batch=cycle_number,
                    priority_score=float(batch_size - rank),  # Higher rank = higher priority
                    status="pending"
                )
                db.add(item)
                items.append(item)
            
            db.commit()
            for item in items:
                db.expunge(item)
            
            logger.info(f"Selected {len(items)} samples using {strategy} strategy")
            return items
    
    def get_pending_items(
        self,
        status: str = "pending",
        limit: int = 50
    ) -> List[LabelQueueItem]:
        """Get items from queue by status."""
        with db_session() as db:
            items = db.query(LabelQueueItem).filter(
                LabelQueueItem.model_id == self.model_id,
                LabelQueueItem.status == status
            ).order_by(LabelQueueItem.priority_score.desc()).limit(limit).all()
            
            for item in items:
                db.expunge(item)
            return items
    
    def submit_label(
        self,
        item_id: UUID,
        label: float,
        source: str,
        confidence: str,
        user_id: UUID,
        notes: Optional[str] = None
    ) -> LabelQueueItem:
        """Submit a label for a queue item."""
        with db_session() as db:
            item = db.query(LabelQueueItem).filter(LabelQueueItem.id == item_id).first()
            if not item:
                raise ValueError(f"Queue item {item_id} not found")
            
            item.label = label
            item.label_source = source
            item.label_confidence = confidence
            item.labeled_by = user_id
            item.labeled_at = datetime.now(timezone.utc)
            item.notes = notes
            item.status = "labeled"
            
            db.commit()
            db.expunge(item)
            
            logger.info(f"Label submitted for item {item_id}")
            return item
    
    def skip_item(self, item_id: UUID, user_id: UUID) -> LabelQueueItem:
        """Mark an item as skipped."""
        with db_session() as db:
            item = db.query(LabelQueueItem).filter(LabelQueueItem.id == item_id).first()
            if not item:
                raise ValueError(f"Queue item {item_id} not found")
            
            item.status = "skipped"
            item.labeled_by = user_id
            item.labeled_at = datetime.now(timezone.utc)
            
            db.commit()
            db.expunge(item)
            return item
    
    def create_cycle(
        self,
        selection_strategy: str,
        batch_size: int
    ) -> ActiveLearningCycle:
        """Create a new active learning cycle."""
        with db_session() as db:
            # Get next cycle number
            last_cycle = db.query(ActiveLearningCycle).filter(
                ActiveLearningCycle.model_id == self.model_id
            ).order_by(ActiveLearningCycle.cycle_number.desc()).first()
            
            cycle_number = (last_cycle.cycle_number + 1) if last_cycle else 1
            
            # Get current model metrics
            model = db.query(MLModel).filter(MLModel.id == self.model_id).first()
            metrics_before = model.metrics if model else {}
            
            cycle = ActiveLearningCycle(
                model_id=self.model_id,
                cycle_number=cycle_number,
                selection_strategy=selection_strategy,
                batch_size=batch_size,
                status="selecting",
                metrics_before=metrics_before
            )
            db.add(cycle)
            db.commit()
            db.expunge(cycle)
            
            logger.info(f"Created AL cycle {cycle_number} for model {self.model_id}")
            return cycle
    
    def get_cycle_stats(self) -> Dict[str, Any]:
        """Get active learning statistics for the model."""
        with db_session() as db:
            cycles = db.query(ActiveLearningCycle).filter(
                ActiveLearningCycle.model_id == self.model_id
            ).order_by(ActiveLearningCycle.cycle_number).all()
            
            pending = db.query(LabelQueueItem).filter(
                LabelQueueItem.model_id == self.model_id,
                LabelQueueItem.status == "pending"
            ).count()
            
            labeled = db.query(LabelQueueItem).filter(
                LabelQueueItem.model_id == self.model_id,
                LabelQueueItem.status == "labeled"
            ).count()
            
            return {
                "model_id": str(self.model_id),
                "total_cycles": len(cycles),
                "pending_items": pending,
                "labeled_items": labeled,
                "cycles": [
                    {
                        "cycle_number": c.cycle_number,
                        "status": c.status,
                        "items_labeled": c.items_labeled,
                        "metrics_before": c.metrics_before,
                        "metrics_after": c.metrics_after
                    }
                    for c in cycles
                ]
            }
    
    def update_cycle_stats(self, cycle_id: UUID) -> Optional[ActiveLearningCycle]:
        """Update cycle statistics based on current queue state."""
        with db_session() as db:
            cycle = db.query(ActiveLearningCycle).filter(ActiveLearningCycle.id == cycle_id).first()
            if not cycle:
                return None
            
            # Count items in this cycle
            labeled = db.query(LabelQueueItem).filter(
                LabelQueueItem.model_id == cycle.model_id,
                LabelQueueItem.selection_batch == cycle.cycle_number,
                LabelQueueItem.status == "labeled"
            ).count()
            
            skipped = db.query(LabelQueueItem).filter(
                LabelQueueItem.model_id == cycle.model_id,
                LabelQueueItem.selection_batch == cycle.cycle_number,
                LabelQueueItem.status == "skipped"
            ).count()
            
            cycle.items_labeled = labeled
            cycle.items_skipped = skipped
            
            # Update status if all items are complete
            if labeled + skipped >= cycle.batch_size:
                cycle.status = "labeling_complete"
            
            db.commit()
            db.expunge(cycle)
            return cycle
    
    def get_labeling_queue(
        self,
        assigned_to: Optional[UUID] = None,
        status: str = "pending",
        limit: int = 20
    ) -> List[Dict[str, Any]]:
        """Get labeling queue with compound details."""
        with db_session() as db:
            query = db.query(LabelQueueItem).join(Compound).filter(
                LabelQueueItem.model_id == self.model_id,
                LabelQueueItem.status == status
            )
            
            if assigned_to:
                query = query.filter(LabelQueueItem.assigned_to == assigned_to)
            
            items = query.order_by(LabelQueueItem.priority_score.desc()).limit(limit).all()
            
            result = []
            for item in items:
                result.append({
                    "id": item.id,
                    "compound_id": item.compound_id,
                    "smiles": item.compound.smiles,
                    "compound_name": item.compound.compound_id,
                    "prediction": item.prediction,
                    "uncertainty": item.uncertainty,
                    "applicability_distance": item.applicability_distance,
                    "priority_score": item.priority_score,
                    "selection_strategy": item.selection_strategy,
                    "selection_batch": item.selection_batch,
                    "status": item.status,
                    "assigned_to": item.assigned_to,
                    "created_at": item.created_at
                })
            
            return result
