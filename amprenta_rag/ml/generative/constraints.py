"""Property constraints for generative chemistry optimization."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class PropertyConstraint:
    """Single property constraint for molecular optimization."""
    
    name: str  # e.g., "logP", "herg", "mw", "tpsa"
    min_value: Optional[float] = None
    max_value: Optional[float] = None
    target_value: Optional[float] = None  # For optimization toward specific target
    weight: float = 1.0  # Importance weight for penalty calculation
    
    def __post_init__(self):
        """Validate constraint parameters."""
        if self.min_value is not None and self.max_value is not None:
            if self.min_value > self.max_value:
                raise ValueError(f"min_value ({self.min_value}) > max_value ({self.max_value})")
        
        if self.weight <= 0:
            raise ValueError(f"weight must be positive, got {self.weight}")
    
    def is_satisfied(self, value: float) -> bool:
        """Check if a value satisfies this constraint.
        
        Args:
            value: Property value to check
            
        Returns:
            True if constraint is satisfied
        """
        if self.min_value is not None and value < self.min_value:
            return False
        
        if self.max_value is not None and value > self.max_value:
            return False
        
        return True
    
    def compute_penalty(self, value: float) -> float:
        """Compute penalty for constraint violation.
        
        Args:
            value: Property value
            
        Returns:
            Penalty value (0 if satisfied, positive if violated)
        """
        penalty = 0.0
        
        # Range violations
        if self.min_value is not None and value < self.min_value:
            penalty += (self.min_value - value) ** 2
        
        if self.max_value is not None and value > self.max_value:
            penalty += (value - self.max_value) ** 2
        
        # Target optimization (if specified)
        if self.target_value is not None:
            penalty += (value - self.target_value) ** 2
        
        return penalty * self.weight
    
    def get_satisfaction_score(self, value: float) -> float:
        """Get satisfaction score (0-1, higher is better).
        
        Args:
            value: Property value
            
        Returns:
            Satisfaction score between 0 and 1
        """
        if self.is_satisfied(value):
            # If target specified, score based on distance to target
            if self.target_value is not None:
                # Use Gaussian-like scoring around target
                distance = abs(value - self.target_value)
                # Scale by reasonable range (use constraint range if available)
                scale = 1.0
                if self.min_value is not None and self.max_value is not None:
                    scale = (self.max_value - self.min_value) / 4  # 4-sigma range
                return max(0.0, 1.0 - (distance / scale) ** 2)
            else:
                return 1.0  # Satisfied range constraint
        else:
            # Penalty-based scoring for violations
            penalty = self.compute_penalty(value)
            return max(0.0, 1.0 - penalty)


@dataclass
class ConstraintSet:
    """Collection of property constraints for optimization."""
    
    constraints: List[PropertyConstraint]
    
    def __post_init__(self):
        """Validate constraint set."""
        if not self.constraints:
            raise ValueError("ConstraintSet must contain at least one constraint")
        
        # Check for duplicate constraint names
        names = [c.name for c in self.constraints]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate constraint names not allowed")
    
    def is_satisfied(self, predictions: Dict[str, float]) -> bool:
        """Check if all constraints are satisfied.
        
        Args:
            predictions: Dictionary of property name -> value
            
        Returns:
            True if all constraints are satisfied
        """
        for constraint in self.constraints:
            if constraint.name not in predictions:
                logger.warning(f"Missing prediction for constraint '{constraint.name}'")
                return False
            
            if not constraint.is_satisfied(predictions[constraint.name]):
                return False
        
        return True
    
    def compute_penalty(self, predictions: Dict[str, float]) -> float:
        """Compute total weighted penalty for all constraint violations.
        
        Args:
            predictions: Dictionary of property name -> value
            
        Returns:
            Total penalty (0 if all satisfied, positive if violations)
        """
        total_penalty = 0.0
        
        for constraint in self.constraints:
            if constraint.name in predictions:
                penalty = constraint.compute_penalty(predictions[constraint.name])
                total_penalty += penalty
            else:
                # Missing prediction is heavily penalized
                total_penalty += 100.0 * constraint.weight
                logger.warning(f"Missing prediction for constraint '{constraint.name}'")
        
        return total_penalty
    
    def compute_score(self, predictions: Dict[str, float]) -> float:
        """Compute overall satisfaction score (0-1, higher is better).
        
        Args:
            predictions: Dictionary of property name -> value
            
        Returns:
            Weighted average satisfaction score
        """
        if not self.constraints:
            return 1.0
        
        total_score = 0.0
        total_weight = 0.0
        
        for constraint in self.constraints:
            if constraint.name in predictions:
                score = constraint.get_satisfaction_score(predictions[constraint.name])
                total_score += score * constraint.weight
                total_weight += constraint.weight
            else:
                # Missing prediction gets score of 0
                total_weight += constraint.weight
                logger.warning(f"Missing prediction for constraint '{constraint.name}'")
        
        return total_score / total_weight if total_weight > 0 else 0.0
    
    def get_violations(self, predictions: Dict[str, float]) -> List[str]:
        """Get list of violated constraint names.
        
        Args:
            predictions: Dictionary of property name -> value
            
        Returns:
            List of constraint names that are violated
        """
        violations = []
        
        for constraint in self.constraints:
            if constraint.name not in predictions:
                violations.append(f"{constraint.name} (missing)")
            elif not constraint.is_satisfied(predictions[constraint.name]):
                violations.append(constraint.name)
        
        return violations
    
    def get_constraint_by_name(self, name: str) -> Optional[PropertyConstraint]:
        """Get constraint by name.
        
        Args:
            name: Constraint name
            
        Returns:
            PropertyConstraint or None if not found
        """
        for constraint in self.constraints:
            if constraint.name == name:
                return constraint
        return None


# Common constraint presets for drug-like molecules
def create_lipinski_constraints() -> ConstraintSet:
    """Create Lipinski Rule of Five constraints."""
    return ConstraintSet([
        PropertyConstraint("mw", min_value=150, max_value=500, weight=1.0),
        PropertyConstraint("logp", min_value=-0.4, max_value=5.6, weight=1.0),
        PropertyConstraint("hbd", min_value=0, max_value=5, weight=1.0),
        PropertyConstraint("hba", min_value=0, max_value=10, weight=1.0),
    ])


def create_lead_like_constraints() -> ConstraintSet:
    """Create lead-like compound constraints (more restrictive than Lipinski)."""
    return ConstraintSet([
        PropertyConstraint("mw", min_value=100, max_value=350, weight=1.0),
        PropertyConstraint("logp", min_value=-1.0, max_value=3.5, weight=1.0),
        PropertyConstraint("tpsa", min_value=20, max_value=90, weight=1.0),
        PropertyConstraint("rotatable_bonds", min_value=0, max_value=7, weight=1.0),
    ])


def create_cns_constraints() -> ConstraintSet:
    """Create CNS drug-like constraints."""
    return ConstraintSet([
        PropertyConstraint("mw", min_value=150, max_value=450, weight=1.0),
        PropertyConstraint("logp", min_value=1.0, max_value=4.0, weight=1.0),
        PropertyConstraint("tpsa", min_value=20, max_value=90, weight=1.0),
        PropertyConstraint("hbd", min_value=0, max_value=3, weight=1.0),
        PropertyConstraint("hba", min_value=0, max_value=7, weight=1.0),
    ])


def create_admet_constraints() -> ConstraintSet:
    """Create ADMET-focused constraints."""
    return ConstraintSet([
        PropertyConstraint("herg", min_value=0.0, max_value=0.3, weight=2.0),  # Low hERG risk
        PropertyConstraint("logs", min_value=-4.0, max_value=-2.0, weight=1.5),  # Good solubility
        PropertyConstraint("logp", min_value=0.0, max_value=4.0, weight=1.0),  # Balanced lipophilicity
    ])
