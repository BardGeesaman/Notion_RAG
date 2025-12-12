"""Data validation utilities."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List
from uuid import UUID

from amprenta_rag.database.models import Experiment, Compound
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class ValidationIssue:
    """Represents a validation issue."""
    entity_type: str
    entity_id: str
    field: str
    issue: str
    severity: str  # "error" or "warning"


def validate_experiment(exp_id: str, db) -> List[ValidationIssue]:
    """
    Validate an experiment.
    
    Args:
        exp_id: UUID of the experiment
        db: Database session
        
    Returns:
        List of ValidationIssue objects
    """
    issues = []
    
    experiment = db.query(Experiment).filter(
        Experiment.id == (UUID(exp_id) if isinstance(exp_id, str) else exp_id)
    ).first()
    
    if not experiment:
        issues.append(ValidationIssue(
            entity_type="experiment",
            entity_id=exp_id,
            field="id",
            issue="Experiment not found",
            severity="error"
        ))
        return issues
    
    # Check name not empty
    if not experiment.name or not experiment.name.strip():
        issues.append(ValidationIssue(
            entity_type="experiment",
            entity_id=str(experiment.id),
            field="name",
            issue="Name is empty",
            severity="error"
        ))
    
    # Check design_type is valid
    valid_design_types = [
        "case_control", "time_course", "intervention", "dose_response",
        "multi_factorial", "observational", None
    ]
    if experiment.design_type and experiment.design_type not in valid_design_types:
        issues.append(ValidationIssue(
            entity_type="experiment",
            entity_id=str(experiment.id),
            field="design_type",
            issue=f"Invalid design_type: {experiment.design_type}",
            severity="warning"
        ))
    
    return issues


def validate_compound(comp_id: str, db) -> List[ValidationIssue]:
    """
    Validate a compound.
    
    Args:
        comp_id: UUID of the compound
        db: Database session
        
    Returns:
        List of ValidationIssue objects
    """
    issues = []
    
    compound = db.query(Compound).filter(
        Compound.id == (UUID(comp_id) if isinstance(comp_id, str) else comp_id)
    ).first()
    
    if not compound:
        issues.append(ValidationIssue(
            entity_type="compound",
            entity_id=comp_id,
            field="id",
            issue="Compound not found",
            severity="error"
        ))
        return issues
    
    # Check SMILES is valid using RDKit
    if compound.smiles:
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(compound.smiles)
            if mol is None:
                issues.append(ValidationIssue(
                    entity_type="compound",
                    entity_id=str(compound.id),
                    field="smiles",
                    issue="Invalid SMILES string",
                    severity="error"
                ))
        except ImportError:
            # RDKit not available, skip SMILES validation
            pass
        except Exception as e:
            issues.append(ValidationIssue(
                entity_type="compound",
                entity_id=str(compound.id),
                field="smiles",
                issue=f"Error validating SMILES: {str(e)}",
                severity="warning"
            ))
    else:
        issues.append(ValidationIssue(
            entity_type="compound",
            entity_id=str(compound.id),
            field="smiles",
            issue="SMILES is missing",
            severity="error"
        ))
    
    # Check molecular_weight is computed
    if compound.smiles and not compound.molecular_weight:
        issues.append(ValidationIssue(
            entity_type="compound",
            entity_id=str(compound.id),
            field="molecular_weight",
            issue="Molecular weight not computed",
            severity="warning"
        ))
    
    return issues


def run_all_validations(db) -> List[ValidationIssue]:
    """
    Run validation checks on all experiments and compounds.
    
    Args:
        db: Database session
        
    Returns:
        List of all ValidationIssue objects
    """
    all_issues = []
    
    # Validate all experiments
    experiments = db.query(Experiment).all()
    logger.info("[VALIDATION] Validating %d experiments", len(experiments))
    for exp in experiments:
        issues = validate_experiment(str(exp.id), db)
        all_issues.extend(issues)
    
    # Validate all compounds
    compounds = db.query(Compound).all()
    logger.info("[VALIDATION] Validating %d compounds", len(compounds))
    for comp in compounds:
        issues = validate_compound(str(comp.id), db)
        all_issues.extend(issues)
    
    logger.info("[VALIDATION] Found %d total issues (%d errors, %d warnings)",
                len(all_issues),
                sum(1 for i in all_issues if i.severity == "error"),
                sum(1 for i in all_issues if i.severity == "warning"))
    
    return all_issues
