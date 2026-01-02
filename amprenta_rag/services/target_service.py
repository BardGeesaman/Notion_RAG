"""Target management service layer."""
from uuid import UUID
from typing import Optional, List
from sqlalchemy.orm import Session
from sqlalchemy import and_

from amprenta_rag.database.models import Target, compound_target, Compound, Experiment
from amprenta_rag.database.session import db_session
import logging

logger = logging.getLogger(__name__)


def create_target(
    name: str,
    gene_symbol: Optional[str] = None,
    uniprot_id: Optional[str] = None,
    target_class: Optional[str] = None,
    target_family: Optional[str] = None,
    description: Optional[str] = None,
    db: Optional[Session] = None
) -> Target:
    """Create a new target."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        target = Target(
            name=name,
            gene_symbol=gene_symbol,
            uniprot_id=uniprot_id,
            target_class=target_class,
            target_family=target_family,
            description=description,
            lifecycle_status="active"
        )
        
        db.add(target)
        db.commit()
        db.refresh(target)
        
        logger.info(f"Created target: {target.name} ({target.id})")
        return target
        
    except Exception as e:
        db.rollback()
        logger.error(f"Failed to create target {name}: {e}")
        raise
    finally:
        if not use_session:
            db.close()


def get_target(target_id: UUID, db: Optional[Session] = None) -> Optional[Target]:
    """Get target by ID."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        target = db.query(Target).filter(Target.id == target_id).first()
        if target:
            db.expunge(target)
        return target
    finally:
        if not use_session:
            db.close()


def list_targets(
    target_class: Optional[str] = None,
    target_family: Optional[str] = None,
    lifecycle_status: str = "active",
    skip: int = 0,
    limit: int = 100,
    db: Optional[Session] = None
) -> List[Target]:
    """List targets with optional filters."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        query = db.query(Target).filter(Target.lifecycle_status == lifecycle_status)
        
        if target_class:
            query = query.filter(Target.target_class == target_class)
        if target_family:
            query = query.filter(Target.target_family == target_family)
            
        targets = query.offset(skip).limit(limit).all()
        
        # Detach from session
        for target in targets:
            db.expunge(target)
            
        return targets
    finally:
        if not use_session:
            db.close()


def update_target(target_id: UUID, updates: dict, db: Optional[Session] = None) -> Optional[Target]:
    """Update target fields."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        target = db.query(Target).filter(Target.id == target_id).first()
        if not target:
            return None
            
        # Update allowed fields
        allowed_fields = {
            'name', 'description', 'gene_symbol', 'uniprot_id', 
            'target_class', 'target_family', 'validation_status',
            'lifecycle_status', 'druggability_score', 'druggability_source',
            'pocket_count', 'validation_evidence', 'is_synthetic',
            'chembl_id', 'ensembl_id'
        }
        
        for field, value in updates.items():
            if field in allowed_fields and hasattr(target, field):
                setattr(target, field, value)
        
        db.commit()
        db.refresh(target)
        db.expunge(target)
        
        logger.info(f"Updated target {target_id}: {list(updates.keys())}")
        return target
        
    except Exception as e:
        db.rollback()
        logger.error(f"Failed to update target {target_id}: {e}")
        raise
    finally:
        if not use_session:
            db.close()


def get_target_assays(target_id: UUID, db: Optional[Session] = None) -> List[dict]:
    """Get experiments/assays linked to target."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        # Find experiments that have compounds linked to this target
        assays = db.query(Experiment).join(
            compound_target, 
            compound_target.c.assay_id == Experiment.id
        ).filter(
            compound_target.c.target_id == target_id
        ).distinct().all()
        
        result = []
        for assay in assays:
            result.append({
                "id": assay.id,
                "name": assay.name,
                "description": assay.description,
                "experiment_type": assay.experiment_type,
                "created_at": assay.created_at
            })
            
        return result
    finally:
        if not use_session:
            db.close()


def get_target_compounds(target_id: UUID, db: Optional[Session] = None) -> List[dict]:
    """Get compounds with activity data for target."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        # Query compound-target associations
        results = db.query(
            Compound.id,
            Compound.smiles,
            Compound.compound_id,
            compound_target.c.activity_type,
            compound_target.c.activity_value,
            compound_target.c.activity_units,
            compound_target.c.assay_id
        ).join(
            compound_target,
            Compound.id == compound_target.c.compound_id
        ).filter(
            compound_target.c.target_id == target_id
        ).all()
        
        compounds = []
        for row in results:
            compounds.append({
                "compound_id": row.id,
                "compound_name": row.compound_id,
                "smiles": row.smiles,
                "activity_type": row.activity_type,
                "activity_value": row.activity_value,
                "activity_units": row.activity_units,
                "assay_id": row.assay_id
            })
            
        return compounds
    finally:
        if not use_session:
            db.close()


def calculate_druggability(target_id: UUID, db: Optional[Session] = None) -> dict:
    """Calculate druggability score for target."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        target = db.query(Target).filter(Target.id == target_id).first()
        if not target:
            return {"error": "Target not found"}
            
        score = 0.0
        factors = {}
        
        # Has UniProt ID: +0.2
        if target.uniprot_id:
            score += 0.2
            factors["uniprot_id"] = 0.2
        
        # Has pocket_count > 0: +0.3
        if target.pocket_count and target.pocket_count > 0:
            score += 0.3
            factors["pockets"] = 0.3
        
        # validation_status == 'validated' or 'clinical': +0.3
        if target.validation_status in ['validated', 'clinical']:
            score += 0.3
            factors["validation"] = 0.3
        
        # Has compounds with activity: +0.2
        compound_count = db.query(compound_target).filter(
            compound_target.c.target_id == target_id
        ).count()
        
        if compound_count > 0:
            score += 0.2
            factors["compounds"] = 0.2
            factors["compound_count"] = compound_count
        
        # Update target with calculated score
        target.druggability_score = min(score, 1.0)  # Cap at 1.0
        target.druggability_source = "calculated"
        db.commit()
        
        return {
            "target_id": target_id,
            "score": target.druggability_score,
            "factors": factors
        }
        
    except Exception as e:
        logger.error(f"Failed to calculate druggability for target {target_id}: {e}")
        return {"error": str(e)}
    finally:
        if not use_session:
            db.close()


def get_target_by_name(name: str, db: Optional[Session] = None) -> Optional[Target]:
    """Get target by name."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        target = db.query(Target).filter(Target.name == name).first()
        if target:
            db.expunge(target)
        return target
    finally:
        if not use_session:
            db.close()


def search_targets(
    query: str,
    limit: int = 20,
    db: Optional[Session] = None
) -> List[Target]:
    """Search targets by name, gene symbol, or description."""
    use_session = db is not None
    
    if not use_session:
        db = next(db_session())
    
    try:
        search_term = f"%{query}%"
        targets = db.query(Target).filter(
            and_(
                Target.lifecycle_status == "active",
                db.or_(
                    Target.name.ilike(search_term),
                    Target.gene_symbol.ilike(search_term),
                    Target.description.ilike(search_term)
                )
            )
        ).limit(limit).all()
        
        # Detach from session
        for target in targets:
            db.expunge(target)
            
        return targets
    finally:
        if not use_session:
            db.close()
