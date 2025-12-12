"""Data import utilities for bulk importing entities."""
from __future__ import annotations

import json
from typing import List, Optional
from uuid import UUID

import pandas as pd

from amprenta_rag.database.models import Experiment, Compound, Sample
from amprenta_rag.chemistry.normalization import normalize_smiles, compute_molecular_descriptors
from amprenta_rag.chemistry.registration import check_duplicate
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


# Required columns for each entity type
REQUIRED_COLUMNS = {
    "experiment": ["name"],
    "compound": ["smiles"],
    "sample": ["name"],
}


def validate_import_data(df: pd.DataFrame, entity_type: str) -> List[str]:
    """
    Validate import data has required columns.
    
    Args:
        df: DataFrame to validate
        entity_type: Type of entity ("experiment", "compound", "sample")
        
    Returns:
        List of error messages (empty if valid)
    """
    errors = []
    
    if df.empty:
        errors.append("DataFrame is empty")
        return errors
    
    required = REQUIRED_COLUMNS.get(entity_type)
    if not required:
        errors.append(f"Unknown entity type: {entity_type}")
        return errors
    
    missing = [col for col in required if col not in df.columns]
    if missing:
        errors.append(f"Missing required columns: {', '.join(missing)}")
    
    # Check for empty required values
    for col in required:
        if col in df.columns:
            empty_count = df[col].isna().sum() + (df[col] == "").sum()
            if empty_count > 0:
                errors.append(f"Column '{col}' has {empty_count} empty values")
    
    return errors


def import_experiments(df: pd.DataFrame, db, user_id: Optional[str] = None) -> dict:
    """
    Import experiments from DataFrame.
    
    Args:
        df: DataFrame with experiment data
        db: Database session
        user_id: User ID for created_by_id
        
    Returns:
        Dict with "created" count and "errors" list
    """
    errors = validate_import_data(df, "experiment")
    if errors:
        return {"created": 0, "errors": errors}
    
    created = 0
    import_errors = []
    
    for idx, row in df.iterrows():
        try:
            # Parse array columns
            def parse_array(value):
                if pd.isna(value) or value == "":
                    return None
                if isinstance(value, str):
                    return [x.strip() for x in value.split(",") if x.strip()]
                return value
            
            experiment = Experiment(
                name=str(row["name"]),
                type=row.get("type"),
                description=row.get("description"),
                disease=parse_array(row.get("disease")),
                matrix=parse_array(row.get("matrix")),
                model_systems=parse_array(row.get("model_systems")),
                targets=parse_array(row.get("targets")),
                modality=parse_array(row.get("modality")),
                stage=row.get("stage"),
                biomarker_role=parse_array(row.get("biomarker_role")),
                treatment_arms=parse_array(row.get("treatment_arms")),
                design_type=row.get("design_type"),
                design_metadata=json.loads(row["design_metadata"]) if pd.notna(row.get("design_metadata")) else None,
                sample_groups=json.loads(row["sample_groups"]) if pd.notna(row.get("sample_groups")) else None,
                created_by_id=UUID(user_id) if user_id and user_id != "test" else None,
            )
            
            db.add(experiment)
            created += 1
            
        except Exception as e:
            import_errors.append(f"Row {idx + 1}: {str(e)}")
            logger.error("[IMPORT] Error importing experiment row %d: %r", idx + 1, e)
    
    try:
        db.commit()
        logger.info("[IMPORT] Imported %d experiments", created)
    except Exception as e:
        db.rollback()
        import_errors.append(f"Commit failed: {str(e)}")
        logger.error("[IMPORT] Commit failed: %r", e)
    
    return {"created": created, "errors": import_errors}


def import_compounds(df: pd.DataFrame, db) -> dict:
    """
    Import compounds from DataFrame.
    
    Args:
        df: DataFrame with compound data (must have 'smiles' column)
        db: Database session
        
    Returns:
        Dict with "created" count and "duplicates" count
    """
    errors = validate_import_data(df, "compound")
    if errors:
        return {"created": 0, "duplicates": 0, "errors": errors}
    
    created = 0
    duplicates = 0
    
    for idx, row in df.iterrows():
        try:
            smiles = str(row["smiles"]).strip()
            if not smiles:
                continue
            
            # Check for duplicate
            existing = check_duplicate(smiles)
            if existing:
                duplicates += 1
                logger.debug("[IMPORT] Duplicate compound: %s", existing)
                continue
            
            # Normalize SMILES and compute descriptors
            canonical, inchi_key, formula = normalize_smiles(smiles)
            descriptors = compute_molecular_descriptors(canonical or smiles)
            
            # Get compound_id from row or generate
            compound_id = row.get("compound_id")
            if not compound_id or pd.isna(compound_id):
                # Generate next corporate ID
                result = db.query(Compound).filter(
                    Compound.compound_id.like("AMP-%")
                ).order_by(Compound.compound_id.desc()).first()
                if result:
                    try:
                        current_num = int(result.compound_id.split("-")[1])
                        next_num = current_num + 1
                    except (IndexError, ValueError):
                        next_num = 1
                else:
                    next_num = 1
                compound_id = f"AMP-{next_num:05d}"
            
            compound = Compound(
                compound_id=str(compound_id),
                smiles=canonical or smiles,
                inchi_key=inchi_key,
                canonical_smiles=canonical,
                molecular_formula=formula,
                molecular_weight=descriptors.get("molecular_weight"),
                logp=descriptors.get("logp"),
                hbd_count=descriptors.get("hbd_count"),
                hba_count=descriptors.get("hba_count"),
                rotatable_bonds=descriptors.get("rotatable_bonds"),
            )
            
            db.add(compound)
            created += 1
            
        except Exception as e:
            logger.error("[IMPORT] Error importing compound row %d: %r", idx + 1, e)
    
    try:
        db.commit()
        logger.info("[IMPORT] Imported %d compounds (%d duplicates skipped)", created, duplicates)
    except Exception as e:
        db.rollback()
        logger.error("[IMPORT] Commit failed: %r", e)
        return {"created": 0, "duplicates": duplicates, "errors": [f"Commit failed: {str(e)}"]}
    
    return {"created": created, "duplicates": duplicates}


def import_samples(df: pd.DataFrame, db, user_id: Optional[str] = None) -> dict:
    """
    Import samples from DataFrame.
    
    Args:
        df: DataFrame with sample data
        db: Database session
        user_id: User ID for created_by_id
        
    Returns:
        Dict with "created" count and "errors" list
    """
    errors = validate_import_data(df, "sample")
    if errors:
        return {"created": 0, "errors": errors}
    
    created = 0
    import_errors = []
    
    for idx, row in df.iterrows():
        try:
            # Resolve storage_location_id if location name provided
            storage_location_id = None
            if pd.notna(row.get("storage_location")):
                from amprenta_rag.database.models import StorageLocation
                location = db.query(StorageLocation).filter(
                    StorageLocation.name == str(row["storage_location"])
                ).first()
                if location:
                    storage_location_id = location.id
                else:
                    import_errors.append(f"Row {idx + 1}: Storage location '{row['storage_location']}' not found")
            
            # Resolve experiment_id if experiment name provided
            experiment_id = None
            if pd.notna(row.get("experiment")):
                experiment = db.query(Experiment).filter(
                    Experiment.name == str(row["experiment"])
                ).first()
                if experiment:
                    experiment_id = experiment.id
                else:
                    import_errors.append(f"Row {idx + 1}: Experiment '{row['experiment']}' not found")
            
            # Resolve parent_sample_id if parent sample name provided
            parent_sample_id = None
            if pd.notna(row.get("parent_sample")):
                parent = db.query(Sample).filter(
                    Sample.name == str(row["parent_sample"])
                ).first()
                if parent:
                    parent_sample_id = parent.id
                else:
                    import_errors.append(f"Row {idx + 1}: Parent sample '{row['parent_sample']}' not found")
            
            sample = Sample(
                name=str(row["name"]),
                sample_type=row.get("sample_type"),
                barcode=row.get("barcode"),
                storage_location_id=storage_location_id,
                position=row.get("position"),
                parent_sample_id=parent_sample_id,
                experiment_id=experiment_id,
                quantity=float(row["quantity"]) if pd.notna(row.get("quantity")) else None,
                unit=row.get("unit"),
                status=row.get("status", "available"),
                created_by_id=UUID(user_id) if user_id and user_id != "test" else None,
                notes=row.get("notes"),
            )
            
            db.add(sample)
            created += 1
            
        except Exception as e:
            import_errors.append(f"Row {idx + 1}: {str(e)}")
            logger.error("[IMPORT] Error importing sample row %d: %r", idx + 1, e)
    
    try:
        db.commit()
        logger.info("[IMPORT] Imported %d samples", created)
    except Exception as e:
        db.rollback()
        import_errors.append(f"Commit failed: {str(e)}")
        logger.error("[IMPORT] Commit failed: %r", e)
    
    return {"created": created, "errors": import_errors}


def get_template(entity_type: str) -> pd.DataFrame:
    """
    Get empty DataFrame template with required columns for entity type.
    
    Args:
        entity_type: Type of entity ("experiment", "compound", "sample")
        
    Returns:
        Empty DataFrame with appropriate columns
    """
    templates = {
        "experiment": pd.DataFrame(columns=[
            "name", "type", "description", "disease", "matrix", "model_systems",
            "targets", "modality", "stage", "biomarker_role", "treatment_arms",
            "design_type", "design_metadata", "sample_groups"
        ]),
        "compound": pd.DataFrame(columns=[
            "compound_id", "smiles", "inchi_key", "canonical_smiles",
            "molecular_formula", "molecular_weight", "logp", "hbd_count",
            "hba_count", "rotatable_bonds"
        ]),
        "sample": pd.DataFrame(columns=[
            "name", "sample_type", "barcode", "storage_location", "position",
            "parent_sample", "experiment", "quantity", "unit", "status", "notes"
        ]),
    }
    
    return templates.get(entity_type, pd.DataFrame())
