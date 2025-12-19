"""Data export utilities for various formats."""
from __future__ import annotations

import io
import json
from typing import List

import pandas as pd

from amprenta_rag.database.models import Experiment, Compound, Signature
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def export_to_csv(df: pd.DataFrame) -> bytes:
    """Export DataFrame to CSV format."""
    buffer = io.BytesIO()
    df.to_csv(buffer, index=False)
    return buffer.getvalue()


def export_to_json(df: pd.DataFrame) -> bytes:
    """Export DataFrame to JSON format."""
    records = df.to_dict(orient="records")
    json_str = json.dumps(records, indent=2, default=str)
    return json_str.encode("utf-8")


def export_to_excel(df: pd.DataFrame) -> bytes:
    """Export DataFrame to Excel format."""
    try:
        import openpyxl  # noqa: F401
    except ImportError:
        raise ImportError("openpyxl is required for Excel export. Install with: pip install openpyxl")

    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Data")
    return buffer.getvalue()


def export_experiments(experiment_ids: List[str], format: str, db) -> bytes:
    """
    Export experiments to specified format.

    Args:
        experiment_ids: List of experiment UUIDs
        format: Export format ("csv", "json", "excel")
        db: Database session

    Returns:
        Exported data as bytes
    """
    experiments = db.query(Experiment).filter(Experiment.id.in_(experiment_ids)).all()

    if not experiments:
        # Return empty DataFrame
        df = pd.DataFrame()
    else:
        data = []
        for exp in experiments:
            data.append({
                "id": str(exp.id),
                "name": exp.name,
                "type": exp.type,
                "description": exp.description,
                "disease": exp.disease,
                "matrix": exp.matrix,
                "model_systems": exp.model_systems,
                "targets": exp.targets,
                "modality": exp.modality,
                "stage": exp.stage,
                "biomarker_role": exp.biomarker_role,
                "treatment_arms": exp.treatment_arms,
                "design_type": exp.design_type,
                "design_metadata": json.dumps(exp.design_metadata) if exp.design_metadata else None,
                "sample_groups": json.dumps(exp.sample_groups) if exp.sample_groups else None,
                "created_at": exp.created_at.isoformat() if exp.created_at else None,
                "updated_at": exp.updated_at.isoformat() if exp.updated_at else None,
            })
        df = pd.DataFrame(data)

    if format == "csv":
        return export_to_csv(df)
    elif format == "json":
        return export_to_json(df)
    elif format == "excel":
        return export_to_excel(df)
    else:
        raise ValueError(f"Unsupported format: {format}")


def export_compounds(compound_ids: List[str], format: str, db) -> bytes:
    """
    Export compounds to specified format.

    Args:
        compound_ids: List of compound IDs (compound_id field, not UUID)
        format: Export format ("csv", "json", "excel")
        db: Database session

    Returns:
        Exported data as bytes
    """
    compounds = db.query(Compound).filter(Compound.compound_id.in_(compound_ids)).all()

    if not compounds:
        df = pd.DataFrame()
    else:
        data = []
        for comp in compounds:
            data.append({
                "compound_id": comp.compound_id,
                "smiles": comp.smiles,
                "canonical_smiles": comp.canonical_smiles,
                "inchi_key": comp.inchi_key,
                "molecular_formula": comp.molecular_formula,
                "molecular_weight": comp.molecular_weight,
                "logp": comp.logp,
                "hbd_count": comp.hbd_count,
                "hba_count": comp.hba_count,
                "rotatable_bonds": comp.rotatable_bonds,
                "external_ids": json.dumps(comp.external_ids) if comp.external_ids else None,
                "created_at": comp.created_at.isoformat() if comp.created_at else None,
                "updated_at": comp.updated_at.isoformat() if comp.updated_at else None,
            })
        df = pd.DataFrame(data)

    if format == "csv":
        return export_to_csv(df)
    elif format == "json":
        return export_to_json(df)
    elif format == "excel":
        return export_to_excel(df)
    else:
        raise ValueError(f"Unsupported format: {format}")


def export_signatures(signature_ids: List[str], format: str, db) -> bytes:
    """
    Export signatures to specified format.

    Args:
        signature_ids: List of signature UUIDs
        format: Export format ("csv", "json", "excel")
        db: Database session

    Returns:
        Exported data as bytes
    """
    signatures = db.query(Signature).filter(Signature.id.in_(signature_ids)).all()

    if not signatures:
        df = pd.DataFrame()
    else:
        data = []
        for sig in signatures:
            data.append({
                "id": str(sig.id),
                "name": sig.name,
                "short_id": sig.short_id,
                "description": sig.description,
                "modalities": sig.modalities,
                "biomarker_role": sig.biomarker_role,
                "phenotype_axes": sig.phenotype_axes,
                "data_ownership": sig.data_ownership,
                "created_at": sig.created_at.isoformat() if sig.created_at else None,
                "updated_at": sig.updated_at.isoformat() if sig.updated_at else None,
            })
        df = pd.DataFrame(data)

    if format == "csv":
        return export_to_csv(df)
    elif format == "json":
        return export_to_json(df)
    elif format == "excel":
        return export_to_excel(df)
    else:
        raise ValueError(f"Unsupported format: {format}")
