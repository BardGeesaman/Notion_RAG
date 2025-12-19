"""Batch operations utilities for bulk actions on entities."""
from __future__ import annotations

from typing import List

from amprenta_rag.database.models import Experiment, Compound, Signature, Dataset
from amprenta_rag.utils.data_export import export_experiments, export_compounds, export_signatures
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def batch_export(entity_type: str, ids: List[str], db, format: str = "csv") -> bytes:
    """
    Export multiple entities of the same type.

    Args:
        entity_type: Type of entity ("experiment", "compound", "signature", "dataset")
        ids: List of entity IDs (UUIDs as strings)
        db: Database session
        format: Export format ("csv", "json", "excel")

    Returns:
        Exported data as bytes
    """
    if entity_type == "experiment":
        return export_experiments(ids, format, db)
    elif entity_type == "compound":
        return export_compounds(ids, format, db)
    elif entity_type == "signature":
        return export_signatures(ids, format, db)
    elif entity_type == "dataset":
        # Use data_export functions for datasets
        from amprenta_rag.utils.data_export import export_to_csv, export_to_json, export_to_excel
        import pandas as pd

        datasets = db.query(Dataset).filter(Dataset.id.in_(ids)).all()
        if not datasets:
            return b""

        data = []
        for ds in datasets:
            data.append({
                "id": str(ds.id),
                "name": ds.name,
                "omics_type": ds.omics_type,
                "description": ds.description,
                "organism": ", ".join(ds.organism) if ds.organism else None,
                "disease": ", ".join(ds.disease) if ds.disease else None,
                "created_at": ds.created_at.isoformat() if ds.created_at else None,
            })
        df = pd.DataFrame(data)

        if format == "csv":
            return export_to_csv(df)
        elif format == "json":
            return export_to_json(df)
        else:
            return export_to_excel(df)
    else:
        raise ValueError(f"Unsupported entity type: {entity_type}")


def batch_delete(entity_type: str, ids: List[str], db) -> int:
    """
    Delete multiple entities of the same type.

    Args:
        entity_type: Type of entity ("experiment", "compound", "signature", "dataset")
        ids: List of entity IDs (UUIDs as strings)
        db: Database session

    Returns:
        Number of entities deleted
    """
    from uuid import UUID

    if not ids:
        return 0

    uuid_ids = [UUID(id_str) if isinstance(id_str, str) else id_str for id_str in ids]

    if entity_type == "experiment":
        deleted = db.query(Experiment).filter(Experiment.id.in_(uuid_ids)).delete(synchronize_session=False)
    elif entity_type == "compound":
        deleted = db.query(Compound).filter(Compound.id.in_(uuid_ids)).delete(synchronize_session=False)
    elif entity_type == "signature":
        deleted = db.query(Signature).filter(Signature.id.in_(uuid_ids)).delete(synchronize_session=False)
    elif entity_type == "dataset":
        deleted = db.query(Dataset).filter(Dataset.id.in_(uuid_ids)).delete(synchronize_session=False)
    else:
        raise ValueError(f"Unsupported entity type: {entity_type}")

    db.commit()
    logger.info("[BATCH_OPS] Deleted %d %s entities", deleted, entity_type)

    return deleted
