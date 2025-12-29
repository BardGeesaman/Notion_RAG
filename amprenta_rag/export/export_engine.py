"""
Data export engine for datasets, experiments, and compounds.

Provides export to CSV, Excel, JSON with manifest generation and packaging.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from datetime import datetime
from typing import Any, Dict, List

import pandas as pd
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Compound, Dataset, Experiment
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def export_dataset(
    db: Session,
    dataset_id: str,
    format: str = "csv",
    include_metadata: bool = True,
) -> bytes:
    """
    Export dataset to specified format.

    Args:
        db: Database session
        dataset_id: Dataset UUID
        format: Output format (csv, excel, json)
        include_metadata: Whether to include metadata header

    Returns:
        Bytes of exported data
    """
    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    
    if not dataset:
        raise ValueError(f"Dataset {dataset_id} not found")
    
    # Build data dict
    data = {
        "dataset_id": str(dataset.id),
        "name": dataset.name,
        "omics_type": dataset.omics_type,
        "feature_count": len(dataset.features) if dataset.features else 0,
    }
    
    # Add features if available
    features_data = []
    for feat in dataset.features or []:
        features_data.append({
            "feature_name": feat.name,
            "value": feat.external_ids.get("value", 1.0) if feat.external_ids else 1.0,
        })
    
    if format == "csv":
        df = pd.DataFrame(features_data) if features_data else pd.DataFrame()
        output = io.BytesIO()
        
        if include_metadata:
            output.write(f"# Dataset: {dataset.name}\n".encode())
            output.write(f"# ID: {dataset_id}\n".encode())
        
        df.to_csv(output, index=False)
        return output.getvalue()
    
    elif format == "excel":
        output = io.BytesIO()
        df = pd.DataFrame(features_data) if features_data else pd.DataFrame()
        df.to_excel(output, index=False, sheet_name="Features")
        return output.getvalue()
    
    elif format == "json":
        export_data = {**data, "features": features_data}
        return json.dumps(export_data, indent=2).encode()
    
    else:
        raise ValueError(f"Unsupported format: {format}")


def export_experiment(
    db: Session,
    experiment_id: str,
    format: str = "csv",
    include_datasets: bool = True,
) -> bytes:
    """
    Export experiment data.

    Args:
        db: Database session
        experiment_id: Experiment UUID
        format: Output format (csv, excel, json)
        include_datasets: Whether to include related datasets

    Returns:
        Bytes of exported data
    """
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()
    
    if not experiment:
        raise ValueError(f"Experiment {experiment_id} not found")
    
    data = {
        "experiment_id": str(experiment.id),
        "name": experiment.name,
        "type": experiment.type,
    }
    
    if format == "json":
        if include_datasets:
            data["datasets"] = [str(ds.id) for ds in experiment.datasets or []]
        return json.dumps(data, indent=2).encode()
    
    elif format == "csv":
        df = pd.DataFrame([data])
        return df.to_csv(index=False).encode()
    
    else:
        raise ValueError(f"Unsupported format: {format}")


def export_compounds(
    db: Session,
    compound_ids: List[str],
    format: str = "csv",
) -> bytes:
    """
    Export compounds list.

    Args:
        db: Database session
        compound_ids: List of compound UUIDs
        format: Output format (csv, excel, json)

    Returns:
        Bytes of exported data
    """
    compounds = db.query(Compound).filter(Compound.id.in_(compound_ids)).all()
    
    data = []
    for comp in compounds:
        data.append({
            "compound_id": comp.compound_id or str(comp.id),
            "smiles": comp.canonical_smiles or comp.smiles,
            "name": comp.name,
        })
    
    if format == "csv":
        df = pd.DataFrame(data)
        return df.to_csv(index=False).encode()
    
    elif format == "json":
        return json.dumps(data, indent=2).encode()
    
    else:
        raise ValueError(f"Unsupported format: {format}")


def generate_manifest(items: List[str], checksums: Dict[str, str]) -> dict:
    """
    Generate export manifest with checksums.

    Args:
        items: List of exported file names
        checksums: Dictionary mapping file names to SHA256 hashes

    Returns:
        Manifest dictionary
    """
    return {
        "export_timestamp": datetime.utcnow().isoformat(),
        "version": "1.0",
        "items": items,
        "checksums": checksums,
    }


def create_export_package(exports: Dict[str, bytes], format: str = "zip") -> bytes:
    """
    Bundle multiple exports into a package.

    Args:
        exports: Dictionary mapping file names to content bytes
        format: Package format (currently only zip supported)

    Returns:
        Bytes of packaged exports
    """
    if format != "zip":
        raise ValueError(f"Unsupported package format: {format}")
    
    output = io.BytesIO()
    
    # Calculate checksums
    checksums = {}
    for filename, content in exports.items():
        sha256 = hashlib.sha256(content).hexdigest()
        checksums[filename] = sha256
    
    # Create manifest
    manifest = generate_manifest(list(exports.keys()), checksums)
    
    # Create ZIP
    with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as zf:
        # Add manifest
        zf.writestr("manifest.json", json.dumps(manifest, indent=2))
        
        # Add exports
        for filename, content in exports.items():
            zf.writestr(filename, content)
    
    return output.getvalue()

