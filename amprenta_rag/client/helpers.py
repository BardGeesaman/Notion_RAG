"""Convenience helpers for notebook write-back operations."""

import os
from urllib.parse import parse_qs, urlparse
from typing import Optional

def get_context_from_url() -> dict:
    """Extract context params passed from Streamlit.
    
    Returns dict with experiment_id, dataset_id, compound_id if present.
    """
    # JupyterHub passes original URL via JUPYTERHUB_OAUTH_CALLBACK_URL or similar
    # For now, check environment variables that might be set
    return {
        "experiment_id": os.environ.get("EXPERIMENT_ID"),
        "dataset_id": os.environ.get("DATASET_ID"),
        "compound_id": os.environ.get("COMPOUND_ID"),
    }

def save_annotation(
    client,
    text: str,
    annotation_type: str = "notebook",
    experiment_id: str = None,
    dataset_id: str = None,
    signature_id: str = None,
    compound_id: str = None,
) -> dict:
    """Save an annotation to an entity.
    
    Args:
        client: RAGClient instance
        text: Annotation text
        annotation_type: Type of annotation (default: "notebook")
        experiment_id: Target experiment UUID
        dataset_id: Target dataset UUID
        signature_id: Target signature UUID
        compound_id: Target compound UUID
        
    Returns:
        Created annotation dict
    """
    if experiment_id:
        return client.experiments.annotate(experiment_id, text, annotation_type)
    elif dataset_id:
        return client.datasets.annotate(dataset_id, text, annotation_type)
    elif signature_id:
        return client.signatures.annotate(signature_id, text, annotation_type)
    elif compound_id:
        return client.compounds.annotate(compound_id, text, annotation_type)
    else:
        raise ValueError("Must provide at least one entity ID")

