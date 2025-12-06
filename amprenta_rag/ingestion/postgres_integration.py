"""
Postgres integration for ingestion pipelines.

This module provides functions for ingestion modules to write directly to
Postgres as the system of record, replacing Notion as the primary database.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.config import get_config
from amprenta_rag.database.base import get_db
from amprenta_rag.database.crud import (
    create_dataset,
    get_or_create_dataset,
    get_or_create_feature,
    get_or_create_program,
    get_or_create_experiment,
    link_dataset_to_features,
    link_dataset_to_program,
    link_dataset_to_experiment,
    bulk_get_or_create_features,
)
from amprenta_rag.database.models import Dataset, Feature, Program, Experiment
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)


def create_or_update_dataset_in_postgres(
    name: str,
    omics_type: OmicsType,
    file_paths: Optional[List[str]] = None,
    description: Optional[str] = None,
    organism: Optional[List[str]] = None,
    sample_type: Optional[List[str]] = None,
    disease: Optional[List[str]] = None,
    notion_page_id: Optional[str] = None,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
    db: Optional[Session] = None,
) -> Dataset:
    """
    Create or update a dataset in Postgres.
    
    Args:
        name: Dataset name
        omics_type: OmicsType enum value
        file_paths: List of file paths
        description: Dataset description
        organism: List of organism names
        sample_type: List of sample types
        disease: List of disease terms
        notion_page_id: Notion page ID for migration compatibility
        program_ids: List of program IDs (Notion or Postgres UUIDs)
        experiment_ids: List of experiment IDs (Notion or Postgres UUIDs)
        db: Database session (optional, will create if not provided)
    
    Returns:
        Dataset object
    """
    # Get or create DB session
    if db is None:
        db = next(get_db())
        close_db = True
    else:
        close_db = False
    
    try:
        # Convert OmicsType enum to string
        omics_type_str = omics_type.value if hasattr(omics_type, 'value') else str(omics_type)
        
        # Get or create dataset
        dataset, created = get_or_create_dataset(
            db,
            name=name,
            omics_type=omics_type_str,
            description=description,
            file_paths=file_paths or [],
            organism=organism or [],
            sample_type=sample_type or [],
            disease=disease or [],
            notion_page_id=notion_page_id,
        )
        
        action = "Created" if created else "Found existing"
        logger.info(
            "[POSTGRES][DATASET] %s dataset: %s (ID: %s)",
            action,
            dataset.name,
            dataset.id,
        )
        
        return dataset
        
    finally:
        if close_db:
            db.close()


def link_features_to_dataset_postgres(
    dataset_id: UUID,
    feature_names: List[str],
    feature_type: str,
    normalized_names: Optional[List[str]] = None,
    db: Optional[Session] = None,
) -> List[Feature]:
    """
    Link features to a dataset in Postgres.
    
    Creates features if they don't exist and links them to the dataset.
    
    Args:
        dataset_id: Dataset UUID
        feature_names: List of feature names
        feature_type: Feature type (gene, protein, metabolite, lipid)
        normalized_names: Optional list of normalized names (parallel to feature_names)
        db: Database session (optional)
    
    Returns:
        List of Feature objects
    """
    if db is None:
        db = next(get_db())
        close_db = True
    else:
        close_db = False
    
    try:
        # Prepare feature data for bulk creation
        features_data = []
        for i, name in enumerate(feature_names):
            norm_name = normalized_names[i] if normalized_names and i < len(normalized_names) else name
            features_data.append({
                'name': name,
                'feature_type': feature_type,
                'normalized_name': norm_name,
            })
        
        # Bulk get or create features
        features = bulk_get_or_create_features(db, features_data)
        
        # Link features to dataset
        feature_ids = [f.id for f in features]
        link_dataset_to_features(db, dataset_id, feature_ids)
        
        logger.info(
            "[POSTGRES][FEATURES] Linked %d %s features to dataset %s",
            len(features),
            feature_type,
            dataset_id,
        )
        
        return features
        
    finally:
        if close_db:
            db.close()


def get_postgres_session() -> Session:
    """Get a Postgres database session."""
    return next(get_db())
