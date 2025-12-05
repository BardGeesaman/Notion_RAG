"""
Postgres-based feature linking.

This module provides functions to create/find features in Postgres and link them
to datasets. This replaces or complements Notion-based feature linking for
better performance and scalability.

Features:
- Create/find features in Postgres
- Link features to datasets via association table
- Batch operations for performance
- Automatic normalization
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple
from uuid import UUID, uuid4

from sqlalchemy.orm import Session
from sqlalchemy import and_, or_

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Feature, Dataset, dataset_feature_assoc
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import FeatureType

logger = get_logger(__name__)


def normalize_feature_name(
    feature_name: str,
    feature_type: FeatureType,
) -> str:
    """
    Normalize feature name based on feature type.
    
    This function delegates to the appropriate normalization function
    for each feature type.
    
    Args:
        feature_name: Raw feature name
        feature_type: Type of feature (gene, protein, metabolite, lipid)
        
    Returns:
        Normalized feature name
    """
    if feature_type == FeatureType.GENE:
        from amprenta_rag.ingestion.transcriptomics.normalization import normalize_gene_identifier
        return normalize_gene_identifier(feature_name)
    elif feature_type == FeatureType.PROTEIN:
        from amprenta_rag.ingestion.proteomics.normalization import normalize_protein_identifier
        return normalize_protein_identifier(feature_name)
    elif feature_type == FeatureType.METABOLITE:
        from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name
        return normalize_metabolite_name(feature_name)
    elif feature_type == FeatureType.LIPID:
        from amprenta_rag.ingestion.signature_matching.species_mapping import map_raw_lipid_to_canonical_species
        result = map_raw_lipid_to_canonical_species(feature_name)
        return result if result else feature_name.strip().lower()
    else:
        # Default: lowercase, strip whitespace
        return feature_name.strip().lower()


def find_or_create_feature_in_postgres(
    name: str,
    feature_type: FeatureType,
    normalized_name: Optional[str] = None,
    db: Optional[Session] = None,
) -> Feature:
    """
    Find or create a feature in Postgres.
    
    Searches by normalized_name first, then by name. If not found, creates
    a new feature record.
    
    Args:
        name: Feature name (raw)
        feature_type: Type of feature
        normalized_name: Optional pre-normalized name (if None, will normalize)
        db: Optional database session (if None, creates new session)
        
    Returns:
        Feature model instance
    """
    # Normalize if not provided
    if normalized_name is None:
        normalized_name = normalize_feature_name(name, feature_type)
    
    # Use provided session or create new one
    if db is None:
        db = next(get_db())
        should_close = True
    else:
        should_close = False
    
    try:
        # Search by normalized name first (preferred)
        feature = db.query(Feature).filter(
            and_(
                Feature.feature_type == feature_type.value,
                or_(
                    Feature.normalized_name == normalized_name,
                    Feature.name == normalized_name,
                ),
            ),
        ).first()
        
        if feature:
            logger.debug(
                "[FEATURE][POSTGRES] Found existing feature: %s (type: %s, id: %s)",
                normalized_name,
                feature_type.value,
                feature.id,
            )
            return feature
        
        # Not found, create new feature
        feature = Feature(
            id=uuid4(),
            name=name,
            feature_type=feature_type.value,
            normalized_name=normalized_name,
        )
        
        db.add(feature)
        db.commit()
        db.refresh(feature)
        
        logger.info(
            "[FEATURE][POSTGRES] Created new feature: %s (type: %s, id: %s)",
            normalized_name,
            feature_type.value,
            feature.id,
        )
        
        return feature
        
    except Exception as e:
        logger.error(
            "[FEATURE][POSTGRES] Error finding/creating feature %s (type: %s): %r",
            normalized_name,
            feature_type.value,
            e,
        )
        if should_close:
            db.rollback()
        raise
    finally:
        if should_close:
            db.close()


def link_feature_to_dataset_in_postgres(
    feature_id: UUID,
    dataset_id: UUID,
    db: Optional[Session] = None,
) -> bool:
    """
    Link a feature to a dataset in Postgres via association table.
    
    Args:
        feature_id: Feature UUID
        dataset_id: Dataset UUID
        db: Optional database session
        
    Returns:
        True if linked (or already linked), False on error
    """
    if db is None:
        db = next(get_db())
        should_close = True
    else:
        should_close = False
    
    try:
        # Check if link already exists
        existing_link = db.execute(
            dataset_feature_assoc.select().where(
                and_(
                    dataset_feature_assoc.c.feature_id == feature_id,
                    dataset_feature_assoc.c.dataset_id == dataset_id,
                ),
            ),
        ).first()
        
        if existing_link:
            logger.debug(
                "[FEATURE][POSTGRES] Feature %s already linked to dataset %s",
                feature_id,
                dataset_id,
            )
            return True
        
        # Create link via association table
        db.execute(
            dataset_feature_assoc.insert().values(
                feature_id=feature_id,
                dataset_id=dataset_id,
            ),
        )
        db.commit()
        
        logger.debug(
            "[FEATURE][POSTGRES] Linked feature %s to dataset %s",
            feature_id,
            dataset_id,
        )
        
        return True
        
    except Exception as e:
        logger.warning(
            "[FEATURE][POSTGRES] Error linking feature %s to dataset %s: %r",
            feature_id,
            dataset_id,
            e,
        )
        if should_close:
            db.rollback()
        return False
    finally:
        if should_close:
            db.close()


def batch_link_features_to_dataset_in_postgres(
    features: List[Tuple[str, FeatureType]],
    dataset_id: UUID,
    max_workers: int = 10,
    db: Optional[Session] = None,
) -> Dict[str, bool]:
    """
    Batch link multiple features to a dataset in Postgres.
    
    Args:
        features: List of (feature_name, feature_type) tuples
        dataset_id: Dataset UUID to link to
        max_workers: Maximum parallel workers (for feature creation)
        db: Optional database session
        
    Returns:
        Dictionary mapping feature_name -> success (True/False)
    """
    logger.info(
        "[FEATURE][POSTGRES][BATCH] Batch linking %d features to dataset %s",
        len(features),
        dataset_id,
    )
    
    results: Dict[str, bool] = {}
    
    # Use provided session or create new one
    if db is None:
        db = next(get_db())
        should_close = True
    else:
        should_close = False
    
    try:
        # Process features
        for feature_name, feature_type in features:
            try:
                # Find or create feature
                feature = find_or_create_feature_in_postgres(
                    name=feature_name,
                    feature_type=feature_type,
                    db=db,
                )
                
                # Link to dataset
                success = link_feature_to_dataset_in_postgres(
                    feature_id=feature.id,
                    dataset_id=dataset_id,
                    db=db,
                )
                
                results[feature_name] = success
                
            except Exception as e:
                logger.warning(
                    "[FEATURE][POSTGRES][BATCH] Error processing feature %s: %r",
                    feature_name,
                    e,
                )
                results[feature_name] = False
        
        logger.info(
            "[FEATURE][POSTGRES][BATCH] Linked %d/%d features to dataset %s",
            sum(1 for v in results.values() if v),
            len(features),
            dataset_id,
        )
        
        return results
        
    finally:
        if should_close and db:
            db.close()


def get_dataset_features_from_postgres(
    dataset_id: UUID,
    feature_type: Optional[FeatureType] = None,
    db: Optional[Session] = None,
) -> List[Feature]:
    """
    Get all features linked to a dataset from Postgres.
    
    Args:
        dataset_id: Dataset UUID
        feature_type: Optional filter by feature type
        db: Optional database session
        
    Returns:
        List of Feature model instances
    """
    if db is None:
        db = next(get_db())
        should_close = True
    else:
        should_close = False
    
    try:
        query = db.query(Feature).join(
            dataset_feature_assoc,
            Feature.id == dataset_feature_assoc.c.feature_id,
        ).filter(
            dataset_feature_assoc.c.dataset_id == dataset_id,
        )
        
        if feature_type:
            query = query.filter(Feature.feature_type == feature_type.value)
        
        features = query.all()
        
        return features
        
    finally:
        if should_close and db:
            db.close()

