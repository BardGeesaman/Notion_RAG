"""
Postgres-based feature extraction for datasets.

Extracts features from Postgres datasets grouped by feature type for use in
signature matching and scoring.
"""

from __future__ import annotations

from typing import Dict, Optional, Set
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset as DatasetModel, Feature
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_dataset_features_by_type_from_postgres(
    dataset_id: UUID,
    omics_type: Optional[str] = None,
) -> Dict[str, Set[str]]:
    """
    Extract features from a Postgres dataset grouped by feature type.
    
    Uses the dataset.features relationship to get linked features from Postgres.
    Groups features by their feature_type (gene, protein, metabolite, lipid).
    
    Args:
        dataset_id: Postgres UUID of the dataset
        omics_type: Optional omics type hint (for logging/debugging)
        
    Returns:
        Dictionary mapping feature_type → set of feature names:
        {
            "gene": {"TP53", "TNF", ...},
            "protein": {"P04637", ...},
            "metabolite": {"Glutamate", ...},
            "lipid": {"Cer(d18:1/16:0)", ...}
        }
    """
    features_by_type: Dict[str, Set[str]] = {
        "gene": set(),
        "protein": set(),
        "metabolite": set(),
        "lipid": set(),
    }
    
    db = next(get_db())
    
    try:
        # Fetch dataset with features
        dataset = (
            db.query(DatasetModel)
            .filter(DatasetModel.id == dataset_id)
            .first()
        )
        
        if not dataset:
            logger.warning(
                "[FEATURE-EXTRACTION] Dataset %s not found in Postgres",
                dataset_id,
            )
            return features_by_type
        
        # Get linked features
        features = dataset.features or []
        
        logger.info(
            "[FEATURE-EXTRACTION] Found %d linked feature(s) for dataset %s",
            len(features),
            dataset_id,
        )
        
        # Group features by type
        for feature in features:
            feature_type = feature.feature_type.lower()
            feature_name = feature.name
            
            # Map feature types to our standard types
            type_mapping = {
                "gene": "gene",
                "g": "gene",
                "protein": "protein",
                "p": "protein",
                "proteomics": "protein",
                "metabolite": "metabolite",
                "m": "metabolite",
                "metabolomics": "metabolite",
                "lipid": "lipid",
                "l": "lipid",
                "lipidomics": "lipid",
                "lipid species": "lipid",
            }
            
            mapped_type = type_mapping.get(feature_type, feature_type)
            
            if mapped_type in features_by_type:
                features_by_type[mapped_type].add(feature_name)
            else:
                logger.debug(
                    "[FEATURE-EXTRACTION] Unknown feature type '%s' for feature '%s'",
                    feature_type,
                    feature_name,
                )
        
        # Log summary
        total_features = sum(len(s) for s in features_by_type.values())
        logger.info(
            "[FEATURE-EXTRACTION] Extracted %d total feature(s) from dataset %s: "
            "genes=%d, proteins=%d, metabolites=%d, lipids=%d",
            total_features,
            dataset_id,
            len(features_by_type["gene"]),
            len(features_by_type["protein"]),
            len(features_by_type["metabolite"]),
            len(features_by_type["lipid"]),
        )
        
        return features_by_type
        
    except Exception as e:
        logger.error(
            "[FEATURE-EXTRACTION] Error extracting features from dataset %s: %r",
            dataset_id,
            e,
        )
        return features_by_type
    finally:
        db.close()


def get_dataset_feature_count_by_type(
    dataset_id: UUID,
) -> Dict[str, int]:
    """
    Get count of features by type for a dataset (without loading all feature names).
    
    Useful for quick checks and summaries.
    
    Args:
        dataset_id: Postgres UUID of the dataset
        
    Returns:
        Dictionary mapping feature_type → count
    """
    from sqlalchemy import func
    
    db = next(get_db())
    
    try:
        # Query feature counts grouped by type
        from amprenta_rag.database.models import dataset_feature_assoc
        
        counts = (
            db.query(Feature.feature_type, func.count(Feature.id))
            .join(dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id)
            .filter(dataset_feature_assoc.c.dataset_id == dataset_id)
            .group_by(Feature.feature_type)
            .all()
        )
        
        result: Dict[str, int] = {
            "gene": 0,
            "protein": 0,
            "metabolite": 0,
            "lipid": 0,
        }
        
        for feature_type, count in counts:
            type_lower = feature_type.lower() if feature_type else ""
            type_mapping = {
                "gene": "gene",
                "protein": "protein",
                "metabolite": "metabolite",
                "lipid": "lipid",
            }
            mapped_type = type_mapping.get(type_lower, type_lower)
            if mapped_type in result:
                result[mapped_type] = count
        
        return result
        
    except Exception as e:
        logger.error(
            "[FEATURE-EXTRACTION] Error getting feature counts for dataset %s: %r",
            dataset_id,
            e,
        )
        return {"gene": 0, "protein": 0, "metabolite": 0, "lipid": 0}
    finally:
        db.close()

