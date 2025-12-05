"""
Postgres-based signature creation.

Creates signatures and signature components directly in Postgres without Notion dependency.
"""

from __future__ import annotations

from typing import List, Optional, Set, Tuple
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import (
    Signature as SignatureModel,
    SignatureComponent as SignatureComponentModel,
    Feature as FeatureModel,
)
from amprenta_rag.ingestion.features.postgres_linking import (
    find_or_create_feature_in_postgres,
)
from amprenta_rag.ingestion.signatures.short_id import generate_signature_short_id
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature, SignatureComponent

logger = get_logger(__name__)


def create_signature_in_postgres(
    signature: Signature,
    signature_type: str = "Literature-derived",
    data_ownership: str = "Public",
    version: Optional[str] = None,
    description: Optional[str] = None,
    biomarker_roles: Optional[List[str]] = None,
    phenotype_axes: Optional[List[str]] = None,
    db: Optional[Session] = None,
) -> SignatureModel:
    """
    Create a signature in Postgres or return existing one.
    
    Args:
        signature: Signature object to create
        signature_type: Type of signature (e.g., "Literature-derived")
        data_ownership: Data ownership classification
        version: Optional version string
        description: Optional description
        biomarker_roles: Optional biomarker roles
        phenotype_axes: Optional phenotype axes
        db: Optional database session
        
    Returns:
        SignatureModel instance
    """
    should_close_db = False
    if db is None:
        db = next(get_db())
        should_close_db = True
    
    try:
        # Generate short ID
        short_id = generate_signature_short_id(signature.name, version)
        
        # Check if signature already exists by short_id or name
        existing = (
            db.query(SignatureModel)
            .filter(
                (SignatureModel.short_id == short_id) | (SignatureModel.name == signature.name)
            )
            .first()
        )
        
        if existing:
            logger.debug(
                "[POSTGRES-SIGNATURE] Found existing signature: %s (ID: %s)",
                existing.name,
                existing.id,
            )
            
            # Update metadata if provided
            updated = False
            if description and not existing.description:
                existing.description = description
                updated = True
            if biomarker_roles and not existing.biomarker_role:
                existing.biomarker_role = biomarker_roles
                updated = True
            if phenotype_axes and not existing.phenotype_axes:
                existing.phenotype_axes = phenotype_axes
                updated = True
            if data_ownership and not existing.data_ownership:
                existing.data_ownership = data_ownership
                updated = True
            if signature.modalities and not existing.modalities:
                existing.modalities = signature.modalities
                updated = True
            if not existing.short_id:
                existing.short_id = short_id
                updated = True
            
            if updated:
                db.commit()
                logger.debug("[POSTGRES-SIGNATURE] Updated signature metadata: %s", existing.name)
            
            return existing
        
        # Create new signature
        new_signature = SignatureModel(
            name=signature.name,
            description=description or signature.description,
            modalities=signature.modalities or [],
            short_id=short_id,
            biomarker_role=biomarker_roles or [],
            phenotype_axes=phenotype_axes or [],
            data_ownership=data_ownership,
        )
        
        db.add(new_signature)
        db.commit()
        db.refresh(new_signature)
        
        logger.info(
            "[POSTGRES-SIGNATURE] Created signature: %s (ID: %s, Short ID: %s)",
            new_signature.name,
            new_signature.id,
            new_signature.short_id,
        )
        
        return new_signature
        
    except Exception as e:
        logger.error(
            "[POSTGRES-SIGNATURE] Error creating signature %s: %r",
            signature.name,
            e,
        )
        if db:
            db.rollback()
        raise
    finally:
        if should_close_db and db:
            db.close()


def create_signature_components_in_postgres(
    signature_model: SignatureModel,
    signature: Signature,
    db: Optional[Session] = None,
) -> Tuple[int, Set[Tuple[str, str]]]:
    """
    Create signature components in Postgres and link to features.
    
    Args:
        signature_model: SignatureModel instance
        signature: Signature object with components
        db: Optional database session
        
    Returns:
        Tuple of (component_count, features_created_set)
        features_created_set contains (feature_type, feature_name) tuples
    """
    should_close_db = False
    if db is None:
        db = next(get_db())
        should_close_db = True
    
    component_count = 0
    features_created: Set[Tuple[str, str]] = set()
    
    try:
        # Check existing components
        existing_components = (
            db.query(SignatureComponentModel)
            .filter(SignatureComponentModel.signature_id == signature_model.id)
            .all()
        )
        
        # Build set of existing component keys
        existing_keys = {
            (comp.feature_name or "", comp.feature_type, comp.direction or "")
            for comp in existing_components
        }
        
        # Import normalization functions
        from amprenta_rag.ingestion.features.postgres_linking import (
            normalize_feature_name,
        )
        from amprenta_rag.models.domain import FeatureType
        
        # Map feature types
        feature_type_map = {
            "gene": FeatureType.GENE,
            "protein": FeatureType.PROTEIN,
            "metabolite": FeatureType.METABOLITE,
            "lipid": FeatureType.LIPID,
        }
        
        for component in signature.components:
            try:
                feature_type_str = getattr(component, "feature_type", "lipid") or "lipid"
                feature_name_raw = getattr(component, "feature_name", component.species) or ""
                direction = getattr(component, "direction", None)
                weight = getattr(component, "weight", 1.0)
                
                if not feature_name_raw:
                    continue
                
                # Normalize feature name
                feature_type_enum = feature_type_map.get(feature_type_str.lower(), FeatureType.LIPID)
                normalized_name = normalize_feature_name(feature_name_raw, feature_type_enum)
                
                # Check if component already exists
                component_key = (feature_name_raw, feature_type_str, direction or "")
                if component_key in existing_keys:
                    logger.debug(
                        "[POSTGRES-SIGNATURE] Component already exists: %s (%s)",
                        feature_name_raw,
                        feature_type_str,
                    )
                    continue
                
                # Find or create feature
                feature_model = find_or_create_feature_in_postgres(
                    name=feature_name_raw,
                    feature_type=feature_type_enum,
                    normalized_name=normalized_name,
                    db=db,
                )
                
                # Create component
                new_component = SignatureComponentModel(
                    signature_id=signature_model.id,
                    feature_id=feature_model.id if feature_model else None,
                    feature_name=feature_name_raw,
                    feature_type=feature_type_str,
                    direction=direction,
                    weight=float(weight) if weight else 1.0,
                )
                
                db.add(new_component)
                component_count += 1
                
                features_created.add((feature_type_str, normalized_name or feature_name_raw))
                
                logger.debug(
                    "[POSTGRES-SIGNATURE] Created component: %s (%s) -> %s",
                    feature_name_raw,
                    feature_type_str,
                    signature_model.name,
                )
                
            except Exception as e:
                feature_name = getattr(component, "feature_name", component.species)
                logger.warning(
                    "[POSTGRES-SIGNATURE] Error creating component %s: %r",
                    feature_name,
                    e,
                )
                continue
        
        db.commit()
        
        logger.info(
            "[POSTGRES-SIGNATURE] Created %d component(s) for signature %s",
            component_count,
            signature_model.name,
        )
        
        return component_count, features_created
        
    except Exception as e:
        logger.error(
            "[POSTGRES-SIGNATURE] Error creating components for signature %s: %r",
            signature_model.id,
            e,
        )
        if db:
            db.rollback()
        return 0, set()
    finally:
        if should_close_db and db:
            db.close()


def create_signature_from_file_in_postgres(
    signature: Signature,
    signature_type: str = "Literature-derived",
    data_ownership: str = "Public",
    version: Optional[str] = None,
    description: Optional[str] = None,
    biomarker_roles: Optional[List[str]] = None,
    phenotype_axes: Optional[List[str]] = None,
    db: Optional[Session] = None,
) -> SignatureModel:
    """
    Create a complete signature in Postgres from a Signature object.
    
    Creates the signature and all its components, linking to features.
    
    Args:
        signature: Signature object from loader
        signature_type: Type of signature
        data_ownership: Data ownership
        version: Optional version
        description: Optional description
        biomarker_roles: Optional biomarker roles
        phenotype_axes: Optional phenotype axes
        db: Optional database session
        
    Returns:
        SignatureModel instance
    """
    should_close_db = False
    if db is None:
        db = next(get_db())
        should_close_db = True
    
    try:
        # Create signature
        signature_model = create_signature_in_postgres(
            signature=signature,
            signature_type=signature_type,
            data_ownership=data_ownership,
            version=version,
            description=description,
            biomarker_roles=biomarker_roles,
            phenotype_axes=phenotype_axes,
            db=db,
        )
        
        # Create components
        component_count, features_created = create_signature_components_in_postgres(
            signature_model=signature_model,
            signature=signature,
            db=db,
        )
        
        logger.info(
            "[POSTGRES-SIGNATURE] Created signature %s with %d component(s) and %d feature(s)",
            signature_model.name,
            component_count,
            len(features_created),
        )
        
        return signature_model
        
    finally:
        if should_close_db and db:
            db.close()


def link_signature_to_postgres_source(
    signature_id: UUID,
    source_type: str,
    source_id: Optional[UUID] = None,
    source_notion_id: Optional[str] = None,
    db: Optional[Session] = None,
) -> bool:
    """
    Link a signature to a source (dataset, experiment, etc.) in Postgres.
    
    Args:
        signature_id: Postgres UUID of the signature
        source_type: Type of source ("dataset", "experiment")
        source_id: Postgres UUID of the source (preferred)
        source_notion_id: Notion page ID of the source (fallback for backward compat)
        db: Optional database session
        
    Returns:
        True if linking succeeded, False otherwise
    """
    should_close_db = False
    if db is None:
        db = next(get_db())
        should_close_db = True
    
    try:
        # Find source by Postgres ID or Notion ID
        if source_id and source_type == "dataset":
            from amprenta_rag.database.models import Dataset as DatasetModel
            
            source = db.query(DatasetModel).filter(DatasetModel.id == source_id).first()
            if source:
                # Link via relationship
                signature = db.query(SignatureModel).filter(SignatureModel.id == signature_id).first()
                if signature:
                    source.signatures.append(signature)
                    db.commit()
                    logger.info(
                        "[POSTGRES-SIGNATURE] Linked signature %s to dataset %s",
                        signature_id,
                        source_id,
                    )
                    return True
        
        # Store source reference in signature's external_ids for now
        # TODO: Create proper source_signature_assoc table if needed
        signature = db.query(SignatureModel).filter(SignatureModel.id == signature_id).first()
        if signature:
            if not signature.external_ids:
                signature.external_ids = {}
            
            source_key = f"{source_type}_source"
            if source_notion_id:
                signature.external_ids[source_key] = source_notion_id
            elif source_id:
                signature.external_ids[f"{source_key}_uuid"] = str(source_id)
            
            db.commit()
            logger.debug(
                "[POSTGRES-SIGNATURE] Stored source reference for signature %s",
                signature_id,
            )
            return True
        
        return False
        
    except Exception as e:
        logger.error(
            "[POSTGRES-SIGNATURE] Error linking signature %s to source: %r",
            signature_id,
            e,
        )
        if db:
            db.rollback()
        return False
    finally:
        if should_close_db and db:
            db.close()

