"""
Postgres-based signature creation.

Creates signatures and signature components directly in Postgres without Notion dependency.
"""

from __future__ import annotations

from typing import List, Optional, Set, Tuple
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import (
    Signature as SignatureModel,
    SignatureComponent as SignatureComponentModel,
)
from amprenta_rag.ingestion.features.postgres_linking import (
    find_or_create_feature_in_postgres,
)
from amprenta_rag.ingestion.signatures.short_id import generate_signature_short_id
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature

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
    if db is None:
        with db_session() as db:
            return _create_signature_in_postgres_impl(
                signature=signature,
                signature_type=signature_type,
                data_ownership=data_ownership,
                version=version,
                description=description,
                biomarker_roles=biomarker_roles,
                phenotype_axes=phenotype_axes,
                db=db,
            )
    else:
        return _create_signature_in_postgres_impl(
            signature=signature,
            signature_type=signature_type,
            data_ownership=data_ownership,
            version=version,
            description=description,
            biomarker_roles=biomarker_roles,
            phenotype_axes=phenotype_axes,
            db=db,
        )


def _create_signature_in_postgres_impl(
    *,
    signature: Signature,
    signature_type: str,
    data_ownership: str,
    version: Optional[str],
    description: Optional[str],
    biomarker_roles: Optional[List[str]],
    phenotype_axes: Optional[List[str]],
    db: Session,
) -> SignatureModel:
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
            setattr(existing, "biomarker_role", biomarker_roles)
            updated = True
        if phenotype_axes and not existing.phenotype_axes:
            setattr(existing, "phenotype_axes", phenotype_axes)
            updated = True
        if data_ownership and not existing.data_ownership:
            existing.data_ownership = data_ownership
            updated = True
        if signature.modalities and not existing.modalities:
            setattr(existing, "modalities", signature.modalities)
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
        modalities=signature.modalities if signature.modalities is not None else [],
        short_id=short_id,
        biomarker_role=biomarker_roles if biomarker_roles is not None else [],
        phenotype_axes=phenotype_axes if phenotype_axes is not None else [],
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
    if db is None:
        with db_session() as db:
            return _create_signature_components_in_postgres_impl(
                signature_model=signature_model,
                signature=signature,
                db=db,
            )
    return _create_signature_components_in_postgres_impl(
        signature_model=signature_model,
        signature=signature,
        db=db,
    )


def _create_signature_components_in_postgres_impl(
    *,
    signature_model: SignatureModel,
    signature: Signature,
    db: Session,
) -> Tuple[int, Set[Tuple[str, str]]]:
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
                feature_name_raw = getattr(component, "feature_name", getattr(component, "species", None)) or ""
                direction = getattr(component, "direction", None)
                weight = getattr(component, "weight", 1.0)

                if not feature_name_raw:
                    continue

                # Normalize feature name
                feature_type_enum = feature_type_map.get(feature_type_str.lower(), FeatureType.LIPID)
                normalized_name = feature_name_raw

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
                    feature_name_raw,
                    feature_type_enum.value if hasattr(feature_type_enum, "value") else feature_type_str,
                    db=db,
                )

                # Create component
                weight_value: float = float(weight) if weight is not None else 1.0
                new_component = SignatureComponentModel(
                    signature_id=signature_model.id,
                    feature_id=feature_model.id if feature_model else None,
                    feature_name=feature_name_raw,
                    feature_type=feature_type_str,
                    direction=direction,
                    weight=weight_value,
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
                feature_name = getattr(component, "feature_name", getattr(component, "species", "unknown"))
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
        db.rollback()
        return 0, set()


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
    if db is None:
        with db_session() as db:
            return _create_signature_from_file_in_postgres_impl(
                signature=signature,
                signature_type=signature_type,
                data_ownership=data_ownership,
                version=version,
                description=description,
                biomarker_roles=biomarker_roles,
                phenotype_axes=phenotype_axes,
                db=db,
            )
    return _create_signature_from_file_in_postgres_impl(
        signature=signature,
        signature_type=signature_type,
        data_ownership=data_ownership,
        version=version,
        description=description,
        biomarker_roles=biomarker_roles,
        phenotype_axes=phenotype_axes,
        db=db,
    )


def _create_signature_from_file_in_postgres_impl(
    *,
    signature: Signature,
    signature_type: str,
    data_ownership: str,
    version: Optional[str],
    description: Optional[str],
    biomarker_roles: Optional[List[str]],
    phenotype_axes: Optional[List[str]],
    db: Session,
) -> SignatureModel:
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
    if db is None:
        with db_session() as db:
            return _link_signature_to_postgres_source_impl(
                signature_id=signature_id,
                source_type=source_type,
                source_id=source_id,
                source_notion_id=source_notion_id,
                db=db,
            )
    return _link_signature_to_postgres_source_impl(
        signature_id=signature_id,
        source_type=source_type,
        source_id=source_id,
        source_notion_id=source_notion_id,
        db=db,
    )


def _link_signature_to_postgres_source_impl(
    *,
    signature_id: UUID,
    source_type: str,
    source_id: Optional[UUID],
    source_notion_id: Optional[str],
    db: Session,
) -> bool:
    try:
        # Find source by Postgres ID or Notion ID
        if source_id and source_type == "dataset":
            from amprenta_rag.database.models import Dataset as DatasetModel

            source = db.query(DatasetModel).filter(DatasetModel.id == source_id).first()
            if source:
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
        signature = db.query(SignatureModel).filter(SignatureModel.id == signature_id).first()
        if signature:
            if not getattr(signature, "external_ids", None):
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
        db.rollback()
        return False

