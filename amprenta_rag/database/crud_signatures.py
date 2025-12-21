"""
Signature CRUD operations (and related linking helpers).

Split out of `amprenta_rag.database.crud` to keep domain operations smaller and
more maintainable. The facade module re-exports these functions for backwards
compatibility.
"""

from __future__ import annotations

from typing import Any, List, Optional, cast
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset, Signature, SignatureComponent
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_signature(
    db: Session,
    name: str,
    description: Optional[str] = None,
    modalities: Optional[List[str]] = None,
    short_id: Optional[str] = None,
    biomarker_role: Optional[List[str]] = None,
    phenotype_axes: Optional[List[str]] = None,
    data_ownership: Optional[str] = None,
    notion_page_id: Optional[str] = None,
    commit: bool = True,
) -> Signature:
    """Create a new signature in the database."""
    modalities_list: List[str] = list(modalities) if modalities else []
    biomarker_list: List[str] = list(biomarker_role) if biomarker_role else []
    phenotype_list: List[str] = list(phenotype_axes) if phenotype_axes else []
    signature = Signature(
        name=name,
        description=description,
        modalities=cast(Any, modalities_list),
        short_id=short_id,
        biomarker_role=cast(Any, biomarker_list),
        phenotype_axes=cast(Any, phenotype_list),
        data_ownership=data_ownership,
        notion_page_id=notion_page_id,
    )

    db.add(signature)

    if commit:
        db.commit()
        db.refresh(signature)
        logger.info(
            "[CRUD][SIGNATURE] Created signature: %s (ID: %s)",
            signature.name,
            signature.id,
        )

    return signature


def create_signature_component(
    db: Session,
    signature_id: UUID,
    feature_id: Optional[UUID],
    feature_name: str,
    feature_type: str,
    direction: Optional[str] = None,
    weight: float = 1.0,
    commit: bool = True,
) -> SignatureComponent:
    """Create a signature component."""
    weight_value: float | None = float(weight) if weight is not None else None
    component = SignatureComponent(
        signature_id=cast(Any, signature_id),
        feature_id=cast(Any, feature_id),
        feature_name=feature_name,
        feature_type=feature_type,
        direction=direction,
        weight=cast(Any, weight_value),
    )

    db.add(component)

    if commit:
        db.commit()
        db.refresh(component)

    return component


def get_signature_by_id(db: Session, signature_id: UUID) -> Optional[Signature]:
    """Get signature by UUID."""
    return db.query(Signature).filter(Signature.id == signature_id).first()


def get_signature_by_name(db: Session, name: str) -> Optional[Signature]:
    """Get signature by name."""
    return db.query(Signature).filter(Signature.name == name).first()


def get_or_create_signature(
    db: Session,
    name: str,
    **kwargs,
) -> tuple[Signature, bool]:
    """
    Get existing signature or create new one.

    Returns:
        Tuple of (Signature, created: bool)
    """
    notion_page_id = kwargs.get("notion_page_id")

    # Try to find by notion_page_id first
    if notion_page_id:
        signature = db.query(Signature).filter(Signature.notion_page_id == notion_page_id).first()
        if signature:
            return signature, False

    # Try to find by name
    signature = get_signature_by_name(db, name)
    if signature:
        return signature, False

    # Create new
    signature = create_signature(db, name=name, **kwargs)
    return signature, True


def link_dataset_to_signature(
    db: Session,
    dataset_id: UUID,
    signature_id: UUID,
    match_score: Optional[float] = None,
    commit: bool = True,
) -> None:
    """Link a dataset to a signature with an optional match score."""
    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    signature = get_signature_by_id(db, signature_id)

    if not dataset or not signature:
        logger.error(
            "[CRUD][LINK] Cannot link dataset %s to signature %s: not found",
            dataset_id,
            signature_id,
        )
        return

    if signature not in dataset.signatures:
        dataset.signatures.append(signature)

        # If match_score is provided, update the junction table
        # This requires direct SQL manipulation - simplified version for now
        if match_score is not None and hasattr(dataset, "signature_match_score"):
            setattr(dataset, "signature_match_score", float(match_score))

        if commit:
            db.commit()
            logger.debug(
                "[CRUD][LINK] Linked dataset %s to signature %s (score: %.3f)",
                dataset.name,
                signature.name,
                match_score or 0.0,
            )


