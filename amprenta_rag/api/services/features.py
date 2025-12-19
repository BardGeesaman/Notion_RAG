"""
CRUD services for Features.
"""

from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Feature as FeatureModel
from amprenta_rag.api.schemas import FeatureCreate, FeatureUpdate
from amprenta_rag.models.domain import FeatureType
import uuid


def create_feature(db: Session, feature: FeatureCreate) -> FeatureModel:
    """Create a new feature."""
    db_feature = FeatureModel(
        id=uuid.uuid4(),
        name=feature.name,
        feature_type=feature.feature_type.value,
        normalized_name=feature.normalized_name,
        aliases=feature.aliases or [],
        external_ids=feature.external_ids or {},
    )
    db.add(db_feature)
    db.commit()
    db.refresh(db_feature)
    return db_feature


def get_feature(db: Session, feature_id: UUID) -> Optional[FeatureModel]:
    """Get a feature by ID."""
    return db.query(FeatureModel).filter(FeatureModel.id == feature_id).first()


def get_feature_by_name(
    db: Session,
    name: str,
    feature_type: Optional[FeatureType] = None,
) -> Optional[FeatureModel]:
    """Get a feature by name (and optionally type)."""
    query = db.query(FeatureModel).filter(FeatureModel.name == name)
    if feature_type:
        query = query.filter(FeatureModel.feature_type == feature_type.value)
    return query.first()


def get_features(
    db: Session,
    skip: int = 0,
    limit: int = 100,
    name_filter: Optional[str] = None,
    feature_type: Optional[FeatureType] = None,
    dataset_id: Optional[UUID] = None,
) -> List[FeatureModel]:
    """Get all features with optional filtering."""
    query = db.query(FeatureModel)

    if name_filter:
        query = query.filter(FeatureModel.name.ilike(f"%{name_filter}%"))

    if feature_type:
        query = query.filter(FeatureModel.feature_type == feature_type.value)

    if dataset_id:
        query = query.filter(FeatureModel.datasets.any(id=dataset_id))

    return query.offset(skip).limit(limit).all()


def update_feature(
    db: Session,
    feature_id: UUID,
    feature: FeatureUpdate,
) -> Optional[FeatureModel]:
    """Update a feature."""
    db_feature = get_feature(db, feature_id)
    if not db_feature:
        return None

    update_data = feature.model_dump(exclude_unset=True)

    for field, value in update_data.items():
        setattr(db_feature, field, value)

    db.commit()
    db.refresh(db_feature)
    return db_feature


def delete_feature(db: Session, feature_id: UUID) -> bool:
    """Delete a feature."""
    db_feature = get_feature(db, feature_id)
    if not db_feature:
        return False

    db.delete(db_feature)
    db.commit()
    return True

