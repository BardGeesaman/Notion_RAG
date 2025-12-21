"""
CRUD services for Datasets.
"""

from typing import List, Optional, Dict, Any, Sequence, cast
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.api.schemas import DatasetCreate, DatasetUpdate
from amprenta_rag.models.domain import FeatureType
import uuid


def create_dataset(db: Session, dataset: DatasetCreate) -> DatasetModel:
    """Create a new dataset."""
    db_dataset = DatasetModel(
        id=cast(Any, uuid.uuid4()),
        name=dataset.name,
        omics_type=dataset.omics_type.value,
        description=dataset.description,
        file_paths=cast(Any, dataset.file_paths if dataset.file_paths is not None else []),
        file_urls=cast(Any, dataset.file_urls if dataset.file_urls is not None else []),
        organism=cast(Any, dataset.organism if dataset.organism is not None else []),
        sample_type=cast(Any, dataset.sample_type if dataset.sample_type is not None else []),
        disease=cast(Any, dataset.disease if dataset.disease is not None else []),
    )

    # Add program relationships
    if dataset.program_ids:
        from amprenta_rag.api.services.programs import get_program
        for program_id in dataset.program_ids:
            program = get_program(db, program_id)
            if program:
                db_dataset.programs.append(program)

    # Add experiment relationships
    if dataset.experiment_ids:
        from amprenta_rag.api.services.experiments import get_experiment
        for experiment_id in dataset.experiment_ids:
            experiment = get_experiment(db, experiment_id)
            if experiment:
                db_dataset.experiments.append(experiment)

    db.add(db_dataset)
    db.commit()
    db.refresh(db_dataset)
    return db_dataset


def get_dataset(db: Session, dataset_id: UUID) -> Optional[DatasetModel]:
    """Get a dataset by ID."""
    return db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()


def get_datasets(
    db: Session,
    skip: int = 0,
    limit: int = 100,
    name_filter: Optional[str] = None,
    omics_type: Optional[str] = None,
    program_id: Optional[UUID] = None,
    experiment_id: Optional[UUID] = None,
) -> List[DatasetModel]:
    """Get all datasets with optional filtering."""
    query = db.query(DatasetModel)

    if name_filter:
        if len(name_filter) > 100:
            raise ValueError("name_filter exceeds maximum length of 100 characters")
        query = query.filter(DatasetModel.name.ilike(f"%{name_filter}%"))

    if omics_type:
        query = query.filter(DatasetModel.omics_type == omics_type)

    if program_id:
        query = query.filter(DatasetModel.programs.any(id=program_id))

    if experiment_id:
        query = query.filter(DatasetModel.experiments.any(id=experiment_id))

    return query.offset(skip).limit(limit).all()


def update_dataset(
    db: Session,
    dataset_id: UUID,
    dataset: DatasetUpdate,
) -> Optional[DatasetModel]:
    """Update a dataset."""
    db_dataset = get_dataset(db, dataset_id)
    if not db_dataset:
        return None

    update_data = dataset.model_dump(exclude_unset=True)
    program_ids = update_data.pop("program_ids", None)
    experiment_ids = update_data.pop("experiment_ids", None)

    # Handle omics_type enum conversion
    if "omics_type" in update_data and update_data["omics_type"]:
        update_data["omics_type"] = update_data["omics_type"].value

    for field, value in update_data.items():
        setattr(db_dataset, field, value)

    # Update program relationships if provided
    if program_ids is not None:
        from amprenta_rag.api.services.programs import get_program
        db_dataset.programs.clear()
        for program_id in program_ids:
            program = get_program(db, program_id)
            if program:
                db_dataset.programs.append(program)

    # Update experiment relationships if provided
    if experiment_ids is not None:
        from amprenta_rag.api.services.experiments import get_experiment
        db_dataset.experiments.clear()
        for experiment_id in experiment_ids:
            experiment = get_experiment(db, experiment_id)
            if experiment:
                db_dataset.experiments.append(experiment)

    db.commit()
    db.refresh(db_dataset)
    return db_dataset


def delete_dataset(db: Session, dataset_id: UUID) -> bool:
    """Delete a dataset."""
    db_dataset = get_dataset(db, dataset_id)
    if not db_dataset:
        return False

    db.delete(db_dataset)
    db.commit()
    return True


def get_dataset_features_by_type(
    db: Session,
    dataset_id: UUID,
) -> Dict[FeatureType, List[UUID]]:
    """Get features grouped by type for a dataset."""
    db_dataset = get_dataset(db, dataset_id)
    if not db_dataset:
        return {}

    features_by_type: Dict[FeatureType, List[UUID]] = {}

    for feature in db_dataset.features:
        try:
            feature_type = FeatureType(feature.feature_type)
            if feature_type not in features_by_type:
                features_by_type[feature_type] = []
            if feature.id is not None:
                features_by_type[feature_type].append(cast(UUID, feature.id))
        except ValueError:
            # Skip invalid feature types
            continue

    return features_by_type

