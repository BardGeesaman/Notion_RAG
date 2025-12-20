"""
RAG builder for Postgres + Notion hybrid architecture.

Builds RAG chunks and metadata using Postgres as the source of truth for
structured data, with Notion providing narrative/documentation context.
"""

from __future__ import annotations

from typing import Dict, Optional, Any
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import (
    Program as ProgramModel,
    Experiment as ExperimentModel,
    Dataset as DatasetModel,
    Feature as FeatureModel,
    Signature as SignatureModel,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata

logger = get_logger(__name__)


def build_dataset_rag_metadata(
    dataset_id: UUID,
    db: Optional[Session] = None,
    include_notion_id: bool = True,
) -> Dict[str, Any]:
    """
    Build RAG metadata for a dataset using Postgres data.

    Args:
        dataset_id: Postgres UUID of the dataset
        db: Database session (if None, creates a new one)
        include_notion_id: Include Notion page ID if available

    Returns:
        Dictionary of metadata for Pinecone
    """
    if db is None:
        with db_session() as db:
            return build_dataset_rag_metadata(dataset_id, db=db, include_notion_id=include_notion_id)

    dataset = db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()
    if not dataset:
        logger.warning("[RAG][POSTGRES] Dataset %s not found", dataset_id)
        return {}

    metadata: Dict[str, Any] = {
        "source_type": "Dataset",
        "dataset_id": str(dataset.id),  # Postgres ID
        "dataset_name": dataset.name,
        "omics_type": dataset.omics_type,
    }

    # Include Notion ID if available and requested
    if include_notion_id and dataset.notion_page_id:
        metadata["notion_page_id"] = dataset.notion_page_id

    # Add structured data from Postgres
    if dataset.organism:
        metadata["organism"] = dataset.organism
    if dataset.sample_type:
        metadata["sample_type"] = dataset.sample_type
    if dataset.disease:
        metadata["disease"] = dataset.disease

    # Add related program/experiment IDs
    if dataset.programs:
        metadata["program_ids"] = [str(p.id) for p in dataset.programs]
    if dataset.experiments:
        metadata["experiment_ids"] = [str(e.id) for e in dataset.experiments]

    # Add signature match score if available
    if dataset.signature_match_score is not None:
        metadata["signature_match_score"] = dataset.signature_match_score

    return sanitize_metadata(metadata)


def build_program_rag_metadata(
    program_id: UUID,
    db: Optional[Session] = None,
    include_notion_id: bool = True,
) -> Dict[str, Any]:
    """Build RAG metadata for a program using Postgres data."""
    if db is None:
        with db_session() as db:
            return build_program_rag_metadata(program_id, db=db, include_notion_id=include_notion_id)

    program = db.query(ProgramModel).filter(ProgramModel.id == program_id).first()
    if not program:
        logger.warning("[RAG][POSTGRES] Program %s not found", program_id)
        return {}

    metadata: Dict[str, Any] = {
        "source_type": "Program",
        "program_id": str(program.id),
        "program_name": program.name,
    }

    if include_notion_id and program.notion_page_id:
        metadata["notion_page_id"] = program.notion_page_id

    if program.disease:
        metadata["disease"] = program.disease

    return sanitize_metadata(metadata)


def build_experiment_rag_metadata(
    experiment_id: UUID,
    db: Optional[Session] = None,
    include_notion_id: bool = True,
) -> Dict[str, Any]:
    """Build RAG metadata for an experiment using Postgres data."""
    if db is None:
        with db_session() as db:
            return build_experiment_rag_metadata(experiment_id, db=db, include_notion_id=include_notion_id)

    experiment = db.query(ExperimentModel).filter(ExperimentModel.id == experiment_id).first()
    if not experiment:
        logger.warning("[RAG][POSTGRES] Experiment %s not found", experiment_id)
        return {}

    metadata: Dict[str, Any] = {
        "source_type": "Experiment",
        "experiment_id": str(experiment.id),
        "experiment_name": experiment.name,
    }

    if include_notion_id and experiment.notion_page_id:
        metadata["notion_page_id"] = experiment.notion_page_id

    if experiment.type:
        metadata["experiment_type"] = experiment.type
    if experiment.disease:
        metadata["disease"] = experiment.disease
    if experiment.matrix:
        metadata["matrix"] = experiment.matrix
    if experiment.model_systems:
        metadata["model_systems"] = experiment.model_systems

    if experiment.programs:
        metadata["program_ids"] = [str(p.id) for p in experiment.programs]

    return sanitize_metadata(metadata)


def build_signature_rag_metadata(
    signature_id: UUID,
    db: Optional[Session] = None,
    include_notion_id: bool = True,
) -> Dict[str, Any]:
    """Build RAG metadata for a signature using Postgres data."""
    if db is None:
        with db_session() as db:
            return build_signature_rag_metadata(signature_id, db=db, include_notion_id=include_notion_id)

    signature = db.query(SignatureModel).filter(SignatureModel.id == signature_id).first()
    if not signature:
        logger.warning("[RAG][POSTGRES] Signature %s not found", signature_id)
        return {}

    metadata: Dict[str, Any] = {
        "source_type": "Signature",
        "signature_id": str(signature.id),
        "signature_name": signature.name,
    }

    if include_notion_id and signature.notion_page_id:
        metadata["notion_page_id"] = signature.notion_page_id

    if signature.modalities:
        metadata["modalities"] = signature.modalities

    if signature.programs:
        metadata["program_ids"] = [str(p.id) for p in signature.programs]

    return sanitize_metadata(metadata)


def build_feature_rag_metadata(
    feature_id: UUID,
    db: Optional[Session] = None,
    include_notion_id: bool = True,
) -> Dict[str, Any]:
    """Build RAG metadata for a feature using Postgres data."""
    if db is None:
        with db_session() as db:
            return build_feature_rag_metadata(feature_id, db=db, include_notion_id=include_notion_id)

    feature = db.query(FeatureModel).filter(FeatureModel.id == feature_id).first()
    if not feature:
        logger.warning("[RAG][POSTGRES] Feature %s not found", feature_id)
        return {}

    metadata: Dict[str, Any] = {
        "source_type": "Feature",
        "feature_id": str(feature.id),
        "feature_name": feature.name,
        "feature_type": feature.feature_type,
    }

    if include_notion_id and feature.notion_page_id:
        metadata["notion_page_id"] = feature.notion_page_id

    return sanitize_metadata(metadata)


def build_dataset_rag_text(
    dataset_id: UUID,
    db: Optional[Session] = None,
) -> str:
    """
    Build RAG text representation for a dataset.

    Args:
        dataset_id: Postgres UUID of the dataset
        db: Database session

    Returns:
        Text representation for embedding
    """
    if db is None:
        with db_session() as db:
            return build_dataset_rag_text(dataset_id, db=db)

    dataset = db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()
    if not dataset:
        return ""

    parts = []

    # Structured data from Postgres
    parts.append(f"Dataset: {dataset.name}")
    if dataset.description:
        parts.append(f"Description: {dataset.description}")
    parts.append(f"Omics Type: {dataset.omics_type}")

    if dataset.organism:
        parts.append(f"Organism: {', '.join(dataset.organism)}")
    if dataset.sample_type:
        parts.append(f"Sample Type: {', '.join(dataset.sample_type)}")
    if dataset.disease:
        parts.append(f"Disease: {', '.join(dataset.disease)}")

    # Related programs
    if dataset.programs:
        program_names = [p.name for p in dataset.programs if p.name]
        parts.append(f"Programs: {', '.join(program_names)}")

    # Related experiments
    if dataset.experiments:
        experiment_names = [e.name for e in dataset.experiments if e.name]
        parts.append(f"Experiments: {', '.join(experiment_names)}")

    # Signature match score
    if dataset.signature_match_score is not None:
        parts.append(f"Signature Match Score: {dataset.signature_match_score:.3f}")

    return "\n".join(parts)


def build_program_rag_text(
    program_id: UUID,
    db: Optional[Session] = None,
) -> str:
    """Build RAG text representation for a program."""
    if db is None:
        with db_session() as db:
            return build_program_rag_text(program_id, db=db)

    program = db.query(ProgramModel).filter(ProgramModel.id == program_id).first()
    if not program:
        return ""

    parts = []
    parts.append(f"Program: {program.name}")
    if program.description:
        parts.append(f"Description: {program.description}")
    if program.disease:
        diseases = [d for d in program.disease if d]
        parts.append(f"Disease Focus: {', '.join(diseases)}")

    # Get related datasets count
    dataset_count = len(program.datasets) if program.datasets else 0
    if dataset_count > 0:
        parts.append(f"Related Datasets: {dataset_count}")

    return "\n".join(parts)

