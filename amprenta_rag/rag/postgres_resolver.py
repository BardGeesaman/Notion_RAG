"""
Postgres ID resolver for RAG queries.

Resolves Postgres IDs from Pinecone metadata and fetches data from Postgres
for RAG context retrieval.
"""

from __future__ import annotations

from typing import Dict, Optional, Any
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import (
    Program as ProgramModel,
    Experiment as ExperimentModel,
    Dataset as DatasetModel,
    Feature as FeatureModel,
    Signature as SignatureModel,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.rag.postgres_builder import (
    build_dataset_rag_text,
    build_program_rag_text,
)

logger = get_logger(__name__)


def resolve_postgres_id_from_metadata(metadata: Dict[str, Any]) -> Optional[UUID]:
    """
    Resolve Postgres ID from Pinecone metadata.
    
    Checks for dataset_id, program_id, experiment_id, signature_id, or feature_id.
    
    Args:
        metadata: Pinecone metadata dictionary
        
    Returns:
        Postgres UUID if found, None otherwise
    """
    # Check for various ID fields
    for id_field in ["dataset_id", "program_id", "experiment_id", "signature_id", "feature_id"]:
        if id_field in metadata:
            try:
                return UUID(metadata[id_field])
            except (ValueError, TypeError):
                continue
    
    return None


def get_entity_type_from_metadata(metadata: Dict[str, Any]) -> Optional[str]:
    """
    Determine entity type from metadata.
    
    Args:
        metadata: Pinecone metadata dictionary
        
    Returns:
        Entity type ("Dataset", "Program", "Experiment", "Signature", "Feature") or None
    """
    source_type = metadata.get("source_type") or metadata.get("source")
    if source_type:
        return source_type
    
    # Infer from ID fields
    if "dataset_id" in metadata:
        return "Dataset"
    if "program_id" in metadata:
        return "Program"
    if "experiment_id" in metadata:
        return "Experiment"
    if "signature_id" in metadata:
        return "Signature"
    if "feature_id" in metadata:
        return "Feature"
    
    return None


def fetch_postgres_context(
    postgres_id: UUID,
    entity_type: str,
    db: Optional[Session] = None,
    include_notion_narrative: bool = True,
) -> Optional[str]:
    """
    Fetch RAG context text from Postgres for a given entity.
    
    Args:
        postgres_id: Postgres UUID
        entity_type: Type of entity ("Dataset", "Program", etc.)
        db: Database session
        
    Returns:
        Text representation for RAG, or None if not found
    """
    if db is None:
        db = next(get_db())
    
    try:
        if entity_type == "Dataset":
            return build_dataset_rag_text(postgres_id, db=db, include_notion_narrative=include_notion_narrative)
        elif entity_type == "Program":
            return build_program_rag_text(postgres_id, db=db, include_notion_narrative=include_notion_narrative)
        elif entity_type == "Experiment":
            # For experiments, build a simple text representation
            experiment = db.query(ExperimentModel).filter(ExperimentModel.id == postgres_id).first()
            if experiment:
                parts = [f"Experiment: {experiment.name}"]
                if experiment.description:
                    parts.append(f"Description: {experiment.description}")
                if experiment.disease:
                    parts.append(f"Disease: {', '.join(experiment.disease)}")
                return "\n".join(parts)
        elif entity_type == "Signature":
            signature = db.query(SignatureModel).filter(SignatureModel.id == postgres_id).first()
            if signature:
                parts = [f"Signature: {signature.name}"]
                if signature.description:
                    parts.append(f"Description: {signature.description}")
                if signature.modalities:
                    parts.append(f"Modalities: {', '.join(signature.modalities)}")
                return "\n".join(parts)
        elif entity_type == "Feature":
            feature = db.query(FeatureModel).filter(FeatureModel.id == postgres_id).first()
            if feature:
                return f"Feature: {feature.name} (Type: {feature.feature_type})"
        
        return None
        
    except Exception as e:
        logger.warning(
            "[RAG][POSTGRES] Error fetching context for %s %s: %r",
            entity_type,
            postgres_id,
            e,
        )
        return None


def get_notion_id_from_postgres(
    postgres_id: UUID,
    entity_type: str,
    db: Optional[Session] = None,
) -> Optional[str]:
    """
    Get Notion page ID from Postgres entity.
    
    Useful for hybrid queries that need to fetch narrative from Notion.
    
    Args:
        postgres_id: Postgres UUID
        entity_type: Type of entity
        db: Database session
        
    Returns:
        Notion page ID if available, None otherwise
    """
    if db is None:
        db = next(get_db())
    
    try:
        if entity_type == "Dataset":
            dataset = db.query(DatasetModel).filter(DatasetModel.id == postgres_id).first()
            return dataset.notion_page_id if dataset else None
        elif entity_type == "Program":
            program = db.query(ProgramModel).filter(ProgramModel.id == postgres_id).first()
            return program.notion_page_id if program else None
        elif entity_type == "Experiment":
            experiment = db.query(ExperimentModel).filter(ExperimentModel.id == postgres_id).first()
            return experiment.notion_page_id if experiment else None
        elif entity_type == "Signature":
            signature = db.query(SignatureModel).filter(SignatureModel.id == postgres_id).first()
            return signature.notion_page_id if signature else None
        elif entity_type == "Feature":
            feature = db.query(FeatureModel).filter(FeatureModel.id == postgres_id).first()
            return feature.notion_page_id if feature else None
        
        return None
        
    except Exception as e:
        logger.warning(
            "[RAG][POSTGRES] Error getting Notion ID for %s %s: %r",
            entity_type,
            postgres_id,
            e,
        )
        return None

