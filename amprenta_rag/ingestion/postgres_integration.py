"""
Postgres integration for ingestion pipelines.

Provides utilities to integrate Postgres as source of truth into ingestion
pipelines, with optional dual-write to Notion during transition.
"""

from __future__ import annotations

from typing import Dict, List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.config import get_config
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.migration.dual_write import DualWriteManager
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)


def create_or_update_dataset_in_postgres(
    name: str,
    omics_type: OmicsType,
    file_paths: Optional[List[str]] = None,
    file_urls: Optional[List[str]] = None,
    description: Optional[str] = None,
    organism: Optional[List[str]] = None,
    sample_type: Optional[List[str]] = None,
    disease: Optional[List[str]] = None,
    program_ids: Optional[List[UUID]] = None,
    experiment_ids: Optional[List[UUID]] = None,
    notion_page_id: Optional[str] = None,
    external_ids: Optional[Dict[str, str]] = None,
    db: Optional[Session] = None,
) -> DatasetModel:
    """
    Create or update a dataset in Postgres.
    
    If notion_page_id is provided and a dataset with that notion_page_id exists,
    it will be updated. Otherwise, a new dataset will be created.
    
    Args:
        name: Dataset name
        omics_type: Omics type
        file_paths: List of file paths
        file_urls: List of file URLs
        description: Dataset description
        organism: List of organisms
        sample_type: List of sample types
        disease: List of diseases
        program_ids: List of program UUIDs
        experiment_ids: List of experiment UUIDs
        notion_page_id: Notion page ID (for linking/updating)
        external_ids: External identifier mappings
        db: Database session (if None, creates a new one)
        
    Returns:
        DatasetModel instance
    """
    if db is None:
        db = next(get_db())
    
    # Check if dataset exists by notion_page_id
    existing_dataset = None
    if notion_page_id:
        existing_dataset = db.query(DatasetModel).filter(
            DatasetModel.notion_page_id == notion_page_id
        ).first()
    
    if existing_dataset:
        # Update existing
        logger.info("[INGEST][POSTGRES] Updating existing dataset: %s", existing_dataset.name)
        existing_dataset.name = name
        existing_dataset.omics_type = omics_type.value
        if description is not None:
            existing_dataset.description = description
        if file_paths is not None:
            existing_dataset.file_paths = file_paths
        if file_urls is not None:
            existing_dataset.file_urls = file_urls
        if organism is not None:
            existing_dataset.organism = organism
        if sample_type is not None:
            existing_dataset.sample_type = sample_type
        if disease is not None:
            existing_dataset.disease = disease
        if external_ids is not None:
            existing_dataset.external_ids = external_ids
        
        db.commit()
        db.refresh(existing_dataset)
        return existing_dataset
    else:
        # Create new
        logger.info("[INGEST][POSTGRES] Creating new dataset: %s", name)
        import uuid
        new_dataset = DatasetModel(
            id=uuid.uuid4(),
            name=name,
            omics_type=omics_type.value,
            description=description,
            file_paths=file_paths or [],
            file_urls=file_urls or [],
            organism=organism or [],
            sample_type=sample_type or [],
            disease=disease or [],
            notion_page_id=notion_page_id,
            external_ids=external_ids or {},
        )
        
        # Add program relationships
        if program_ids:
            from amprenta_rag.database.models import Program as ProgramModel
            for program_id in program_ids:
                program = db.query(ProgramModel).filter(ProgramModel.id == program_id).first()
                if program:
                    new_dataset.programs.append(program)
        
        # Add experiment relationships
        if experiment_ids:
            from amprenta_rag.database.models import Experiment as ExperimentModel
            for experiment_id in experiment_ids:
                experiment = db.query(ExperimentModel).filter(ExperimentModel.id == experiment_id).first()
                if experiment:
                    new_dataset.experiments.append(experiment)
        
        db.add(new_dataset)
        db.commit()
        db.refresh(new_dataset)
        return new_dataset


def embed_dataset_with_postgres_metadata(
    dataset_id: UUID,
    dataset_name: str,
    species_or_features: List[str],
    omics_type: OmicsType,
    signature_matches: Optional[List] = None,
    notion_page_id: Optional[str] = None,
) -> None:
    """
    Embed dataset into Pinecone using Postgres metadata.
    
    This uses Postgres as the source of truth for metadata while maintaining
    backward compatibility with Notion IDs.
    
    Args:
        dataset_id: Postgres UUID of the dataset
        dataset_name: Dataset name
        species_or_features: List of species/features in the dataset
        omics_type: Omics type
        signature_matches: Optional list of signature match results
        notion_page_id: Optional Notion page ID for backward compatibility
    """
    from amprenta_rag.clients.pinecone_client import get_pinecone_index
    from amprenta_rag.config import get_config
    from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
    from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
    from amprenta_rag.rag.postgres_builder import build_dataset_rag_metadata
    
    cfg = get_config()
    
    try:
        # Build text representation
        text_parts = [
            f"Dataset: {dataset_name}",
            f"Omics Type: {omics_type.value}",
            f"Data Origin: Internal – Amprenta",
            "",
            f"Contains {len(species_or_features)} features:",
        ]
        
        # Add features list (truncate if too long)
        features_list = sorted(species_or_features)
        if len(features_list) > 50:
            text_parts.append(", ".join(features_list[:50]))
            text_parts.append(f"... and {len(features_list) - 50} more features")
        else:
            text_parts.append(", ".join(features_list))
        
        # Add signature matches if available
        if signature_matches:
            text_parts.append("")
            text_parts.append("Signature Matches:")
            for match in signature_matches[:5]:  # Top 5
                match_name = getattr(match, "signature_name", "Unknown")
                match_score = getattr(match, "score", 0.0)
                match_overlap = getattr(match, "overlap_fraction", 0.0)
                text_parts.append(
                    f"- {match_name}: score {match_score:.3f}, overlap {match_overlap:.2f}"
                )
        
        dataset_text = "\n".join(text_parts)
        
        # Chunk and embed
        chunks = chunk_text(dataset_text, max_chars=2000)
        if not chunks:
            logger.warning(
                "[INGEST][POSTGRES] No chunks generated for dataset %s",
                dataset_id,
            )
            return
        
        embeddings = embed_texts(chunks)
        
        # Build metadata using Postgres builder
        base_metadata = build_dataset_rag_metadata(
            dataset_id=dataset_id,
            include_notion_id=True,  # Include Notion ID if available
        )
        
        # Upsert to Pinecone
        index = get_pinecone_index()
        vectors = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{str(dataset_id).replace('-', '')}_chunk_{order:03d}"
            
            meta = base_metadata.copy()
            meta.update({
                "snippet": chunk[:300],
            })
            
            vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(meta),
                }
            )
        
        # Batch upsert
        batch_size = 100
        for i in range(0, len(vectors), batch_size):
            batch = vectors[i : i + batch_size]
            index.upsert(vectors=batch, namespace=cfg.pinecone.namespace)
        
        logger.info(
            "[INGEST][POSTGRES] Embedded dataset %s to Pinecone (%d vectors)",
            dataset_id,
            len(vectors),
        )
        
    except Exception as e:
        logger.warning(
            "[INGEST][POSTGRES] Error embedding dataset %s: %r",
            dataset_id,
            e,
        )
        # Don't raise - embedding is non-critical

