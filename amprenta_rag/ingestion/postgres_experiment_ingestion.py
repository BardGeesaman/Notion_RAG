"""
Postgres-only experiment ingestion module.

Handles ingestion of experiments directly from Postgres into Pinecone.
No Notion dependencies - works entirely with Postgres data.
"""

from __future__ import annotations

import textwrap
from typing import Any, Dict, List, cast
from uuid import UUID
from uuid import UUID

from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Experiment as ExperimentModel
from amprenta_rag.ingestion.feature_extraction import extract_features_from_text
from amprenta_rag.ingestion.features.postgres_linking import (
    batch_link_features_to_dataset_in_postgres,
)
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import FeatureType

logger = get_logger(__name__)


def build_experiment_text_content(experiment: ExperimentModel) -> str:
    """
    Build text content for embedding from Postgres experiment fields.

    Args:
        experiment: ExperimentModel instance from Postgres

    Returns:
        Combined text content for embedding
    """
    parts: List[str] = []

    # Title
    if experiment.name:
        parts.append(f"Experiment: {experiment.name}")

    # Type
    if experiment.type:
        parts.append(f"Type: {experiment.type}")

    # Description
    if experiment.description:
        parts.append(f"\nDescription: {experiment.description}")

    # Disease
    if experiment.disease:
        diseases = [d for d in experiment.disease if d]
        if diseases:
            parts.append(f"Disease: {', '.join(diseases)}")

    # Matrix
    if experiment.matrix:
        matrices = [m for m in experiment.matrix if m]
        if matrices:
            parts.append(f"Matrix: {', '.join(matrices)}")

    # Model systems
    if experiment.model_systems:
        systems = [m for m in experiment.model_systems if m]
        if systems:
            parts.append(f"Model Systems: {', '.join(systems)}")

    # Additional experiment fields
    if experiment.targets:
        targets = [t for t in experiment.targets if t]
        if targets:
            parts.append(f"Targets: {', '.join(targets)}")
    if experiment.modality:
        mods = [m for m in experiment.modality if m]
        if mods:
            parts.append(f"Modality: {', '.join(mods)}")
    if experiment.stage:
        parts.append(f"Stage: {experiment.stage}")
    if experiment.biomarker_role:
        biomarker_roles = [b for b in experiment.biomarker_role if b]
        if biomarker_roles:
            parts.append(f"Biomarker Role: {', '.join(biomarker_roles)}")
    if experiment.treatment_arms:
        arms = [a for a in experiment.treatment_arms if a]
        if arms:
            parts.append(f"Treatment Arms: {', '.join(arms)}")

    # Linked programs
    if experiment.programs:
        program_names = [str(p.name) for p in experiment.programs if getattr(p, "name", None)]
        if program_names:
            parts.append(f"Programs: {', '.join(program_names)}")

    # Linked datasets
    if experiment.datasets:
        dataset_names = [d.name for d in experiment.datasets if d.name]
        parts.append(f"Related Datasets: {', '.join(dataset_names)}")
        parts.append(f"({len(dataset_names)} dataset(s) linked)")

    return "\n".join(parts)


def get_experiment_metadata_from_postgres(experiment: ExperimentModel) -> Dict[str, Any]:
    """
    Extract metadata from Postgres experiment for Pinecone embedding.

    Args:
        experiment: ExperimentModel instance

    Returns:
        Metadata dictionary for Pinecone
    """
    metadata: Dict[str, Any] = {
        "source": "Experiment",
        "source_type": "Experiment",
        "experiment_id": str(experiment.id),
        "title": experiment.name or "Untitled Experiment",
    }

    # Add array fields
    if experiment.disease:
        metadata["diseases"] = [d for d in experiment.disease if d]
    if experiment.matrix:
        metadata["matrix"] = [m for m in experiment.matrix if m]
    if experiment.model_systems:
        metadata["model_systems"] = [m for m in experiment.model_systems if m]

    # Add type
    if experiment.type:
        metadata["experiment_type"] = experiment.type

    # Add new experiment fields
    if experiment.targets:
        metadata["targets"] = experiment.targets
    if experiment.modality:
        metadata["modality"] = experiment.modality
    if experiment.stage:
        metadata["stage"] = experiment.stage
    if experiment.biomarker_role:
        metadata["biomarker_role"] = experiment.biomarker_role
    if experiment.treatment_arms:
        metadata["treatment_arms"] = experiment.treatment_arms

    return metadata


def ingest_experiment_from_postgres(
    experiment_id: UUID,
    force: bool = False,
    update_notion: bool = False,
) -> None:
    """
    Ingest an experiment directly from Postgres into Pinecone.

    This function performs the complete experiment ingestion pipeline using only
    Postgres data - no Notion API calls required:

    1. Fetches experiment from Postgres
    2. Builds text content from Postgres fields
    3. Extracts metadata from Postgres fields
    4. Chunks and embeds the text
    5. Upserts vectors to Pinecone
    6. Links features to Postgres
    7. Detects signatures (if enabled)

    Args:
        experiment_id: Postgres experiment UUID
        force: If True, re-ingest even if already embedded
        update_notion: If True, optionally update Notion (disabled by default)

    Raises:
        ValueError: If experiment not found
        Exception: If ingestion fails at any step
    """
    logger.info(
        "[INGEST][POSTGRES-EXPERIMENT] Starting Postgres-only ingestion for experiment %s",
        experiment_id,
    )

    with db_session() as db:
        experiment = (
            db.query(ExperimentModel)
            .filter(ExperimentModel.id == experiment_id)
            .first()
        )

        if not experiment:
            raise ValueError(f"Experiment {experiment_id} not found in Postgres")

    logger.info(
        "[INGEST][POSTGRES-EXPERIMENT] Found experiment: %s (%s)",
        experiment.name,
        experiment.type,
    )

    # Build text content from Postgres fields
    full_text = build_experiment_text_content(experiment)

    if not full_text or len(full_text.strip()) < 50:
        logger.warning(
            "[INGEST][POSTGRES-EXPERIMENT] Experiment %s has very little text content; skipping ingestion",
            experiment_id,
        )
        return

    # Extract metadata from Postgres
    base_meta = get_experiment_metadata_from_postgres(experiment)

    # Chunk and embed
    chunks = chunk_text(full_text)
    if not chunks:
        logger.warning(
            "[INGEST][POSTGRES-EXPERIMENT] No chunks produced for experiment %s; skipping",
            experiment_id,
        )
        return

    logger.info(
        "[INGEST][POSTGRES-EXPERIMENT] Generated %d chunk(s) for experiment %s",
        len(chunks),
        experiment_id,
    )

    try:
        embeddings = embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][POSTGRES-EXPERIMENT] Error embedding chunks for experiment %s: %r",
            experiment_id,
            e,
        )
        raise

    # Prepare vectors for Pinecone
    index = get_pinecone_index()
    cfg = get_config()

    vectors: List[Dict[str, Any]] = []
    embedding_ids: List[str] = []

    experiment_id_str = str(experiment_id).replace("-", "")

    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{experiment_id_str}_chunk_{order:03d}"
        embedding_ids.append(chunk_id)

        snippet = textwrap.shorten(chunk, width=300)

        meta: Dict[str, Any] = {
            **base_meta,
            "snippet": snippet,
        }

        vectors.append(
            {
                "id": chunk_id,
                "values": emb,
                "metadata": sanitize_metadata(meta),
            }
        )

    if not vectors:
        logger.warning(
            "[INGEST][POSTGRES-EXPERIMENT] No vectors to upsert for experiment %s; skipping",
            experiment_id,
        )
        return

    logger.info(
        "[INGEST][POSTGRES-EXPERIMENT] Upserting %d vectors into Pinecone for experiment %s",
        len(vectors),
        experiment_id,
    )

    # Batch upserts to avoid Pinecone size limits
    batch_size = 100
    try:
        for i in range(0, len(vectors), batch_size):
            batch = vectors[i : i + batch_size]
            batch_num = (i // batch_size) + 1
            total_batches = (len(vectors) + batch_size - 1) // batch_size

            logger.debug(
                "[INGEST][POSTGRES-EXPERIMENT] Upserting batch %d/%d (%d vectors) for experiment %s",
                batch_num,
                total_batches,
                len(batch),
                experiment_id,
            )

            index.upsert(vectors=batch, namespace=cfg.pinecone.namespace)

            if batch_num % 10 == 0 or batch_num == total_batches:
                logger.info(
                    "[INGEST][POSTGRES-EXPERIMENT] Completed batch %d/%d for experiment %s",
                    batch_num,
                    total_batches,
                    experiment_id,
                )
    except Exception as e:
        logger.error(
            "[INGEST][POSTGRES-EXPERIMENT] Error upserting vectors for experiment %s: %r",
            experiment_id,
            e,
        )
        raise

    # Extract and link features from text
    try:
        feature_names = extract_features_from_text(full_text)
        if feature_names:
            # Link features to experiment's datasets (if any)
            # Note: Features are linked to datasets, not directly to experiments
            # But we can extract them from experiment text for search
            logger.info(
                "[INGEST][POSTGRES-EXPERIMENT] Extracted %d feature(s) from experiment text",
                len(feature_names),
            )

            # If experiment has linked datasets, we could link features to those datasets
            # For now, just log that features were found
            for dataset in experiment.datasets:
                # Link features to dataset if they exist
                features_to_link = [(str(name), str(FeatureType.METABOLITE)) for name in feature_names]
                try:
                    batch_link_features_to_dataset_in_postgres(
                        features=features_to_link,
                        dataset_id=cast(UUID, dataset.id) if dataset.id is not None else UUID(int=0),  # type: ignore[arg-type]
                        db=db,
                    )
                    logger.info(
                        "[INGEST][POSTGRES-EXPERIMENT] Linked %d features to dataset %s",
                        len(feature_names),
                        dataset.id,
                    )
                except Exception as e:
                    logger.debug(
                        "[INGEST][POSTGRES-EXPERIMENT] Could not link features to dataset %s: %r",
                        dataset.id,
                        e,
                    )
    except Exception as e:
        logger.warning(
            "[INGEST][POSTGRES-EXPERIMENT] Error extracting/linking features for experiment %s: %r",
            experiment_id,
            e,
        )

    # Detect and ingest signatures from content
    try:
        source_metadata = {
            "diseases": base_meta.get("diseases", []),
            "matrix": base_meta.get("matrix", []),
            "model_systems": base_meta.get("model_systems", []),
        }

        # Signature detection (Postgres-based)
        try:
            from amprenta_rag.ingestion.postgres_signature_detection import (
                detect_and_ingest_signatures_from_postgres_content,
            )

            detect_and_ingest_signatures_from_postgres_content(
                all_text_content=full_text,
                attachment_paths=[],
                source_type="experiment",
                source_metadata=source_metadata,
                source_name=experiment.name or "Untitled Experiment",
                source_experiment_id=experiment_id,
            )
        except Exception as e:
            logger.debug(
                "[INGEST][POSTGRES-EXPERIMENT] Signature detection skipped or failed (optional feature): %r",
                e,
            )
    except Exception as e:
        logger.warning(
            "[INGEST][POSTGRES-EXPERIMENT] Error detecting/ingesting signatures for experiment %s: %r",
            experiment_id,
            e,
        )

    logger.info(
        "[INGEST][POSTGRES-EXPERIMENT] Postgres-only ingestion complete for experiment %s",
        experiment_id,
    )

