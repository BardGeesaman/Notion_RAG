"""
Postgres-only dataset ingestion module.

Handles ingestion of experimental datasets directly from Postgres into Pinecone.
No Notion dependencies - works entirely with Postgres data.

This module provides fast, scalable ingestion without slow Notion API calls.
"""

from __future__ import annotations

import json
import textwrap
from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.clients.vector_store import get_vector_store
from amprenta_rag.config import get_config
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.ingestion.feature_extraction import extract_features_from_mwtab
from amprenta_rag.ingestion.features.postgres_linking import (
    batch_link_features_to_dataset_in_postgres,
)
from amprenta_rag.ingestion.metadata.postgres_semantic_extraction import (
    enhance_metadata_with_semantic_extraction,
)
from amprenta_rag.ingestion.mwtab_extraction import (
    extract_metadata_from_mwtab,
    fetch_mwtab_from_api,
)
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.postgres_signature_matching import (
    find_matching_signatures_for_postgres_dataset,
)
from amprenta_rag.ingestion.postgres_signature_linking import (
    link_signature_to_dataset_in_postgres,
)
from amprenta_rag.ingestion.signature_matching import (
    map_raw_lipid_to_canonical_species,
)
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import FeatureType, OmicsType

logger = get_logger(__name__)


def build_dataset_text_content(dataset: DatasetModel) -> str:
    """
    Build text content for embedding from Postgres dataset fields.

    Args:
        dataset: DatasetModel instance from Postgres

    Returns:
        Combined text content for embedding
    """
    parts: List[str] = []

    # Title
    if dataset.name:
        parts.append(f"Dataset: {dataset.name}")

    # Description
    if dataset.description:
        parts.append(f"\nDescription: {dataset.description}")

    # Omics type
    if dataset.omics_type:
        parts.append(f"\nOmics Type: {dataset.omics_type}")

    # Organism/Model systems
    if dataset.organism:
        orgs = [str(o) for o in dataset.organism if isinstance(o, str) and o]
        if orgs:
            parts.append(f"\nOrganism: {', '.join(orgs)}")

    # Sample type/Matrix
    if dataset.sample_type:
        sample_types = [str(s) for s in dataset.sample_type if isinstance(s, str) and s]
        if sample_types:
            parts.append(f"\nSample Type: {', '.join(sample_types)}")

    # Disease
    if dataset.disease:
        diseases = [str(d) for d in dataset.disease if isinstance(d, str) and d]
        if diseases:
            parts.append(f"\nDisease: {', '.join(diseases)}")

    # File URLs (repository links)
    if dataset.file_urls:
        urls = [u for u in dataset.file_urls if isinstance(u, str)]
        if urls:
            parts.append(f"\nData Files: {', '.join(urls)}")

    # External IDs (repository study IDs, DOI, etc.)
    if dataset.external_ids:
        external_info = []
        for key, value in dataset.external_ids.items():
            if not value:
                continue
            if key == "doi":
                external_info.append(f"DOI: {value}")
            elif key.endswith("_study_id"):
                repo_name = key.replace("_study_id", "").upper()
                external_info.append(f"{repo_name} Study ID: {value}")
            else:
                external_info.append(f"{key}: {value}")
        if external_info:
            parts.append(f"\nExternal Identifiers: {'; '.join(external_info)}")

    # Linked programs
    if dataset.programs:
        program_names = [str(p.name) for p in dataset.programs if getattr(p, "name", None)]
        if program_names:
            parts.append(f"\nPrograms: {', '.join(program_names)}")

    # Linked experiments
    if dataset.experiments:
        experiment_names = [str(e.name) for e in dataset.experiments if getattr(e, "name", None)]
        if experiment_names:
            parts.append(f"\nExperiments: {', '.join(experiment_names)}")

    return "\n".join(parts)


def get_dataset_metadata_from_postgres(dataset: DatasetModel) -> Dict[str, Any]:
    """
    Extract metadata from Postgres dataset for Pinecone embedding.

    Args:
        dataset: DatasetModel instance

    Returns:
        Metadata dictionary for Pinecone
    """
    metadata: Dict[str, Any] = {
        "source": "Dataset",
        "source_type": "Dataset",
        "dataset_id": str(dataset.id),
        "title": dataset.name or "Untitled Dataset",
        "omics_type": dataset.omics_type,
    }

    # Add array fields
    if dataset.disease:
        metadata["diseases"] = dataset.disease
    if dataset.sample_type:
        metadata["matrix"] = dataset.sample_type
    if dataset.organism:
        metadata["model_systems"] = dataset.organism

    # Add scientific metadata fields
    if dataset.methods:
        metadata["methods"] = dataset.methods
    if dataset.summary:
        metadata["summary"] = dataset.summary
    if dataset.results:
        metadata["results"] = dataset.results
    if dataset.conclusions:
        metadata["conclusions"] = dataset.conclusions
    if dataset.dataset_source_type:
        metadata["dataset_source_type"] = dataset.dataset_source_type
    if dataset.data_origin:
        metadata["data_origin"] = dataset.data_origin

    # Add external IDs for filtering
    if dataset.external_ids:
        for key, value in dataset.external_ids.items():
            if key == "doi" and value:
                metadata["doi"] = value
            elif key.endswith("_study_id") and value:
                metadata[f"study_id_{key}"] = value
            elif key == "source_url" and value:
                metadata["source_url"] = value

    return metadata


def ingest_dataset_from_postgres(
    dataset_id: UUID,
    force: bool = False,
    update_notion: bool = False,
) -> None:
    """
    Ingest a dataset directly from Postgres into Pinecone.

    This function performs the complete dataset ingestion pipeline using only
    Postgres data - no Notion API calls required:

    1. Fetches dataset from Postgres
    2. Builds text content from Postgres fields
    3. Fetches mwTab data from repository API (if available)
    4. Extracts metadata from Postgres fields
    5. Chunks and embeds the text
    6. Upserts vectors to Pinecone
    7. Links features to Postgres
    8. Matches signatures (if enabled)

    Args:
        dataset_id: Postgres dataset UUID
        force: If True, re-ingest even if already embedded
        update_notion: If True, optionally update Notion (disabled by default)

    Raises:
        ValueError: If dataset not found
        Exception: If ingestion fails at any step
    """
    logger.info(
        "[INGEST][POSTGRES] Starting Postgres-only ingestion for dataset %s",
        dataset_id,
    )

    # Legacy arg kept for backwards compatibility; Notion sync is generally disabled in Postgres-first mode.
    if update_notion:
        logger.debug("[INGEST][POSTGRES] update_notion requested (currently ignored)")

    with db_session() as db:
        dataset = db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()

        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found in Postgres")

        logger.info(
            "[INGEST][POSTGRES] Found dataset: %s (%s)",
            dataset.name,
            dataset.omics_type,
        )

        full_text = build_dataset_text_content(dataset)

        if not full_text or len(full_text.strip()) < 50:
            logger.warning(
                "[INGEST][POSTGRES] Dataset %s has very little text content; skipping ingestion",
                dataset_id,
            )
            return

        base_meta = get_dataset_metadata_from_postgres(dataset)

        try:
            base_meta = enhance_metadata_with_semantic_extraction(base_meta, full_text)
        except Exception as e:
            logger.debug(
                "[INGEST][POSTGRES] Error enhancing metadata with semantic extraction: %r",
                e,
            )

        mwtab_data: Optional[Dict[str, Any]] = None
        study_id: Optional[str] = None

        if dataset.external_ids:
            for key, value in dataset.external_ids.items():
                if key == "mw_study_id" or (key.endswith("_study_id") and "mw" in key.lower()):
                    study_id = value
                    break

            if study_id:
                try:
                    logger.info(
                        "[INGEST][POSTGRES] Fetching mwTab data for study %s",
                        study_id,
                    )
                    mwtab_data = fetch_mwtab_from_api(study_id)

                    if mwtab_data:
                        mwtab_text = json.dumps(mwtab_data, indent=2)
                        full_text = f"{full_text}\n\nmwTab Data:\n{mwtab_text}"
                        logger.info(
                            "[INGEST][POSTGRES] Successfully fetched mwTab data for study %s",
                            study_id,
                        )

                        try:
                            scientific_metadata = extract_metadata_from_mwtab(mwtab_data)

                            if scientific_metadata.get("methods"):
                                dataset.methods = scientific_metadata["methods"]
                            if scientific_metadata.get("summary"):
                                dataset.summary = scientific_metadata["summary"]
                            if scientific_metadata.get("results"):
                                dataset.results = scientific_metadata["results"]
                            if scientific_metadata.get("conclusions"):
                                dataset.conclusions = scientific_metadata["conclusions"]
                            if scientific_metadata.get("dataset_source_type"):
                                dataset.dataset_source_type = scientific_metadata["dataset_source_type"]
                            if scientific_metadata.get("data_origin"):
                                dataset.data_origin = scientific_metadata["data_origin"]

                            if scientific_metadata.get("model_systems"):
                                existing_organism = set(dataset.organism or [])
                                existing_organism.update(scientific_metadata["model_systems"])
                                dataset.organism = sorted(list(existing_organism))

                            if scientific_metadata.get("disease_terms"):
                                existing_disease = set(dataset.disease or [])
                                existing_disease.update(scientific_metadata["disease_terms"])
                                dataset.disease = sorted(list(existing_disease))

                            if scientific_metadata.get("matrix_terms"):
                                existing_sample_type = set(dataset.sample_type or [])
                                existing_sample_type.update(scientific_metadata["matrix_terms"])
                                dataset.sample_type = sorted(list(existing_sample_type))

                            if scientific_metadata.get("source_url"):
                                if not dataset.external_ids:
                                    dataset.external_ids = {}
                                dataset.external_ids["source_url"] = scientific_metadata["source_url"]

                            db.commit()
                            logger.info(
                                "[INGEST][POSTGRES] Updated dataset %s with scientific metadata from mwTab",
                                dataset_id,
                            )
                        except Exception as e:
                            logger.warning(
                                "[INGEST][POSTGRES] Error extracting/storing scientific metadata for dataset %s: %r",
                                dataset_id,
                                e,
                            )
                except Exception as e:
                    logger.warning(
                        "[INGEST][POSTGRES] Could not fetch mwTab data for study %s: %r",
                        study_id,
                        e,
                    )

    logger.info(
        "[INGEST][POSTGRES] Found dataset: %s (%s)",
        dataset.name,
        dataset.omics_type,
    )

    # Build text content from Postgres fields
    full_text = build_dataset_text_content(dataset)

    if not full_text or len(full_text.strip()) < 50:
        logger.warning(
            "[INGEST][POSTGRES] Dataset %s has very little text content; skipping ingestion",
            dataset_id,
        )
        return

    # Extract metadata from Postgres
    base_meta = get_dataset_metadata_from_postgres(dataset)

    # Enhance metadata with semantic extraction from text content
    # This extracts diseases, targets, signatures from text without Notion dependency
    try:
        base_meta = enhance_metadata_with_semantic_extraction(base_meta, full_text)
    except Exception as e:
        logger.debug(
            "[INGEST][POSTGRES] Error enhancing metadata with semantic extraction: %r",
            e,
        )
        # Non-blocking - continue with base metadata

    # Try to fetch mwTab data if we have a Metabolomics Workbench study ID
    # reuse mwtab_data/study_id if needed; no redeclarations

    if dataset.external_ids:
        # Check for MW study ID
        for key, value in dataset.external_ids.items():
            if key == "mw_study_id" or (key.endswith("_study_id") and "mw" in key.lower()):
                study_id = value
                break

        # Fetch mwTab from API
        if study_id:
            try:
                logger.info(
                    "[INGEST][POSTGRES] Fetching mwTab data for study %s",
                    study_id,
                )
                mwtab_data = fetch_mwtab_from_api(study_id)

                if mwtab_data:
                    # Add mwTab JSON to text content
                    mwtab_text = json.dumps(mwtab_data, indent=2)
                    full_text = f"{full_text}\n\nmwTab Data:\n{mwtab_text}"
                    logger.info(
                        "[INGEST][POSTGRES] Successfully fetched mwTab data for study %s",
                        study_id,
                    )

                    # Extract and store scientific metadata from mwTab
                    try:
                        scientific_metadata = extract_metadata_from_mwtab(mwtab_data)

                        # Update dataset with extracted metadata
                        if scientific_metadata.get("methods"):
                            dataset.methods = scientific_metadata["methods"]
                        if scientific_metadata.get("summary"):
                            dataset.summary = scientific_metadata["summary"]
                        if scientific_metadata.get("results"):
                            dataset.results = scientific_metadata["results"]
                        if scientific_metadata.get("conclusions"):
                            dataset.conclusions = scientific_metadata["conclusions"]
                        if scientific_metadata.get("dataset_source_type"):
                            dataset.dataset_source_type = scientific_metadata["dataset_source_type"]
                        if scientific_metadata.get("data_origin"):
                            dataset.data_origin = scientific_metadata["data_origin"]

                        # Merge model systems, disease terms, matrix terms from mwTab with existing arrays
                        if scientific_metadata.get("model_systems"):
                            existing_organism = set(dataset.organism or [])
                            existing_organism.update(scientific_metadata["model_systems"])
                            dataset.organism = sorted(list(existing_organism))

                        if scientific_metadata.get("disease_terms"):
                            existing_disease = set(dataset.disease or [])
                            existing_disease.update(scientific_metadata["disease_terms"])
                            dataset.disease = sorted(list(existing_disease))

                        if scientific_metadata.get("matrix_terms"):
                            existing_sample_type = set(dataset.sample_type or [])
                            existing_sample_type.update(scientific_metadata["matrix_terms"])
                            dataset.sample_type = sorted(list(existing_sample_type))

                        # Update source URL if available
                        if scientific_metadata.get("source_url"):
                            if not dataset.external_ids:
                                dataset.external_ids = {}
                            dataset.external_ids["source_url"] = scientific_metadata["source_url"]

                        # Commit metadata updates
                        db.commit()
                        logger.info(
                            "[INGEST][POSTGRES] Updated dataset %s with scientific metadata from mwTab",
                            dataset_id,
                        )
                    except Exception as e:
                        logger.warning(
                            "[INGEST][POSTGRES] Error extracting/storing scientific metadata for dataset %s: %r",
                            dataset_id,
                            e,
                        )
            except Exception as e:
                logger.warning(
                    "[INGEST][POSTGRES] Could not fetch mwTab data for study %s: %r",
                    study_id,
                    e,
                )

    # Chunk and embed
    chunks = chunk_text(full_text)
    if not chunks:
        logger.warning(
            "[INGEST][POSTGRES] No chunks produced for dataset %s; skipping",
            dataset_id,
        )
        return

    logger.info(
        "[INGEST][POSTGRES] Generated %d chunk(s) for dataset %s",
        len(chunks),
        dataset_id,
    )

    try:
        embeddings = embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][POSTGRES] Error embedding chunks for dataset %s: %r",
            dataset_id,
            e,
        )
        raise

    # Prepare vectors for Pinecone
    cfg = get_config()
    store = get_vector_store()

    vectors: List[Dict[str, Any]] = []
    embedding_ids: List[str] = []

    dataset_id_str = str(dataset_id).replace("-", "")

    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{dataset_id_str}_chunk_{order:03d}"
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
            "[INGEST][POSTGRES] No vectors to upsert for dataset %s; skipping",
            dataset_id,
        )
        return

    logger.info(
        "[INGEST][POSTGRES] Upserting %d vectors into Pinecone for dataset %s",
        len(vectors),
        dataset_id,
    )

    # Batch upserts to avoid Pinecone size limits
    batch_size = 100
    try:
        for i in range(0, len(vectors), batch_size):
            batch = vectors[i : i + batch_size]
            batch_num = (i // batch_size) + 1
            total_batches = (len(vectors) + batch_size - 1) // batch_size

            logger.debug(
                "[INGEST][POSTGRES] Upserting batch %d/%d (%d vectors) for dataset %s",
                batch_num,
                total_batches,
                len(batch),
                dataset_id,
            )

            store.upsert(vectors=batch, namespace=cfg.pinecone.namespace)

            if batch_num % 10 == 0 or batch_num == total_batches:
                logger.info(
                    "[INGEST][POSTGRES] Completed batch %d/%d for dataset %s",
                    batch_num,
                    total_batches,
                    dataset_id,
                )
    except Exception as e:
        logger.error(
            "[INGEST][POSTGRES] Error upserting vectors for dataset %s: %r",
            dataset_id,
            e,
        )
        raise

    # Extract and link features from mwTab (if available)
    if mwtab_data:
        try:
            feature_names = extract_features_from_mwtab(mwtab_data)
            if feature_names:
                feature_type_map: Dict[OmicsType, str] = {
                    OmicsType.METABOLOMICS: "metabolite",
                    OmicsType.LIPIDOMICS: "lipid",
                    OmicsType.PROTEOMICS: "protein",
                    OmicsType.TRANSCRIPTOMICS: "gene",
                }
                feature_type = feature_type_map.get(
                    OmicsType(dataset.omics_type), "metabolite"
                )

                # Link features to Postgres dataset
                features_to_link = [(name, feature_type) for name in feature_names]
                batch_link_features_to_dataset_in_postgres(
                    features=features_to_link,
                    dataset_id=dataset_id,
                    db=db,
                )
                logger.info(
                    "[INGEST][POSTGRES] Linked %d features to dataset %s",
                    len(feature_names),
                    dataset_id,
                )
        except Exception as e:
            logger.warning(
                "[INGEST][POSTGRES] Error linking features for dataset %s: %r",
                dataset_id,
                e,
            )

    # Match dataset against signatures (if enabled) - Postgres-first approach
    cfg = get_config()
    if cfg.pipeline.enable_signature_scoring:
        try:
            logger.info(
                "[INGEST][POSTGRES] Matching dataset %s against signatures",
                dataset_id,
            )

            # Extract legacy species from mwTab for backward compatibility (if mwTab available)
            dataset_species_set: set[str] = set()
            if mwtab_data:
                metabolite_sections = [
                    "MS_METABOLITE_DATA",
                    "GC_METABOLITE_DATA",
                    "LC_METABOLITE_DATA",
                    "METABOLITE_DATA",
                ]
                for section_key in metabolite_sections:
                    if section_key in mwtab_data:
                        section = mwtab_data.get(section_key, {})
                        data_array = section.get("Data", [])
                        if isinstance(data_array, list):
                            for row in data_array:
                                if isinstance(row, dict):
                                    for key in row.keys():
                                        if key.lower() in [
                                            "metabolite",
                                            "metabolite_name",
                                            "compound",
                                            "name",
                                        ]:
                                            raw_name = row.get(key)
                                            if raw_name and isinstance(raw_name, str):
                                                raw_name = raw_name.strip()
                                                if raw_name:
                                                    canonical = (
                                                        map_raw_lipid_to_canonical_species(
                                                            raw_name
                                                        )
                                                    )
                                                    dataset_species_set.add(
                                                        canonical if canonical else raw_name
                                                    )

            # Use Postgres-based signature matching
            matches = find_matching_signatures_for_postgres_dataset(
                dataset_id=dataset_id,
                dataset_species=dataset_species_set if dataset_species_set else None,
                overlap_threshold=cfg.pipeline.signature_overlap_threshold,
                omics_type=dataset.omics_type,
            )

            if matches:
                logger.info(
                    "[INGEST][POSTGRES] Found %d matching signature(s) for dataset %s",
                    len(matches),
                    dataset_id,
                )

                # Link matched signatures to dataset in Postgres
                for match in matches:
                    try:
                        # Extract signature_id from match
                        # In Postgres mode, signature_page_id contains the UUID string
                        from uuid import UUID as UUIDType

                        try:
                            signature_id = UUIDType(match.signature_page_id)
                        except (ValueError, AttributeError, TypeError):
                            logger.warning(
                                "[INGEST][POSTGRES] Could not parse signature ID from match: %s",
                                match.signature_page_id,
                            )
                            continue

                        # Link signature to dataset
                        link_signature_to_dataset_in_postgres(
                            signature_id=signature_id,
                            dataset_id=dataset_id,
                            match_score=match.overlap_fraction,
                            db=db,
                        )
                    except Exception as e:
                        logger.warning(
                            "[INGEST][POSTGRES] Error linking signature match to dataset: %r",
                            e,
                        )
            else:
                logger.info(
                    "[INGEST][POSTGRES] No signature matches found for dataset %s",
                    dataset_id,
                )
        except Exception as e:
            logger.warning(
                "[INGEST][POSTGRES] Error matching signatures for dataset %s: %r",
                dataset_id,
                e,
            )

    # Detect and ingest signatures from content (Postgres-based)
    # Signature detection is enabled by default, can be disabled via config if needed
    try:
        try:
            source_metadata = {
                "diseases": base_meta.get("diseases", []),
                "matrix": base_meta.get("matrix", []),
                "model_systems": base_meta.get("model_systems", []),
            }

            # Use Postgres-based signature detection
            from amprenta_rag.ingestion.postgres_signature_detection import (
                detect_and_ingest_signatures_from_postgres_content,
            )

            detect_and_ingest_signatures_from_postgres_content(
                all_text_content=full_text,
                attachment_paths=[],
                source_type="dataset",
                source_metadata=source_metadata,
                source_name=dataset.name or "Untitled Dataset",
                source_dataset_id=dataset_id,
            )
        except Exception as e:
            logger.debug(
                "[INGEST][POSTGRES] Signature detection skipped or failed (optional feature): %r",
                e,
            )
    except Exception as e:
        logger.debug(
            "[INGEST][POSTGRES] Signature detection skipped or failed (optional feature): %r",
            e,
        )

    logger.info(
        "[INGEST][POSTGRES] Postgres-only ingestion complete for dataset %s",
        dataset_id,
    )

