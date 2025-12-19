"""
Omics Service - Main Orchestrator for Multi-Omics Dataset Ingestion

This module provides the primary entry point for ingesting omics datasets into the
Amprenta platform. It handles routing to appropriate omics-specific parsers and
coordinates the complete Postgres-first ingestion pipeline.

**Architecture**: Postgres-First (100% Complete)
- All data stored directly in PostgreSQL (sole system of record)
- No Notion dependency for core functionality
- Notion sync optional via ENABLE_NOTION_SYNC flag

**Supported Omics Types**:
- Lipidomics
- Metabolomics
- Proteomics
- Transcriptomics

**Complete Ingestion Flow**:
1. Parse file → Extract features → Normalize names
2. Create/update Dataset record in Postgres
3. Link features (find/create Feature records)
4. Match signatures (score against all signatures)
5. Build text content from Postgres fields
6. Create RAGChunk records
7. Embed and upsert to Pinecone

See: docs/INGESTION_ARCHITECTURE.md for complete data flow documentation
"""

from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.domain.omics import OmicsDatasetIngestRequest
from amprenta_rag.config import AUTO_LINK_ENABLED, AUTO_LINK_MIN_CONFIDENCE
from amprenta_rag.ingestion.auto_linking import (
    infer_experiment_from_metadata,
    infer_program_from_metadata,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def ingest_dataset_from_file(req: OmicsDatasetIngestRequest) -> UUID:
    """
    Main entry point for omics dataset ingestion.

    This function orchestrates the complete ingestion pipeline:
    1. Routes to appropriate omics-specific parser based on req.omics_type
    2. Creates/updates Dataset record in Postgres
    3. Delegates to postgres_dataset_ingestion for complete pipeline:
       - Feature linking
       - Signature matching
       - RAGChunk creation
       - Pinecone embedding
    4. Returns the dataset UUID

    Args:
        req: OmicsDatasetIngestRequest containing:
            - omics_type: "lipidomics" | "metabolomics" | "proteomics" | "transcriptomics"
            - name: Dataset name
            - file_path: Path to CSV/TSV file
            - description: Optional description
            - disease: Optional disease list
            - sample_type: Optional sample type list
            - organism: Optional organism list
            - And other metadata fields

    Returns:
        UUID: The Postgres UUID of the created/updated dataset

    Raises:
        ValueError: If omics_type is not recognized
        Exception: If parsing or ingestion fails (ingestion_status set to "failed")

    Example:
        >>> req = OmicsDatasetIngestRequest(
        ...     omics_type="lipidomics",
        ...     name="ALS CSF Lipidomics",
        ...     file_path="data.csv",
        ...     disease=["ALS"]
        ... )
        >>> dataset_uuid = ingest_dataset_from_file(req)
        >>> print(f"Ingested dataset: {dataset_uuid}")

    Note:
        - Uses Postgres as sole system of record
        - Notion sync optional (via ENABLE_NOTION_SYNC config)
        - Ingestion status tracked in Dataset.ingestion_status field
    """
    with db_session() as db:
        dataset = (
            db.query(DatasetModel).filter(DatasetModel.id == req.dataset_path).first()
            if hasattr(req, "dataset_path")
            else None
        )

        if dataset:
            dataset.ingestion_status = "in_progress"
            db.commit()

        try:
            if req.omics_type == "lipidomics":
                parsed = _import_lipidomics(req)
            elif req.omics_type == "metabolomics":
                parsed = _import_metabolomics(req)
            elif req.omics_type == "proteomics":
                parsed = _import_proteomics(req)
            elif req.omics_type == "transcriptomics":
                parsed = _import_transcriptomics(req)
            else:
                raise ValueError(f"Unknown omics_type: {req.omics_type}")

            from amprenta_rag.ingestion.postgres_dataset_ingestion import ingest_dataset_from_postgres

            dataset_uuid = ingest_dataset_from_postgres(parsed)

            if AUTO_LINK_ENABLED:
                dataset_obj = db.query(DatasetModel).get(dataset_uuid)
                if dataset_obj:
                    if not getattr(parsed, "program_ids", None) and not dataset_obj.programs:
                        prog_id, prog_conf = infer_program_from_metadata(
                            diseases=getattr(dataset_obj, "disease", []) or [],
                            keywords=[dataset_obj.name] if dataset_obj.name else [],
                            filename=dataset_obj.file_paths[0] if getattr(dataset_obj, "file_paths", None) else None,
                            min_confidence=AUTO_LINK_MIN_CONFIDENCE,
                        )
                        if prog_id:
                            program = db.query(type(dataset_obj.programs.property.mapper.class_)).get(prog_id)  # type: ignore
                            if program and program not in dataset_obj.programs:
                                dataset_obj.programs.append(program)
                                db.commit()
                                logger.info(
                                    "[INGEST][AUTO-LINK] Linked dataset %s to program %s (conf=%.2f)",
                                    dataset_uuid,
                                    prog_id,
                                    prog_conf,
                                )

                    if not getattr(parsed, "experiment_ids", None) and not dataset_obj.experiments:
                        exp_id, exp_conf = infer_experiment_from_metadata(
                            diseases=getattr(dataset_obj, "disease", []) or [],
                            matrix=getattr(dataset_obj, "matrix", []) if hasattr(dataset_obj, "matrix") else [],
                            model_systems=getattr(dataset_obj, "model_systems", []) if hasattr(dataset_obj, "model_systems") else [],
                            min_confidence=AUTO_LINK_MIN_CONFIDENCE,
                        )
                        if exp_id:
                            experiment = db.query(type(dataset_obj.experiments.property.mapper.class_)).get(exp_id)  # type: ignore
                            if experiment and experiment not in dataset_obj.experiments:
                                dataset_obj.experiments.append(experiment)
                                db.commit()
                                logger.info(
                                    "[INGEST][AUTO-LINK] Linked dataset %s to experiment %s (conf=%.2f)",
                                    dataset_uuid,
                                    exp_id,
                                    exp_conf,
                                )

            if dataset:
                dataset.ingestion_status = "complete"
                db.commit()

            return dataset_uuid

        except Exception:
            if dataset:
                dataset.ingestion_status = "failed"
                db.commit()
            raise


def reingest_dataset_from_postgres(dataset_id: UUID, force: bool = False) -> None:
    """
    Re-run ingestion pipeline for an existing dataset.

    Useful for:
    - Re-embedding after text content changes
    - Re-scoring against new signatures
    - Fixing ingestion errors

    Args:
        dataset_id: UUID of existing dataset in Postgres
        force: If True, re-ingest even if already embedded

    Note:
        This does NOT re-parse the original file - it works with existing
        Dataset and Feature records in Postgres.
    """
    from amprenta_rag.ingestion.postgres_dataset_ingestion import ingest_dataset_from_postgres

    ingest_dataset_from_postgres(dataset_id, force=force)


# Stubs for per-omics normalizer bindings (to be imported from each per-omics submodule)
def _import_lipidomics(req):
    from amprenta_rag.ingestion.lipidomics.ingestion import parse_and_normalize

    return parse_and_normalize(req)


def _import_metabolomics(req):
    from amprenta_rag.ingestion.metabolomics.ingestion import parse_and_normalize

    return parse_and_normalize(req)


def _import_proteomics(req):
    from amprenta_rag.ingestion.proteomics.ingestion import parse_and_normalize

    return parse_and_normalize(req)


def _import_transcriptomics(req):
    from amprenta_rag.ingestion.transcriptomics.ingestion import parse_and_normalize

    return parse_and_normalize(req)
