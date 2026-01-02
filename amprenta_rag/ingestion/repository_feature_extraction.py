"""
Unified feature extraction for all repository types.

Extracts features (genes, proteins, metabolites, lipids) from repository datasets
and links them to Postgres datasets.

Supports:
- GEO: Genes from expression matrices (using GEOparse library)
- PRIDE: Proteins from identification files
- MetaboLights: Metabolites from ISA-Tab data files
- MW: Metabolites from REST API /data endpoint
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Set, Tuple
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset
from amprenta_rag.ingestion.features.postgres_linking import (
    batch_link_features_to_dataset_in_postgres,
)
from amprenta_rag.ingestion.repositories.geo import extract_geo_features_with_geoparse
from amprenta_rag.ingestion.repositories.pride import extract_pride_proteins_from_data_files
from amprenta_rag.ingestion.repositories.metabolights import extract_metabolights_metabolites_from_isa_tab
from amprenta_rag.ingestion.repositories.mw import extract_mw_metabolites_from_data_endpoint
from amprenta_rag.ingestion.genomics.pipeline import extract_gene_counts_from_salmon
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import FeatureType

logger = get_logger(__name__)


def extract_ena_features(
    study_id: str,
    download_dir: Path = None,
) -> Set[str]:
    """
    Extract gene features from ENA study.
    
    Uses existing Salmon quantification output if available.
    
    Args:
        study_id: ENA run accession (e.g., ERR123456)
        download_dir: Directory containing quantification results
        
    Returns:
        Set of gene identifiers from quantification
    """
    if download_dir is None:
        download_dir = Path("./quants")
    
    # Look for Salmon quant.sf output
    quant_file = download_dir / study_id / "quant.sf"
    
    if not quant_file.exists():
        # Try alternate location
        quant_file = download_dir / f"{study_id}_quant" / "quant.sf"
    
    if not quant_file.exists():
        logger.warning(
            "[FEATURE-EXTRACT] No quantification file found for ENA study %s",
            study_id,
        )
        return set()
    
    # Extract gene counts using existing function
    gene_counts = extract_gene_counts_from_salmon(quant_file)
    
    # Return gene IDs with non-zero expression
    gene_set = set(gene_counts.keys())
    
    logger.info(
        "[FEATURE-EXTRACT] Extracted %d genes from ENA study %s",
        len(gene_set),
        study_id,
    )
    
    return gene_set


def extract_features_from_repository_dataset(
    dataset_id: UUID,
    repository: str,
    study_id: Optional[str] = None,
    download_dir: Path = None,
) -> int:
    """
    Extract features from a repository dataset and link them to Postgres.

    Supports:
    - GEO: Genes from expression data (using GEOparse)
    - PRIDE: Proteins from TSV/CSV protein tables
    - MetaboLights: Metabolites from ISA-Tab data files
    - MW: Metabolites from REST API /data endpoint

    Args:
        dataset_id: Postgres dataset UUID
        repository: Repository name (GEO, PRIDE, MetaboLights, MW)
        study_id: Optional repository-specific study ID
        download_dir: Directory for temporary files

    Returns:
        Number of features linked
    """
    with db_session() as db:
        try:
            dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if not dataset:
                logger.error("[FEATURE-EXTRACT] Dataset %s not found", dataset_id)
                return 0

            if not study_id:
                if hasattr(dataset, "external_ids") and dataset.external_ids:
                    if isinstance(dataset.external_ids, dict):
                        study_id = dataset.external_ids.get("study_id") or dataset.external_ids.get("accession")

            if not study_id:
                logger.warning(
                    "[FEATURE-EXTRACT] No study_id provided and couldn't extract from dataset %s",
                    dataset_id,
                )
                return 0

            logger.info(
                "[FEATURE-EXTRACT] Extracting features from %s study %s for dataset %s",
                repository,
                study_id,
                dataset_id,
            )

            features: List[Tuple[str, str]] = []

            if repository.upper() == "GEO":
                gene_set = extract_geo_features_with_geoparse(
                    study_id=study_id,
                    download_dir=download_dir,
                )
                features = [(gene, FeatureType.GENE.value) for gene in gene_set]

            elif repository.upper() == "PRIDE":
                protein_set = extract_pride_proteins_from_data_files(
                    study_id=study_id,
                    download_dir=download_dir,
                )
                features = [(protein, FeatureType.PROTEIN.value) for protein in protein_set]

            elif repository.upper() == "METABOLIGHTS":
                metabolite_set = extract_metabolights_metabolites_from_isa_tab(
                    study_id=study_id,
                    download_dir=download_dir,
                )
                features = [(metabolite, FeatureType.METABOLITE.value) for metabolite in metabolite_set]

            elif repository.upper() in ["MW", "METABOLOMICS WORKBENCH"]:
                metabolite_set = extract_mw_metabolites_from_data_endpoint(
                    study_id=study_id,
                )
                features = [(metabolite, FeatureType.METABOLITE.value) for metabolite in metabolite_set]

            elif repository.upper() == "ENA":
                gene_set = extract_ena_features(
                    study_id=study_id,
                    download_dir=download_dir,
                )
                features = [(gene, FeatureType.GENE.value) for gene in gene_set]

            else:
                logger.warning("[FEATURE-EXTRACT] Unknown repository: %s", repository)
                return 0

            if not features:
                logger.warning("[FEATURE-EXTRACT] No features extracted from %s study %s", repository, study_id)
                return 0

            logger.info(
                "[FEATURE-EXTRACT] Linking %d features to dataset %s",
                len(features),
                dataset_id,
            )

            results = batch_link_features_to_dataset_in_postgres(
                features=features,
                dataset_id=dataset_id,
                db=db,
            )

            linked_count = sum(1 for v in results.values() if v)
            logger.info(
                "[FEATURE-EXTRACT] Successfully linked %d/%d features to dataset %s",
                linked_count,
                len(features),
                dataset_id,
            )

            return linked_count

        except Exception as e:
            logger.error(
                "[FEATURE-EXTRACT] Error extracting features for dataset %s: %r",
                dataset_id,
                e,
            )
            return 0
