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

import io
import os
import re
import time
from pathlib import Path
from typing import List, Optional, Set, Tuple, Dict, Any, cast
from uuid import UUID

import pandas as pd
import requests

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset
from amprenta_rag.ingestion.features.postgres_linking import (
    batch_link_features_to_dataset_in_postgres,
)
from amprenta_rag.ingestion.repositories.geo import extract_geo_features_with_geoparse
from amprenta_rag.ingestion.repositories.pride import extract_pride_proteins_from_data_files
from amprenta_rag.ingestion.repositories.metabolights import extract_metabolights_metabolites_from_isa_tab
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import FeatureType

logger = get_logger(__name__)


def extract_mw_metabolites_from_data_endpoint(
    study_id: str,
) -> Set[str]:
    """
    Extract metabolite features from Metabolomics Workbench using the cleaner /data endpoint.

    This is much more robust than MetaboLights because:
    - Returns clean JSON directly (no file parsing)
    - No brittle API endpoints
    - More stable (NIH-hosted)

    Args:
        study_id: MW study ID (e.g., "ST000001")

    Returns:
        Set of normalized metabolite identifiers
    """
    metabolite_set: Set[str] = set()

    # Ensure study_id format
    study_id = study_id.upper()
    if not study_id.startswith("ST"):
        logger.error("[FEATURE-EXTRACT][MW] Invalid study ID format: %s", study_id)
        return metabolite_set

    logger.info("[FEATURE-EXTRACT][MW] Extracting metabolites from %s using /data endpoint", study_id)

    try:
        # Use the cleaner REST API endpoint
        base_url = "https://www.metabolomicsworkbench.org/rest"
        url = f"{base_url}/study/study_id/{study_id}/data"

        logger.info("[FEATURE-EXTRACT][MW] Fetching data from: %s", url)

        from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
        headers = {"User-Agent": REPOSITORY_USER_AGENT}
        response = requests.get(url, headers=headers, timeout=30)

        # CASE 1: Success
        if response.status_code == 200:
            data = response.json()

            # MW API returns a dictionary where keys are row numbers (e.g., '1', '2', '3')
            # Each value is a dict with metabolite information including:
            # - 'metabolite_name': The metabolite name
            # - 'refmet_name': Reference metabolite name
            # - 'metabolite_id': Metabolite ID
            # - 'DATA': Quantification values

            from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name

            if isinstance(data, dict):
                # Iterate through all rows in the data dictionary
                for row_key, metabolite_data in data.items():
                    if isinstance(metabolite_data, dict):
                        # Extract metabolite name from various possible fields
                        metabolite_id = (
                            metabolite_data.get('refmet_name') or  # Reference metabolite name (preferred)
                            metabolite_data.get('metabolite_name') or  # Metabolite name
                            metabolite_data.get('metabolite_identification') or
                            metabolite_data.get('database_identifier') or
                            metabolite_data.get('name') or
                            metabolite_data.get('chemical_name')
                        )

                        if metabolite_id:
                            # Normalize metabolite identifier
                            normalized = normalize_metabolite_name(str(metabolite_id))
                            if normalized and len(normalized) > 1:
                                metabolite_set.add(normalized)

            elif isinstance(data, list):
                # Alternative structure: list of metabolite dictionaries
                for item in data:
                    if isinstance(item, dict):
                        metabolite_id = (
                            item.get('refmet_name') or
                            item.get('metabolite_name') or
                            item.get('metabolite_identification') or
                            item.get('database_identifier') or
                            item.get('name')
                        )
                        if metabolite_id:
                            from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name
                            normalized = normalize_metabolite_name(str(metabolite_id))
                            if normalized and len(normalized) > 1:
                                metabolite_set.add(normalized)

        # CASE 2: Known error codes
        elif response.status_code in [403, 404]:
            logger.warning("[FEATURE-EXTRACT][MW] Study %s not found or is private (403/404)", study_id)

        # CASE 3: Server error
        elif response.status_code >= 500:
            logger.warning("[FEATURE-EXTRACT][MW] Server error (500) for study %s - skipping", study_id)

        else:
            logger.warning("[FEATURE-EXTRACT][MW] Unexpected status %d for study %s", response.status_code, study_id)

        logger.info(
            "[FEATURE-EXTRACT][MW] Extracted %d unique metabolites from %s",
            len(metabolite_set),
            study_id,
        )

        return metabolite_set

    except requests.exceptions.RequestException as e:
        logger.error("[FEATURE-EXTRACT][MW] Connection error for study %s: %r", study_id, e)
        return metabolite_set
    except Exception as e:
        logger.error("[FEATURE-EXTRACT][MW] Error extracting metabolites from %s: %r", study_id, e)
        return metabolite_set


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
