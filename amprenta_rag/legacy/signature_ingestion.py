# amprenta_rag/ingestion/signature_ingestion.py

"""
External lipid signature ingestion into Notion databases.

This module handles ingesting external lipid signatures (from TSV, CSV, JSON, etc.)
into the Notion knowledge graph:
- Lipid Signatures database
- Lipid Signature Components database
- Lipid Species database (canonical ontology)

Creates the full relation graph: Signature → Components → Species
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from amprenta_rag.ingestion.signature_embedding import embed_signature
from amprenta_rag.ingestion.signature_linking import (
    link_component_to_feature,
    link_component_to_lipid_species,
    link_component_to_metabolite_feature,
    link_signature_to_source,
)
from amprenta_rag.ingestion.signatures import (
    find_or_create_component_page,
    find_or_create_lipid_species_page,
    find_or_create_signature_page,
    generate_signature_short_id,
    update_lipid_species_synonyms,
    update_signature_modalities,
    update_signature_page_if_needed,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature, SignatureComponent, load_signature_from_tsv

logger = get_logger(__name__)


# Backward compatibility wrappers for functions moved to other modules
# These allow existing code to continue using the old function names


def _generate_short_id(signature_name: str, version: Optional[str] = None) -> str:
    """Backward compatibility wrapper for generate_signature_short_id."""
    return generate_signature_short_id(signature_name, version)


def _find_or_create_signature_page(
    signature: Signature,
    signature_type: str = "Literature-derived",
    data_ownership: str = "Public",
    version: Optional[str] = None,
    description: Optional[str] = None,
) -> Optional[str]:
    """Backward compatibility wrapper for find_or_create_signature_page."""
    return find_or_create_signature_page(signature, signature_type, data_ownership, version, description)


def _update_signature_page_if_needed(
    page_id: str,
    signature: Signature,
    signature_type: str,
    data_ownership: str,
    version: Optional[str],
    description: Optional[str],
    short_id: str,
) -> None:
    """Backward compatibility wrapper for update_signature_page_if_needed."""
    return update_signature_page_if_needed(
        page_id, signature, signature_type, data_ownership, version, description, short_id
    )


def _find_or_create_component_page(
    component: SignatureComponent,
    signature_page_id: str,
    disease_context: Optional[List[str]] = None,
    matrix: Optional[List[str]] = None,
) -> Optional[str]:
    """Backward compatibility wrapper for find_or_create_component_page."""
    return find_or_create_component_page(component, signature_page_id, disease_context, matrix)


def _find_or_create_lipid_species_page(
    lipid_name: str,
) -> Optional[str]:
    """Backward compatibility wrapper for find_or_create_lipid_species_page."""
    return find_or_create_lipid_species_page(lipid_name)


def _update_lipid_species_synonyms(page_id: str, new_name: str) -> None:
    """Backward compatibility wrapper for update_lipid_species_synonyms."""
    return update_lipid_species_synonyms(page_id, new_name)


def _link_component_to_lipid_species(
    component_page_id: str,
    lipid_species_page_id: str,
) -> None:
    """Backward compatibility wrapper for link_component_to_lipid_species."""
    return link_component_to_lipid_species(component_page_id, lipid_species_page_id)


def link_signature_to_source(
    signature_page_id: str,
    source_page_id: str,
    source_type: str,
) -> None:
    """Backward compatibility wrapper for link_signature_to_source."""
    from amprenta_rag.ingestion.signature_linking import link_signature_to_source as link_func

    return link_func(signature_page_id, source_page_id, source_type)


def link_component_to_metabolite_feature(
    component_page_id: str,
    lipid_species_page_id: str,
) -> None:
    """Backward compatibility wrapper for link_component_to_metabolite_feature."""
    return link_component_to_metabolite_feature(component_page_id, lipid_species_page_id)


def ingest_signature_from_file(
    signature_path: Path,
    signature_type: str = "Literature-derived",
    data_ownership: str = "Public",
    version: Optional[str] = None,
    description: Optional[str] = None,
    disease_context: Optional[List[str]] = None,
    matrix: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Ingest a single signature from a file (TSV, CSV) into Notion.

    Creates the full relation graph:
    - Lipid Signature page
    - Lipid Signature Component pages (one per component)
    - Lipid Species pages (canonical ontology)

    Links: Signature → Components → Species

    Args:
        signature_path: Path to signature file (TSV or CSV)
        signature_type: "Consortium", "Literature-derived", "Open Dataset", or "Other"
        data_ownership: "Public" or appropriate value
        version: Optional version string
        description: Optional description text
        disease_context: Optional list of disease context strings for all components
        matrix: Optional list of matrix strings for all components

    Returns:
        Dictionary with ingestion results:
        - signature_page_id: Notion page ID of signature
        - component_count: Number of components created
        - species_count: Number of lipid species created/linked
        - warnings: List of warning messages
    """
    logger.info(
        "[INGEST][SIGNATURES] Loading signature from file: %s",
        signature_path,
    )

    # Load signature from file
    try:
        signature = load_signature_from_tsv(signature_path)
    except Exception as e:
        logger.error(
            "[INGEST][SIGNATURES] Error loading signature from %s: %r",
            signature_path,
            e,
        )
        raise

    logger.info(
        "[INGEST][SIGNATURES] Loaded signature '%s' with %d components",
        signature.name,
        len(signature.components),
    )

    warnings: List[str] = []

    # Create or find signature page
    signature_page_id = _find_or_create_signature_page(
        signature=signature,
        signature_type=signature_type,
        data_ownership=data_ownership,
        version=version,
        description=description,
    )

    if not signature_page_id:
        error_msg = f"Failed to create/find signature page for {signature.name}"
        logger.error(f"[INGEST][SIGNATURES] {error_msg}")
        warnings.append(error_msg)
        return {
            "signature_page_id": None,
            "component_count": 0,
            "species_count": 0,
            "warnings": warnings,
        }

    # Create component pages and link to feature pages (multi-omics support)
    component_count = 0
    features_created: Set[tuple] = set()  # (feature_type, feature_name) tuples

    # Import normalization functions
    from amprenta_rag.ingestion.lipidomics import normalize_lipid_species
    from amprenta_rag.ingestion.metabolomics import normalize_metabolite_name
    from amprenta_rag.ingestion.proteomics import normalize_protein_identifier
    from amprenta_rag.ingestion.transcriptomics import normalize_gene_identifier

    for component in signature.components:
        try:
            # Get feature type and name
            feature_type = getattr(component, "feature_type", "lipid")  # Default to lipid for backward compat
            feature_name_raw = getattr(component, "feature_name", component.species)

            # Normalize feature name based on type
            feature_name_normalized = feature_name_raw
            if feature_type == "gene":
                feature_name_normalized = normalize_gene_identifier(feature_name_raw)
            elif feature_type == "protein":
                feature_name_normalized = normalize_protein_identifier(feature_name_raw)
            elif feature_type == "metabolite":
                feature_name_normalized = normalize_metabolite_name(feature_name_raw)
            elif feature_type == "lipid":
                normalized = normalize_lipid_species(feature_name_raw)
                if normalized:
                    feature_name_normalized = normalized
                else:
                    feature_name_normalized = feature_name_raw

            # Create or find component page
            component_page_id = _find_or_create_component_page(
                component=component,
                signature_page_id=signature_page_id,
                disease_context=disease_context,
                matrix=matrix,
            )

            if not component_page_id:
                warning = f"Failed to create component page for {feature_name_raw} ({feature_type})"
                logger.warning(f"[INGEST][SIGNATURES] {warning}")
                warnings.append(warning)
                continue

            component_count += 1

            # Link component to appropriate feature page
            try:
                link_component_to_feature(
                    component_page_id=component_page_id,
                    feature_type=feature_type,
                    feature_name=feature_name_normalized,
                )

                feature_key = (feature_type, feature_name_normalized)
                if feature_key not in features_created:
                    features_created.add(feature_key)

                logger.debug(
                    "[INGEST][SIGNATURES] Linked component %s to %s feature '%s'",
                    component_page_id,
                    feature_type,
                    feature_name_normalized,
                )
            except Exception as e:
                warning = f"Error linking component to {feature_type} feature '{feature_name_normalized}': {e}"
                logger.warning(f"[INGEST][SIGNATURES] {warning}")
                warnings.append(warning)
                # Non-blocking - continue

        except Exception as e:
            feature_name = getattr(component, "feature_name", component.species)
            warning = f"Error processing component {feature_name}: {e}"
            logger.warning(f"[INGEST][SIGNATURES] {warning}")
            warnings.append(warning)

    # Embed signature into Pinecone for RAG queries
    try:
        embed_signature(
            signature_page_id=signature_page_id,
            signature=signature,
        )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error embedding signature '%s': %r",
            signature.name,
            e,
        )
        # Non-blocking - continue

    # Update signature page with Modalities field if available
    try:
        if signature.modalities:
            update_signature_modalities(signature_page_id, signature.modalities)
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error updating modalities for signature '%s': %r",
            signature.name,
            e,
        )
        # Non-blocking

    logger.info(
        "[INGEST][SIGNATURES] Ingestion complete for signature '%s': " "%d components, %d features (modalities: %s)",
        signature.name,
        component_count,
        len(features_created),
        ", ".join(signature.modalities) if signature.modalities else "none",
    )

    return {
        "signature_page_id": signature_page_id,
        "component_count": component_count,
        "species_count": len(features_created),  # Renamed for backward compatibility
        "feature_count": len(features_created),
        "modalities": signature.modalities or [],
        "warnings": warnings,
    }
