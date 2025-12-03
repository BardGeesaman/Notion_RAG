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
    link_component_to_lipid_species, link_component_to_metabolite_feature,
    link_signature_to_source)
from amprenta_rag.ingestion.signature_notion_crud import (
    find_or_create_component_page, find_or_create_lipid_species_page,
    find_or_create_signature_page, generate_signature_short_id,
    update_lipid_species_synonyms, update_signature_page_if_needed)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import (Signature,
                                                      SignatureComponent,
                                                      load_signature_from_tsv)

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
    return find_or_create_signature_page(
        signature, signature_type, data_ownership, version, description
    )


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
    return find_or_create_component_page(
        component, signature_page_id, disease_context, matrix
    )


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

    # Create component pages and link to lipid species
    component_count = 0
    species_created: Set[str] = set()

    for component in signature.components:
        try:
            # Create or find component page
            component_page_id = _find_or_create_component_page(
                component=component,
                signature_page_id=signature_page_id,
                disease_context=disease_context,
                matrix=matrix,
            )

            if not component_page_id:
                warning = f"Failed to create component page for {component.species}"
                logger.warning(f"[INGEST][SIGNATURES] {warning}")
                warnings.append(warning)
                continue

            component_count += 1

            # Create or find lipid species page
            lipid_species_page_id = _find_or_create_lipid_species_page(
                lipid_name=component.species,
            )

            if lipid_species_page_id:
                # Link component to lipid species
                link_component_to_lipid_species(
                    component_page_id=component_page_id,
                    lipid_species_page_id=lipid_species_page_id,
                )

                # Link component to metabolite feature (cross-link)
                try:
                    link_component_to_metabolite_feature(
                        component_page_id=component_page_id,
                        lipid_species_page_id=lipid_species_page_id,
                    )
                except Exception as e:
                    logger.warning(
                        "[INGEST][SIGNATURES] Error linking component to metabolite feature: %r",
                        e,
                    )
                    # Non-blocking - continue

                if component.species not in species_created:
                    species_created.add(component.species)
            else:
                warning = f"Failed to create/link lipid species for {component.species}"
                logger.warning(f"[INGEST][SIGNATURES] {warning}")
                warnings.append(warning)

        except Exception as e:
            warning = f"Error processing component {component.species}: {e}"
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

    logger.info(
        "[INGEST][SIGNATURES] Ingestion complete for signature '%s': "
        "%d components, %d species",
        signature.name,
        component_count,
        len(species_created),
    )

    return {
        "signature_page_id": signature_page_id,
        "component_count": component_count,
        "species_count": len(species_created),
        "warnings": warnings,
    }
