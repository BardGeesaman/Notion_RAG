"""
Internal lipidomics dataset ingestion pipeline.

This module handles ingestion of internal Amprenta lipidomics files (CSV/TSV),
including species normalization, signature scoring, Notion integration, and RAG embedding.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd
import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.dataset_notion_utils import (
    fetch_dataset_page, update_dataset_embedding_metadata,
    update_dataset_scientific_metadata)
from amprenta_rag.ingestion.omics_ingestion_utils import (
    attach_file_to_page, create_omics_dataset_page,
    link_to_programs_and_experiments)
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.signature_matching import (
    find_matching_signatures_for_dataset, map_raw_lipid_to_canonical_species,
    update_dataset_with_signature_matches)
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

__all__ = [
    "ingest_lipidomics_file",
    "normalize_lipid_species",
    "extract_species_from_file",
]


def normalize_lipid_species(raw_name: str) -> Optional[str]:
    """
    Normalize raw lipid name to canonical species format.

    Extended version of map_raw_lipid_to_canonical_species with support for
    more vendor formats and edge cases.

    Examples:
        "CER 16:0" → "Cer(d18:1/16:0)"
        "Cer 16:0" → "Cer(d18:1/16:0)"
        "Cer d18:1/16:0" → "Cer(d18:1/16:0)"
        "SM 24:1;O2" → "SM(d18:1/24:1)"
        "SM(d18:1_24:1)" → "SM(d18:1/24:1)"
        "hex_cer_24_0" → "HexCer(d18:1/24:0)"
        "CER(d18:1/16:0)+H" → "Cer(d18:1/16:0)"

    Args:
        raw_name: Raw lipid name from file

    Returns:
        Canonical species name or None if cannot be normalized
    """
    if not raw_name or not isinstance(raw_name, str):
        return None

    raw = raw_name.strip()
    if not raw:
        return None

    # If already in canonical format, return as-is (with cleanup)
    if "(" in raw and ")" in raw:
        # Clean up: remove adducts, fix separators
        canonical = re.sub(r"[\+\-].*$", "", raw)  # Remove adducts (+H, -H2O, etc.)
        canonical = canonical.replace("_", "/")  # Fix underscore separators
        canonical = canonical.replace(" ", "")  # Remove spaces
        # Ensure proper format: Class(d18:1/chain) or Class(d18:1_chain) -> Class(d18:1/chain)
        if re.match(r"^[A-Za-z]+\(d\d+:\d+[/_]\d+:\d+\)$", canonical):
            # Normalize separator to /
            canonical = canonical.replace("_", "/")
            logger.debug(
                "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
                raw_name,
                canonical,
            )
            return canonical

    raw_lower = raw.lower().strip()

    # Remove adducts and charge info
    raw_lower = re.sub(r"[\+\-].*$", "", raw_lower)
    raw_lower = re.sub(r";.*$", "", raw_lower)  # Remove modifications like ;O2

    # Class mapping
    class_map = {
        "cer": "Cer",
        "ceramide": "Cer",
        "sm": "SM",
        "sphingomyelin": "SM",
        "hexcer": "HexCer",
        "hexosylceramide": "HexCer",
        "laccer": "LacCer",
        "lactosylceramide": "LacCer",
        "glccer": "GlcCer",
        "glucosylceramide": "GlcCer",
    }

    # Try to match various formats
    # Format 1: CER 16:0, SM 24:1, etc.
    simple_match = re.match(
        r"([a-z]+)\s*([cC]?)(\d+):(\d+)", raw_lower
    )
    if simple_match:
        class_name = simple_match.group(1)
        chain = f"{simple_match.group(3)}:{simple_match.group(4)}"
        canonical_class = class_map.get(class_name, class_name.capitalize())
        result = f"{canonical_class}(d18:1/{chain})"
        logger.info(
            "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
            raw_name,
            result,
        )
        return result

    # Format 2: Cer d18:1/16:0, SM d18:1/24:1, etc.
    with_backbone = re.match(
        r"([a-z]+)\s*d(\d+):(\d+)[_/](\d+):(\d+)", raw_lower
    )
    if with_backbone:
        class_name = with_backbone.group(1)
        backbone = f"d{with_backbone.group(2)}:{with_backbone.group(3)}"
        chain = f"{with_backbone.group(4)}:{with_backbone.group(5)}"
        canonical_class = class_map.get(class_name, class_name.capitalize())
        result = f"{canonical_class}({backbone}/{chain})"
        logger.info(
            "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
            raw_name,
            result,
        )
        return result

    # Format 3: hex_cer_24_0, sm_18_1, etc.
    underscore_format = re.match(
        r"([a-z_]+)_(\d+)_(\d+)", raw_lower
    )
    if underscore_format:
        class_part = underscore_format.group(1)
        chain = f"{underscore_format.group(2)}:{underscore_format.group(3)}"
        # Extract class name - check for hex, lac, glc prefixes first
        if "hex" in class_part:
            result = f"HexCer(d18:1/{chain})"
        elif "lac" in class_part:
            result = f"LacCer(d18:1/{chain})"
        elif "glc" in class_part:
            result = f"GlcCer(d18:1/{chain})"
        else:
            # Check other class mappings
            for key, canonical_class in class_map.items():
                if key in class_part:
                    result = f"{canonical_class}(d18:1/{chain})"
                    break
            else:
                # Default to Cer if no match
                result = f"Cer(d18:1/{chain})"
        logger.info(
            "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
            raw_name,
            result,
        )
        return result

    # Try the existing map_raw_lipid_to_canonical_species as fallback
    from amprenta_rag.ingestion.signature_matching import (
        map_raw_lipid_to_canonical_species)
    result = map_raw_lipid_to_canonical_species(raw_name)
    if result:
        logger.info(
            "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
            raw_name,
            result,
        )
        return result

    logger.warning(
        "[INGEST][LIPIDOMICS] WARNING: Could not normalize raw lipid '%s'",
        raw_name,
    )
    return None


def extract_species_from_file(file_path: str) -> tuple[Set[str], int]:
    """
    Extract and normalize lipid species from a CSV/TSV file.

    Args:
        file_path: Path to the lipidomics file

    Returns:
        Tuple of (set of normalized species, total number of raw rows)
    """
    try:
        # Try to detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
    except Exception as e:
        logger.error(
            "[INGEST][LIPIDOMICS] Error reading file %s: %r",
            file_path,
            e,
        )
        raise

    try:
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
    except Exception as e:
        logger.error(
            "[INGEST][LIPIDOMICS] Error parsing CSV/TSV file %s: %r",
            file_path,
            e,
        )
        raise

    if df.empty:
        logger.warning("[INGEST][LIPIDOMICS] File %s is empty", file_path)
        return set(), 0

    # Detect lipid identity column
    lipid_column = None
    candidate_names = ["species", "lipid", "Lipid", "Name", "Molecule", "compound", "metabolite"]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            lipid_column = col
            break

    if not lipid_column:
        # Try case-insensitive match
        for col in df.columns:
            if any(cand.lower() in col.lower() for cand in candidate_names):
                lipid_column = col
                break

    if not lipid_column:
        logger.error(
            "[INGEST][LIPIDOMICS] Could not find lipid identity column in file %s. "
            "Available columns: %s",
            file_path,
            list(df.columns),
        )
        raise ValueError(f"Could not find lipid identity column in {file_path}")

    logger.info(
        "[INGEST][LIPIDOMICS] Using column '%s' for lipid identity in file %s",
        lipid_column,
        file_path,
    )

    # Extract and normalize species
    species_set: Set[str] = set()
    total_rows = len(df)

    for idx, row in df.iterrows():
        raw_name = row.get(lipid_column)
        if pd.isna(raw_name):
            continue

        raw_name_str = str(raw_name).strip()
        if not raw_name_str:
            continue

        normalized = normalize_lipid_species(raw_name_str)
        if normalized:
            species_set.add(normalized)
        else:
            # Keep raw name if normalization fails (for now)
            species_set.add(raw_name_str)

    logger.info(
        "[INGEST][LIPIDOMICS] Extracted %d unique species from %d rows in file %s",
        len(species_set),
        total_rows,
        file_path,
    )

    return species_set, total_rows


def create_lipidomics_dataset_page(
    file_path: str,
    species_count: int,
    raw_rows: int,
) -> str:
    """
    Create a new Experimental Data Asset page for an internal lipidomics dataset.

    Args:
        file_path: Path to the lipidomics file
        species_count: Number of unique normalized species
        raw_rows: Number of raw rows in the file

    Returns:
        Created Notion page ID
    """
    return create_omics_dataset_page(
        file_path=file_path,
        omics_type="Lipidomics",
        entity_count=species_count,
        raw_rows=raw_rows,
        entity_name="species",
    )




def embed_lipidomics_dataset(
    page_id: str,
    dataset_name: str,
    species: Set[str],
    signature_matches: Optional[List[Any]] = None,
) -> None:
    """
    Embed lipidomics dataset into Pinecone for RAG.

    Args:
        page_id: Notion page ID
        dataset_name: Dataset name
        species: Set of normalized species
        signature_matches: Optional list of signature match results
    """
    cfg = get_config()

    try:
        # Build text representation
        text_parts = [
            f"Internal Lipidomics Dataset: {dataset_name}",
            f"Data Origin: Internal – Amprenta",
            f"Dataset Source Type: Processed table",
            "",
            f"Contains {len(species)} normalized lipid species:",
        ]

        # Add species list (truncate if too long)
        species_list = sorted(list(species))
        if len(species_list) > 50:
            text_parts.append(", ".join(species_list[:50]))
            text_parts.append(f"... and {len(species_list) - 50} more species")
        else:
            text_parts.append(", ".join(species_list))

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
                "[INGEST][LIPIDOMICS] No chunks generated for dataset %s",
                page_id,
            )
            return

        embeddings = embed_texts(chunks)

        # Upsert to Pinecone
        index = get_pinecone_index()
        vectors = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{page_id.replace('-', '')}_lipidomics_chunk_{order:03d}"

            meta = {
                "source": "Dataset",
                "source_type": "Dataset",
                "dataset_page_id": page_id.replace("-", ""),
                "dataset_name": dataset_name,
                "data_origin": "Internal – Amprenta",
                "snippet": chunk[:300],
            }

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

        # Update Notion with embedding metadata
        embedding_ids = [v["id"] for v in vectors]
        update_dataset_embedding_metadata(
            page_id=page_id,
            embedding_ids=embedding_ids,
            embedding_count=len(embedding_ids),
        )

        logger.info(
            "[INGEST][LIPIDOMICS] Embedded dataset %s to Pinecone (%d vectors)",
            page_id,
            len(vectors),
        )

    except Exception as e:
        logger.warning(
            "[INGEST][LIPIDOMICS] Error embedding dataset %s: %r",
            page_id,
            e,
        )
        # Don't raise - embedding is non-critical


def ingest_lipidomics_file(
    file_path: str,
    notion_page_id: Optional[str] = None,
    create_page: bool = False,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Ingest an internal lipidomics dataset file.

    Args:
        file_path: Path to the lipidomics file (CSV/TSV)
        notion_page_id: Existing Experimental Data Assets page ID to attach to
        create_page: If True and notion_page_id is None, create a new page
        program_ids: Optional list of Program page IDs to link
        experiment_ids: Optional list of Experiment page IDs to link

    Returns:
        The Notion page ID of the Experimental Data Asset
    """
    logger.info(
        "[INGEST][LIPIDOMICS] Starting ingestion of lipidomics file: %s",
        file_path,
    )

    # Validate file exists
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Lipidomics file not found: {file_path}")

    # Extract species from file
    species_set, raw_rows = extract_species_from_file(file_path)

    if not species_set:
        raise ValueError(f"No species extracted from file: {file_path}")

    # Determine page ID
    if notion_page_id:
        # Use existing page
        page_id = notion_page_id
        logger.info(
            "[INGEST][LIPIDOMICS] Using existing dataset page %s",
            page_id,
        )
        # Update summary
        try:
            page = fetch_dataset_page(page_id)
            page_props = page.get("properties", {}) or {}
            summary_prop = page_props.get("Summary", {}) or {}
            existing_summary = "".join(
                rt.get("plain_text", "")
                for rt in summary_prop.get("rich_text", [])
            )
            new_note = (
                f"\n\n---\n\n"
                f"Lipidomics file ingested: {file_path}\n"
                f"Contains {raw_rows} raw entries, {len(species_set)} normalized species."
            )
            combined_summary = existing_summary + new_note

            # Update page
            url = f"{get_config().notion.base_url}/pages/{page_id}"
            resp = requests.patch(
                url,
                headers=notion_headers(),
                json={"properties": {"Summary": {"rich_text": [{"text": {"content": combined_summary}}]}}},
                timeout=30,
            )
            resp.raise_for_status()
        except Exception as e:
            logger.warning(
                "[INGEST][LIPIDOMICS] Could not update summary for page %s: %r",
                page_id,
                e,
            )
    elif create_page:
        # Create new page
        page_id = create_lipidomics_dataset_page(
            file_path=file_path,
            species_count=len(species_set),
            raw_rows=raw_rows,
        )
    else:
        raise ValueError(
            "Either notion_page_id must be provided or create_page must be True"
        )

    # Attach file (noted in summary for now)
    attach_file_to_page(page_id, file_path, "Lipidomics")

    # Link to Programs and Experiments
    if program_ids or experiment_ids:
        link_to_programs_and_experiments(
            dataset_page_id=page_id,
            omics_type="Lipidomics",
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    # Link lipid species to Lipid Species DB
    try:
        from amprenta_rag.ingestion.feature_extraction import link_feature

        logger.info(
            "[INGEST][LIPIDOMICS] Linking %d lipid species to Lipid Species DB",
            len(species_set),
        )
        linked_count = 0
        for lipid in species_set:
            try:
                link_feature("lipid", lipid, page_id)
                linked_count += 1
            except Exception as e:
                logger.warning(
                    "[INGEST][LIPIDOMICS] Error linking lipid species '%s': %r",
                    lipid,
                    e,
                )
        logger.info(
            "[INGEST][LIPIDOMICS] Linked %d/%d lipid species to Lipid Species DB",
            linked_count,
            len(species_set),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][LIPIDOMICS] Feature linking skipped (error): %r",
            e,
        )

    # Score against signatures
    cfg = get_config()
    signature_matches = []
    if cfg.pipeline.enable_signature_scoring:
        logger.info(
            "[INGEST][LIPIDOMICS] Scoring dataset %s against signatures (%d species)",
            page_id,
            len(species_set),
        )

        matches = find_matching_signatures_for_dataset(
            dataset_species=species_set,  # Legacy support
            overlap_threshold=cfg.pipeline.signature_overlap_threshold,
            dataset_page_id=page_id,  # Multi-omics support
            omics_type="Lipidomics",
        )

        if matches:
            logger.info(
                "[INGEST][LIPIDOMICS] Found %d matching signature(s) for dataset %s",
                len(matches),
                page_id,
            )
            update_dataset_with_signature_matches(
                dataset_page_id=page_id,
                matches=matches,
            )
            signature_matches = matches
        else:
            logger.info(
                "[INGEST][LIPIDOMICS] No signature matches found for dataset %s",
                page_id,
            )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem
    embed_lipidomics_dataset(
        page_id=page_id,
        dataset_name=dataset_name,
        species=species_set,
        signature_matches=signature_matches,
    )

    # Auto-ingest linked Experiment pages
    if experiment_ids:
        from amprenta_rag.ingestion.experiments_ingestion import ingest_experiment

        for exp_id in experiment_ids:
            try:
                logger.info(
                    "[INGEST][EXPERIMENT] Auto-ingesting Experiment %s after lipidomics dataset %s created",
                    exp_id,
                    page_id,
                )
                ingest_experiment(
                    exp_page_id=exp_id,
                    parent_type="Experiment",
                    force=False,
                )
            except Exception as e:
                logger.warning(
                    "[INGEST][EXPERIMENT] Error auto-ingesting Experiment %s: %r",
                    exp_id,
                    e,
                )
                # Non-blocking - continue with other experiments

    logger.info(
        "[INGEST][LIPIDOMICS] Completed ingestion of file %s -> dataset page %s",
        file_path,
        page_id,
    )

    return page_id

