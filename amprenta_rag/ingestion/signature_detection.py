# amprenta_rag/ingestion/signature_detection.py

"""
Signature detection and extraction from various content sources.

This module detects lipid signature definitions in:
- Text content (literature, emails, experiments)
- Attached files (TSV/CSV)
- Links to signature repositories
- Embedded tables
- mwTab datasets

Disease-agnostic and source-agnostic.
"""

from __future__ import annotations

import csv
import io
import re
from typing import Dict, Any, List, Optional, Set, Tuple
from pathlib import Path

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Signature-related keywords (disease-agnostic)
SIGNATURE_KEYWORDS = [
    "signature",
    "panel",
    "lipid panel",
    "ceramide signature",
    "sphingolipid signature",
    "metabolic signature",
    "biomarker panel",
    "core ceramides",
    "six ceramide",
    "ceramide profile",
    "lipid profile",
]

# Patterns that indicate signature definitions
SIGNATURE_PATTERNS = [
    r'signature\s+components?',
    r'lipid\s+panel',
    r'ceramide\s+signature',
    r'sphingolipid\s+signature',
    r'biomarker\s+panel',
]


def detect_signature_keywords(text: str) -> bool:
    """
    Detect if text contains signature-related keywords.
    
    Args:
        text: Text content to scan
        
    Returns:
        True if signature keywords found, False otherwise
    """
    if not text:
        return False
    
    text_lower = text.lower()
    
    for keyword in SIGNATURE_KEYWORDS:
        if keyword in text_lower:
            return True
    
    return False


def extract_embedded_signature_table(text: str) -> Optional[List[Dict[str, str]]]:
    """
    Extract signature-like tables from text content.
    
    Looks for tabular structures with lipid names and directions/weights.
    
    Args:
        text: Text content containing potential signature tables
        
    Returns:
        List of dictionaries with 'species', 'direction', 'weight' keys, or None
    """
    if not text:
        return None
    
    lines = text.split('\n')
    table_rows: List[Dict[str, str]] = []
    
    # Look for table-like structures (tab-separated or multiple spaces)
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        
        # Check if line looks like a table row with lipid + direction
        # Pattern: "Cer(d18:1/16:0) ↑ 1.0" or "Ceramide\t↑\t1.0" or "Cer(d18:1/16:0) ↑"
        # Split on tabs, multiple spaces, or single spaces (more flexible)
        parts = re.split(r'\s{2,}|\t|\s', line)
        parts = [p.strip() for p in parts if p.strip()]  # Remove empty parts
        
        if len(parts) >= 2:
            # Check if first part looks like a lipid name
            first_part = parts[0].strip()
            if any(keyword in first_part.lower() for keyword in ['cer', 'sm', 'hexcer', 'lacc', 'sphing']):
                row: Dict[str, str] = {"species": first_part}
                
                # Look for direction in remaining parts
                for part in parts[1:]:
                    part = part.strip()
                    if part in ['↑', '↓', 'up', 'down', 'increased', 'decreased']:
                        row["direction"] = part
                    elif re.match(r'^\d+\.?\d*$', part):
                        row["weight"] = part
                
                if len(row) >= 2:  # At least species + direction/weight
                    table_rows.append(row)
    
    if len(table_rows) >= 2:  # Need at least 2 rows to be considered a signature
        logger.debug(
            "[INGEST][SIGNATURES] Extracted %d rows from embedded table",
            len(table_rows),
        )
        return table_rows
    
    return None


def find_attached_signature_files(content_paths: List[Path]) -> List[Path]:
    """
    Find attached signature TSV/CSV files in a list of file paths.
    
    Args:
        content_paths: List of file paths (e.g., from attachments)
        
    Returns:
        List of paths to potential signature files
    """
    signature_files: List[Path] = []
    
    for path in content_paths:
        if not path.exists() or not path.is_file():
            continue
        
        # Check file extension
        if path.suffix.lower() in ['.tsv', '.csv']:
            # Check filename for signature keywords
            filename_lower = path.name.lower()
            if any(keyword in filename_lower for keyword in ['signature', 'panel', 'ceramide', 'lipid']):
                signature_files.append(path)
                logger.debug(
                    "[INGEST][SIGNATURES] Found potential signature file: %s",
                    path.name,
                )
    
    return signature_files


def extract_signature_from_text_table(text: str) -> Optional[Dict[str, Any]]:
    """
    Extract signature definition from text table content.
    
    Attempts to parse tabular data that looks like a signature.
    
    Args:
        text: Text content with potential signature table
        
    Returns:
        Dictionary with signature components or None
    """
    # Try to parse as TSV/CSV
    for delimiter in ['\t', ',']:
        try:
            reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
            fieldnames = reader.fieldnames
            if not fieldnames:
                continue
            
            # Check if this looks like a signature table
            # Should have columns like: species, metabolite, lipid + direction, weight
            has_species_col = any(
                col.lower() in ['species', 'metabolite', 'lipid', 'name', 'compound']
                for col in fieldnames
            )
            has_direction_col = any(
                col.lower() in ['direction', 'dir', 'change', 'trend']
                for col in fieldnames
            )
            
            if not has_species_col:
                continue
            
            components: List[Dict[str, str]] = []
            
            for row in reader:
                # Find species column
                species = None
                for col in fieldnames:
                    col_lower = col.lower()
                    if col_lower in ['species', 'metabolite', 'lipid', 'name', 'compound']:
                        species = row.get(col, "").strip()
                        break
                
                if not species:
                    continue
                
                component: Dict[str, str] = {"species": species}
                
                # Find direction column
                for col in fieldnames:
                    col_lower = col.lower()
                    if col_lower in ['direction', 'dir', 'change', 'trend']:
                        direction = row.get(col, "").strip()
                        if direction:
                            component["direction"] = direction
                        break
                
                # Find weight column
                for col in fieldnames:
                    col_lower = col.lower()
                    if col_lower in ['weight', 'w', 'importance']:
                        weight = row.get(col, "").strip()
                        if weight:
                            component["weight"] = weight
                        break
                
                components.append(component)
            
            if len(components) >= 2:  # Need at least 2 components
                logger.debug(
                    "[INGEST][SIGNATURES] Extracted signature from text table with %d components",
                    len(components),
                )
                return {"components": components}
        
        except Exception:
            continue
    
    return None


def detect_signature_urls(text: str) -> List[str]:
    """
    Detect URLs pointing to signature repositories or files.
    
    Args:
        text: Text content to scan
        
    Returns:
        List of URLs that might point to signature files
    """
    urls: List[str] = []
    
    # URL pattern
    url_pattern = re.compile(
        r'https?://[^\s<>"{}|\\^`\[\]]+\.(?:tsv|csv|json|txt)',
        re.IGNORECASE,
    )
    
    matches = url_pattern.findall(text)
    for url in matches:
        # Check if URL looks signature-related
        url_lower = url.lower()
        if any(keyword in url_lower for keyword in ['signature', 'panel', 'ceramide', 'lipid']):
            urls.append(url)
            logger.debug(
                "[INGEST][SIGNATURES] Found potential signature URL: %s",
                url,
            )
    
    return urls


def infer_signature_metadata_from_source(
    source_type: str,
    source_metadata: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Infer signature metadata from source page metadata.
    
    Disease-agnostic: extracts disease, matrix, etc. from source metadata.
    
    Args:
        source_type: "literature", "dataset", "email", "experiment"
        source_metadata: Metadata dict from source page
        
    Returns:
        Dictionary with inferred metadata (disease_context, matrix, etc.)
    """
    metadata: Dict[str, Any] = {}
    
    # Extract disease context
    diseases = source_metadata.get("diseases", []) or []
    if diseases:
        metadata["disease_context"] = diseases
    
    # Extract matrix
    matrix = source_metadata.get("matrix", []) or []
    if matrix:
        metadata["matrix"] = matrix
    
    # Infer signature type from source type
    if source_type == "dataset":
        metadata["signature_type"] = "Open Dataset"
    elif source_type == "literature":
        metadata["signature_type"] = "Literature-derived"
    else:
        metadata["signature_type"] = "Literature-derived"
    
    return metadata


def save_extracted_signature_to_file(
    components: List[Dict[str, str]],
    output_dir: Path,
    signature_name: str,
) -> Optional[Path]:
    """
    Save extracted signature components to a temporary TSV file for ingestion.
    
    Args:
        components: List of component dictionaries
        output_dir: Directory to save the file
        signature_name: Base name for the signature file
        
    Returns:
        Path to saved file or None if save failed
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate filename
    safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', signature_name)[:50]
    output_path = output_dir / f"{safe_name}.tsv"
    
    try:
        with output_path.open("w", encoding="utf-8", newline="") as f:
            # Determine columns
            all_keys: Set[str] = set()
            for comp in components:
                all_keys.update(comp.keys())
            
            # Standardize column names
            fieldnames = ["species"]
            if "direction" in all_keys:
                fieldnames.append("direction")
            if "weight" in all_keys:
                fieldnames.append("weight")
            
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            
            for comp in components:
                writer.writerow(comp)
        
        logger.info(
            "[INGEST][SIGNATURES] Saved extracted signature to: %s",
            output_path,
        )
        return output_path
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error saving signature to file: %r",
            e,
        )
        return None

