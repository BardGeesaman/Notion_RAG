"""
Feature type inference for multi-omics signatures.

Automatically detects whether a feature name represents a gene, protein,
metabolite, or lipid species using pattern matching and heuristics.
"""

import re

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Common metabolite names (expandable)
COMMON_METABOLITES = {
    "glucose",
    "lactate",
    "pyruvate",
    "citrate",
    "succinate",
    "fumarate",
    "malate",
    "oxaloacetate",
    "glutamate",
    "glutamine",
    "serine",
    "alanine",
    "aspartate",
    "asparagine",
    "lysine",
    "arginine",
    "proline",
    "valine",
    "leucine",
    "isoleucine",
    "methionine",
    "phenylalanine",
    "tyrosine",
    "tryptophan",
    "cysteine",
    "threonine",
    "histidine",
    "glycine",
    "choline",
    "creatine",
    "carnitine",
    "atp",
    "adp",
    "amp",
    "gtp",
    "gdp",
    "nad",
    "nadh",
    "nadp",
    "nadph",
}

# Lipid class prefixes
LIPID_CLASSES = [
    "cer",
    "ceramide",
    "sm",
    "sphingomyelin",
    "hexcer",
    "hexosylceramide",
    "laccer",
    "lactosylceramide",
    "glccer",
    "glucosylceramide",
    "pc",
    "pe",
    "ps",
    "pi",
    "pg",
    "pa",
    "tg",
    "dg",
    "mg",
    "lpc",
    "lpe",
    "lps",
    "lpi",
    "lpg",
    "lpa",
    "fa",
    "ffa",
]


def infer_feature_type(feature_name: str) -> str:
    """
    Infer feature type from feature name using heuristics.

    Detection order:
    1. Lipid species formats (Cer(, SM(, etc.)
    2. Metabolite names (common names, HMDB/KEGG IDs)
    3. Protein identifiers (UniProt IDs, FASTA formats)
    4. Gene identifiers (gene symbols, Ensembl IDs)
    5. Default fallback: gene

    Args:
        feature_name: Raw feature name/identifier

    Returns:
        Feature type: "lipid", "metabolite", "protein", or "gene"
    """
    if not feature_name or not isinstance(feature_name, str):
        logger.warning(
            "[SIGNATURE][FEATURE-TYPE] Empty or invalid feature name, defaulting to 'gene'"
        )
        return "gene"

    name = feature_name.strip()
    if not name:
        return "gene"

    name_lower = name.lower()

    # 1. Check for lipid species formats
    # Patterns: Cer(d18:1/16:0), SM(d18:1/24:1), LPC(16:0), etc.
    lipid_patterns = [
        r"^[A-Za-z]+\s*\(d?\d+:\d+[/_]\d+:\d+\)",  # Cer(d18:1/16:0)
        r"^[A-Za-z]+\s*\(\d+:\d+\)",  # LPC(16:0)
        r"^[A-Za-z]+\s*\d+:\d+",  # CER 16:0, SM 24:1
    ]

    for pattern in lipid_patterns:
        if re.match(pattern, name, re.IGNORECASE):
            logger.debug(
                "[SIGNATURE][FEATURE-TYPE] Detected lipid species format: '%s'",
                name,
            )
            return "lipid"

    # Check for lipid class prefixes
    name_words = re.split(r"[\s_\-\(\)]+", name_lower)
    for word in name_words:
        if word in LIPID_CLASSES:
            logger.debug(
                "[SIGNATURE][FEATURE-TYPE] Detected lipid class keyword: '%s' in '%s'",
                word,
                name,
            )
            return "lipid"

    # 2. Check for metabolite identifiers
    # HMDB: HMDB0000123, KEGG: C00001, CHEBI: CHEBI:12345
    metabolite_id_patterns = [
        r"^HMDB\d+",
        r"^C\d{5}",  # KEGG compound IDs
        r"^CHEBI:\d+",
        r"^PubChem:\d+",
        r"^CAS:\d+",
    ]

    for pattern in metabolite_id_patterns:
        if re.match(pattern, name, re.IGNORECASE):
            logger.debug(
                "[SIGNATURE][FEATURE-TYPE] Detected metabolite ID format: '%s'",
                name,
            )
            return "metabolite"

    # Check common metabolite names
    if name_lower in COMMON_METABOLITES:
        logger.debug(
            "[SIGNATURE][FEATURE-TYPE] Detected common metabolite name: '%s'",
            name,
        )
        return "metabolite"

    # 3. Check for protein identifiers (BEFORE gene symbols to avoid conflicts)
    # FASTA format: sp|P12345|TP53_HUMAN or tr|... (check this first)
    if re.match(r"^(sp|tr|ref)\|", name, re.IGNORECASE):
        logger.debug(
            "[SIGNATURE][FEATURE-TYPE] Detected FASTA format: '%s'",
            name,
        )
        return "protein"

    # Protein with _HUMAN suffix: TP53_HUMAN, ACTB_HUMAN
    if re.search(r"_[A-Z]+$", name) and len(name.split("_")) == 2:
        last_part = name.split("_")[-1].upper()
        if last_part in ["HUMAN", "MOUSE", "RAT", "BOVIN", "PIG", "CHICK"]:
            logger.debug(
                "[SIGNATURE][FEATURE-TYPE] Detected protein with species suffix: '%s'",
                name,
            )
            return "protein"

    # UniProt IDs: P04637, Q9Y6K9 (typically 6-10 alphanumeric characters)
    # Pattern: Starts with letter, then number, then alphanumeric (6-10 total chars)
    # Common format: [A-Z][0-9][A-Z0-9]{4,8}
    if len(name) >= 6 and len(name) <= 10:
        # Check if it matches UniProt pattern (letter-number-alphanumeric)
        if re.match(r"^[A-Z][0-9][A-Z0-9]{4,8}$", name.upper()):
            # Additional check: should have at least one more letter after the number
            if re.search(r"[A-Z]", name.upper()[2:]):
                logger.debug(
                    "[SIGNATURE][FEATURE-TYPE] Detected UniProt ID format: '%s'",
                    name,
                )
                return "protein"

    # 4. Check for gene identifiers
    # Ensembl IDs: ENSG00000139618, ENSMUSG00000000001
    if re.match(r"^ENS[A-Z]*G\d+", name, re.IGNORECASE):
        logger.debug(
            "[SIGNATURE][FEATURE-TYPE] Detected Ensembl ID format: '%s'",
            name,
        )
        return "gene"

    # Gene symbols: typically 1-10 uppercase letters/numbers (TP53, ACTB, etc.)
    # But be careful - some metabolites/proteins can look like this too
    # Check for UniProt-like patterns first (6-10 chars, alphanumeric)
    is_uniprot_like = (
        len(name) >= 6
        and len(name) <= 10
        and re.match(r"^[A-Z0-9]+$", name.upper())
        and name_lower not in COMMON_METABOLITES
    )

    # If it looks like a UniProt ID (already checked above), skip gene symbol check
    # Otherwise, check for gene symbol pattern
    if not is_uniprot_like and re.match(r"^[A-Z]{1,10}\d*[A-Z]*$", name.upper()) and len(name) <= 10:
        # Could be gene symbol, but check if it's a known metabolite first
        if name_lower not in COMMON_METABOLITES:
            logger.debug(
                "[SIGNATURE][FEATURE-TYPE] Detected gene symbol pattern: '%s'",
                name,
            )
            return "gene"

    # 5. Default fallback: gene
    logger.debug(
        "[SIGNATURE][FEATURE-TYPE] Could not infer feature type for '%s', defaulting to 'gene'",
        name,
    )
    return "gene"

