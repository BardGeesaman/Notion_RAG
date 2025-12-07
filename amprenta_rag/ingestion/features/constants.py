"""
Constants for feature extraction.

Contains metabolite synonyms, amino acids, and nucleotide lists used
for text scanning and normalization.
"""

from __future__ import annotations

from typing import Dict, List

# Common metabolite synonyms mapping (expandable)
METABOLITE_SYNONYMS: Dict[str, str] = {
    "l-glutamate": "glutamate",
    "l-glutamic acid": "glutamate",
    "d-glutamate": "glutamate",
    "glutamic acid": "glutamate",
    "l-glutamine": "glutamine",
    "d-glutamine": "glutamine",
    "gln": "glutamine",
    "glu": "glutamate",
}

# Common amino acids (for text scanning)
AMINO_ACIDS: List[str] = [
    "alanine",
    "arginine",
    "asparagine",
    "aspartic acid",
    "cysteine",
    "glutamine",
    "glutamate",
    "glutamic acid",
    "glycine",
    "histidine",
    "isoleucine",
    "leucine",
    "lysine",
    "methionine",
    "phenylalanine",
    "proline",
    "serine",
    "threonine",
    "tryptophan",
    "tyrosine",
    "valine",
]

# Common nucleotides (for text scanning)
NUCLEOTIDES: List[str] = [
    "adenosine",
    "adenine",
    "guanosine",
    "guanine",
    "cytidine",
    "cytosine",
    "thymidine",
    "thymine",
    "uracil",
    "uridine",
    "atp",
    "adp",
    "amp",
    "gtp",
    "gdp",
    "gmp",
    "ctp",
    "cdp",
    "cmp",
    "utp",
    "udp",
    "ump",
]
