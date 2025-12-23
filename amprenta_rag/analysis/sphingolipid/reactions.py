"""Hardcoded sphingolipid pathway reactions and key entities.

This is an MVP knowledge base to support imbalance scoring in lipidomics datasets.
"""

from __future__ import annotations

from typing import Dict, List


SPHINGOLIPID_ENZYMES: List[str] = [
    # Ceramide synthases
    "CERS1",
    "CERS2",
    "CERS3",
    "CERS4",
    "CERS5",
    "CERS6",
    # Sphingomyelinases
    "SMPD1",
    "SMPD3",
    "SMPD4",
    # Ceramidase
    "ASAH1",
    # Sphingomyelin synthases
    "SMS1",
    "SMS2",
    # Sphingosine kinases
    "SPHK1",
    "SPHK2",
    # S1P degradation
    "SGPL1",
    # Glycosylation
    "UGCG",
]


SPHINGOLIPID_METABOLITES: List[str] = [
    "Sphingosine",
    "S1P",
    "Cer",  # Ceramides (class)
    "SM",  # Sphingomyelin (class)
    "HexCer",  # Hexosylceramide (class)
]


# Simple reaction graph. Keys are reaction IDs.
SPHINGOLIPID_REACTIONS: Dict[str, Dict[str, object]] = {
    "ceramide_synthase": {
        "substrates": ["Sphingosine"],
        "products": ["Cer"],
        "enzymes": ["CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6"],
        "notes": "De novo ceramide synthesis (proxy: CERS expression).",
    },
    "sphingomyelinase": {
        "substrates": ["SM"],
        "products": ["Cer"],
        "enzymes": ["SMPD1", "SMPD3", "SMPD4"],
        "notes": "SM -> Cer via sphingomyelinases.",
    },
    "ceramidase": {
        "substrates": ["Cer"],
        "products": ["Sphingosine"],
        "enzymes": ["ASAH1"],
        "notes": "Cer -> Sphingosine via ceramidase.",
    },
    "glycosylation": {
        "substrates": ["Cer"],
        "products": ["HexCer"],
        "enzymes": ["UGCG"],
        "notes": "Cer -> GlcCer/HexCer via UGCG (class-level).",
    },
    "sm_synthase": {
        "substrates": ["Cer"],
        "products": ["SM"],
        "enzymes": ["SMS1", "SMS2"],
        "notes": "Cer -> SM via sphingomyelin synthases.",
    },
    "sphingosine_kinase": {
        "substrates": ["Sphingosine"],
        "products": ["S1P"],
        "enzymes": ["SPHK1", "SPHK2"],
        "notes": "Sphingosine -> S1P via sphingosine kinases.",
    },
    "s1p_lyase": {
        "substrates": ["S1P"],
        "products": ["degradation"],
        "enzymes": ["SGPL1"],
        "notes": "S1P degradation via SGPL1.",
    },
}


__all__ = ["SPHINGOLIPID_REACTIONS", "SPHINGOLIPID_ENZYMES", "SPHINGOLIPID_METABOLITES"]


