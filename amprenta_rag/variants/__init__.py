"""Variant interpretation utilities (VEP parsing, ClinVar integration)."""

from .vep_parser import parse_vep_tsv
from .clinvar_loader import download_clinvar, parse_clinvar
from .clinvar_matcher import match_variants_to_clinvar
from .gene_burden import compute_gene_burden
from .ingestion_service import ingest_vep_tsv

__all__ = [
    "parse_vep_tsv",
    "download_clinvar",
    "parse_clinvar",
    "match_variants_to_clinvar",
    "compute_gene_burden",
    "ingest_vep_tsv",
]


