"""VCF file parser for creating GeneticVariant records."""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Max variants to parse in single operation
MAX_VARIANTS = 10000


def parse_vcf_file(vcf_path: Path) -> list[dict]:
    """
    Parse VCF file and return variant dictionaries.
    
    Args:
        vcf_path: Path to VCF or VCF.gz file
        
    Returns:
        List of variant dictionaries with keys:
        - chromosome, position, reference_allele, alternate_allele
        - quality, vcf_filter, gene (from INFO if available)
    """
    try:
        from pysam import VariantFile
    except ImportError:
        raise ImportError("pysam required for VCF parsing: pip install pysam")
    
    variants = []
    
    with VariantFile(str(vcf_path)) as vcf:
        for i, record in enumerate(vcf):
            if i >= MAX_VARIANTS:
                logger.warning(f"Truncated at {MAX_VARIANTS} variants")
                break
            
            # Extract gene from INFO field if present
            gene = None
            if "GENE" in record.info:
                gene = record.info["GENE"]
            elif "ANN" in record.info:
                # SnpEff annotation format
                ann = record.info["ANN"]
                if ann and isinstance(ann, (list, tuple)) and len(ann) > 0:
                    parts = str(ann[0]).split("|")
                    if len(parts) > 3:
                        gene = parts[3]
            
            for alt in record.alts or []:
                variants.append({
                    "chromosome": record.chrom,
                    "position": record.pos,
                    "reference_allele": record.ref,
                    "alternate_allele": alt,
                    "quality": record.qual,
                    "vcf_filter": ";".join(record.filter.keys()) if record.filter else "PASS",
                    "gene": gene,
                    "variant": f"{record.chrom}:{record.pos}{record.ref}>{alt}",
                })
    
    logger.info(f"Parsed {len(variants)} variants from {vcf_path}")
    return variants


def parse_vcf_bytes(content: bytes) -> list[dict]:
    """Parse VCF content from bytes (for file uploads)."""
    import tempfile
    
    # Write to temp file (pysam requires file path)
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as f:
        f.write(content)
        temp_path = Path(f.name)
    
    try:
        return parse_vcf_file(temp_path)
    finally:
        temp_path.unlink()
