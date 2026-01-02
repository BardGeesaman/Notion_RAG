"""BAM/CRAM file parser for alignment viewing."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

from pydantic import BaseModel

logger = logging.getLogger(__name__)

# Storage configuration
ALIGNMENT_STORAGE_PATH = os.getenv("ALIGNMENT_STORAGE_PATH", "./data/alignments")


class AlignmentStats(BaseModel):
    """Statistics for an alignment file."""
    total_reads: int
    mapped_reads: int
    unmapped_reads: int
    duplicate_rate: float
    mean_coverage: Optional[float] = None


class ReadRecord(BaseModel):
    """Single read from alignment file."""
    query_name: str
    flag: int
    reference_name: str
    reference_start: int
    mapping_quality: int
    cigar: str
    sequence: str
    is_reverse: bool
    is_duplicate: bool
    is_paired: bool


def check_index_exists(bam_path: Path) -> tuple[bool, Optional[Path]]:
    """
    Check if index file exists for BAM/CRAM.
    
    Returns:
        (has_index, index_path)
    """
    # Check common index patterns
    for suffix in [".bai", ".crai", f"{bam_path.suffix}.bai", f"{bam_path.suffix}.crai"]:
        index_path = bam_path.with_suffix(bam_path.suffix + suffix.lstrip(bam_path.suffix))
        if index_path.exists():
            return True, index_path
    
    # Also check {filename}.bai pattern
    bai_path = Path(str(bam_path) + ".bai")
    if bai_path.exists():
        return True, bai_path
    
    return False, None


def parse_bam_header(bam_path: Path) -> dict:
    """
    Extract header info from BAM/CRAM file.
    
    Returns:
        Dict with: reference_sequences, read_groups, programs, reference_genome
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam required: pip install pysam")
    
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        header = bam.header.to_dict()
        
        # Extract reference sequences
        references = []
        for sq in header.get("SQ", []):
            references.append({
                "name": sq.get("SN"),
                "length": sq.get("LN"),
            })
        
        # Extract read groups
        read_groups = []
        for rg in header.get("RG", []):
            read_groups.append({
                "id": rg.get("ID"),
                "sample": rg.get("SM"),
                "library": rg.get("LB"),
                "platform": rg.get("PL"),
            })
        
        # Try to detect reference genome from header
        reference_genome = None
        for sq in header.get("SQ", []):
            if "AS" in sq:  # Assembly field
                reference_genome = sq["AS"]
                break
        
        return {
            "reference_sequences": references,
            "read_groups": read_groups,
            "num_references": len(references),
            "reference_genome": reference_genome,
        }


def get_alignment_stats(bam_path: Path) -> AlignmentStats:
    """
    Calculate alignment statistics.
    
    Note: For large files, consider using samtools flagstat for speed.
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam required: pip install pysam")
    
    total = 0
    mapped = 0
    unmapped = 0
    duplicates = 0
    
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            total += 1
            if read.is_unmapped:
                unmapped += 1
            else:
                mapped += 1
            if read.is_duplicate:
                duplicates += 1
    
    duplicate_rate = duplicates / total if total > 0 else 0.0
    
    return AlignmentStats(
        total_reads=total,
        mapped_reads=mapped,
        unmapped_reads=unmapped,
        duplicate_rate=duplicate_rate,
    )


def fetch_reads(
    bam_path: Path,
    region: str,
    offset: int = 0,
    limit: int = 100,
) -> list[ReadRecord]:
    """
    Fetch reads from a region with pagination.
    
    Args:
        bam_path: Path to BAM/CRAM file
        region: Region in chr:start-end format
        offset: Number of reads to skip
        limit: Maximum reads to return (max 1000)
    
    Returns:
        List of ReadRecord objects
    
    Raises:
        ValueError: If no index file or invalid region
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam required: pip install pysam")
    
    has_index, _ = check_index_exists(bam_path)
    if not has_index:
        raise ValueError("Index file required for region queries")
    
    limit = min(limit, 1000)  # Hard cap
    
    reads = []
    count = 0
    
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(region=region):
            if count < offset:
                count += 1
                continue
            
            if len(reads) >= limit:
                break
            
            reads.append(ReadRecord(
                query_name=read.query_name or "",
                flag=read.flag,
                reference_name=read.reference_name or "",
                reference_start=read.reference_start,
                mapping_quality=read.mapping_quality,
                cigar=read.cigarstring or "",
                sequence=read.query_sequence or "",
                is_reverse=read.is_reverse,
                is_duplicate=read.is_duplicate,
                is_paired=read.is_paired,
            ))
            count += 1
    
    return reads


def get_coverage(
    bam_path: Path,
    region: str,
    bin_size: int = 100,
) -> list[dict]:
    """
    Calculate coverage histogram for a region.
    
    Args:
        bam_path: Path to BAM/CRAM file
        region: Region in chr:start-end format
        bin_size: Size of each coverage bin
    
    Returns:
        List of {position, coverage} dicts
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam required: pip install pysam")
    
    has_index, _ = check_index_exists(bam_path)
    if not has_index:
        raise ValueError("Index file required for coverage queries")
    
    # Parse region
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))
    
    coverage_data = []
    
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for pos in range(start, end, bin_size):
            bin_end = min(pos + bin_size, end)
            
            # Count coverage in this bin
            count = 0
            for pileup in bam.pileup(chrom, pos, bin_end, truncate=True):
                count += pileup.n
            
            avg_coverage = count / bin_size if bin_size > 0 else 0
            
            coverage_data.append({
                "position": pos,
                "coverage": round(avg_coverage, 2),
            })
    
    return coverage_data
