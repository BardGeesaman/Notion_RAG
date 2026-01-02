"""Ensembl VEP REST API client for variant annotation."""

from __future__ import annotations

import logging
import time
from typing import Optional

import requests
from pydantic import BaseModel
from tenacity import retry, wait_exponential, stop_after_attempt, retry_if_exception_type

logger = logging.getLogger(__name__)

VEP_BASE_URL = "https://rest.ensembl.org"
RATE_LIMIT_DELAY = 0.1  # 10 requests/second max


class VEPAnnotation(BaseModel):
    """Parsed VEP annotation result."""
    consequence: Optional[str] = None
    impact: Optional[str] = None
    symbol: Optional[str] = None
    gene_id: Optional[str] = None
    transcript_id: Optional[str] = None
    biotype: Optional[str] = None
    amino_acids: Optional[str] = None
    codons: Optional[str] = None
    protein_position: Optional[int] = None
    sift_prediction: Optional[str] = None
    sift_score: Optional[float] = None
    polyphen_prediction: Optional[str] = None
    polyphen_score: Optional[float] = None
    clin_sig: Optional[str] = None


@retry(
    wait=wait_exponential(multiplier=1, min=4, max=60),
    stop=stop_after_attempt(3),
    retry=retry_if_exception_type(requests.HTTPError),
)
def _vep_request(url: str, method: str = "GET", **kwargs) -> requests.Response:
    """Make VEP request with automatic retry on rate limit (429)."""
    time.sleep(RATE_LIMIT_DELAY)
    
    if method == "GET":
        response = requests.get(url, **kwargs)
    else:
        response = requests.post(url, **kwargs)
    
    response.raise_for_status()
    return response


def annotate_variant(
    chromosome: str,
    position: int,
    ref: str,
    alt: str,
    assembly: str = "GRCh38",
) -> Optional[VEPAnnotation]:
    """
    Annotate a single variant via VEP REST API.
    
    Args:
        chromosome: Chromosome name (e.g., "1", "chr1", "X")
        position: Genomic position (1-based)
        ref: Reference allele
        alt: Alternate allele
        assembly: Genome assembly (GRCh37 or GRCh38)
    
    Returns:
        VEPAnnotation object or None if annotation failed
    """
    # Normalize chromosome (remove 'chr' prefix if present)
    chrom = chromosome.replace("chr", "")
    
    # Build region string: chr:pos:pos/alt
    region = f"{chrom}:{position}:{position}/{alt}"
    
    url = f"{VEP_BASE_URL}/vep/human/region/{region}"
    headers = {"Content-Type": "application/json"}
    params = {
        "assembly": assembly,
        "SIFT": "b",  # Include SIFT with score
        "PolyPhen": "b",  # Include PolyPhen with score
    }
    
    try:
        response = _vep_request(url, method="GET", headers=headers, params=params, timeout=30)
        data = response.json()
        
        if data and len(data) > 0:
            return _parse_vep_response(data[0])
            
    except requests.HTTPError as e:
        logger.warning(f"VEP annotation failed for {region}: HTTP {e.response.status_code}")
    except Exception as e:
        logger.warning(f"VEP annotation failed for {region}: {e}")
    
    return None


def annotate_variants_batch(
    variants: list[dict],
    assembly: str = "GRCh38",
) -> list[Optional[VEPAnnotation]]:
    """
    Annotate multiple variants via VEP POST endpoint (max 200 per request).
    
    Args:
        variants: List of dicts with keys: chromosome, position, ref, alt
        assembly: Genome assembly (GRCh37 or GRCh38)
    
    Returns:
        List of VEPAnnotation objects (None for failed annotations)
    """
    if not variants:
        return []
    
    # VEP max 200 variants per request
    batch = variants[:200]
    
    url = f"{VEP_BASE_URL}/vep/human/region"
    headers = {"Content-Type": "application/json"}
    
    # Build HGVS-style notation for each variant
    hgvs_list = []
    for v in batch:
        chrom = str(v.get("chromosome", "")).replace("chr", "")
        pos = v.get("position")
        alt = v.get("alt") or v.get("alternate_allele")
        
        if chrom and pos and alt:
            hgvs_list.append(f"{chrom} {pos} {pos} {alt} 1")  # VEP POST format
    
    if not hgvs_list:
        return [None] * len(variants)
    
    payload = {"variants": hgvs_list}
    params = {
        "assembly": assembly,
        "SIFT": "b",
        "PolyPhen": "b",
    }
    
    try:
        response = _vep_request(
            url, 
            method="POST", 
            headers=headers, 
            json=payload, 
            params=params, 
            timeout=60
        )
        data = response.json()
        
        # Parse each result
        results = []
        for item in data:
            results.append(_parse_vep_response(item))
        
        # Pad with None if we got fewer results than input
        while len(results) < len(batch):
            results.append(None)
        
        return results
        
    except requests.HTTPError as e:
        logger.error(f"VEP batch annotation failed: HTTP {e.response.status_code}")
    except Exception as e:
        logger.error(f"VEP batch annotation failed: {e}")
    
    return [None] * len(batch)


def _parse_vep_response(data: dict) -> VEPAnnotation:
    """Parse VEP JSON response into VEPAnnotation object."""
    # Get transcript consequences (most detailed)
    tc = data.get("transcript_consequences", [])
    if not tc:
        tc = data.get("intergenic_consequences", [{}])
    
    # Use first (most severe) consequence
    most_severe = tc[0] if tc else {}
    
    # Extract SIFT info
    sift_pred = most_severe.get("sift_prediction")
    sift_score = most_severe.get("sift_score")
    
    # Extract PolyPhen info
    polyphen_pred = most_severe.get("polyphen_prediction")
    polyphen_score = most_severe.get("polyphen_score")
    
    return VEPAnnotation(
        consequence=data.get("most_severe_consequence"),
        impact=most_severe.get("impact"),
        symbol=most_severe.get("gene_symbol"),
        gene_id=most_severe.get("gene_id"),
        transcript_id=most_severe.get("transcript_id"),
        biotype=most_severe.get("biotype"),
        amino_acids=most_severe.get("amino_acids"),
        codons=most_severe.get("codons"),
        protein_position=most_severe.get("protein_start"),
        sift_prediction=sift_pred,
        sift_score=sift_score,
        polyphen_prediction=polyphen_pred,
        polyphen_score=polyphen_score,
        clin_sig=most_severe.get("clin_sig"),
    )


def get_vep_version() -> Optional[str]:
    """Get VEP REST API version for source tracking."""
    try:
        url = f"{VEP_BASE_URL}/info/rest"
        response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
        response.raise_for_status()
        data = response.json()
        return data.get("release")
    except Exception as e:
        logger.warning(f"Failed to get VEP version: {e}")
        return None
