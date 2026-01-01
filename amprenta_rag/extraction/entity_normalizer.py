"""Entity normalization service for standardizing extracted entities."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)

# API endpoints
PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
UNIPROT_API = "https://rest.uniprot.org"

# Rate limiting
MIN_REQUEST_INTERVAL = 0.5
REQUEST_TIMEOUT = 30

# In-memory cache
_CACHE: Dict[str, dict] = {}


@dataclass
class GeneInfo:
    """Normalized gene information."""
    symbol: str
    name: Optional[str]
    entrez_id: Optional[str]
    uniprot_id: Optional[str]
    organism: Optional[str]


@dataclass
class CompoundInfo:
    """Normalized compound information."""
    name: str
    cid: Optional[int]  # PubChem CID
    inchi_key: Optional[str]
    smiles: Optional[str]
    molecular_formula: Optional[str]
    molecular_weight: Optional[float]


@dataclass
class DiseaseInfo:
    """Normalized disease information."""
    name: str
    mesh_id: Optional[str]
    doid: Optional[str]  # Disease Ontology ID
    synonyms: List[str]


class EntityNormalizer:
    """Normalize extracted entities to standard identifiers."""

    def __init__(self):
        """Initialize entity normalizer with rate limiting."""
        self._last_request_time: Dict[str, float] = {}

    def _rate_limit(self, service: str) -> None:
        """Apply rate limiting per service."""
        if service in self._last_request_time:
            elapsed = time.time() - self._last_request_time[service]
            if elapsed < MIN_REQUEST_INTERVAL:
                time.sleep(MIN_REQUEST_INTERVAL - elapsed)
        self._last_request_time[service] = time.time()

    def _get_cached(self, cache_key: str) -> Optional[dict]:
        """Get cached result if available."""
        return _CACHE.get(cache_key)

    def _set_cached(self, cache_key: str, value: dict) -> None:
        """Cache a result."""
        _CACHE[cache_key] = value

    def normalize_gene(self, name: str) -> Optional[GeneInfo]:
        """
        Normalize gene symbol via UniProt.
        
        Args:
            name: Gene symbol or name to normalize
            
        Returns:
            GeneInfo with standardized identifiers or None if not found
        """
        cache_key = f"gene:{name.lower()}"
        cached = self._get_cached(cache_key)
        if cached:
            return GeneInfo(**cached) if cached.get("symbol") else None

        try:
            self._rate_limit("uniprot")
            
            # Search UniProt
            url = f"{UNIPROT_API}/uniprotkb/search"
            params = {
                "query": f"gene:{name} AND organism_id:9606",  # Human
                "format": "json",
                "size": 1
            }
            
            response = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
            response.raise_for_status()
            data = response.json()
            
            results = data.get("results", [])
            if not results:
                self._set_cached(cache_key, {})
                return None
            
            entry = results[0]
            genes = entry.get("genes", [{}])
            gene_data = genes[0] if genes else {}
            
            info = GeneInfo(
                symbol=gene_data.get("geneName", {}).get("value", name),
                name=entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
                entrez_id=None,  # Would need NCBI lookup
                uniprot_id=entry.get("primaryAccession"),
                organism="Homo sapiens"
            )
            
            self._set_cached(cache_key, info.__dict__)
            logger.info(f"Normalized gene '{name}' -> {info.symbol} ({info.uniprot_id})")
            return info
            
        except Exception as e:
            logger.warning(f"Failed to normalize gene '{name}': {e}")
            self._set_cached(cache_key, {})
            return None

    def normalize_compound(self, name: str) -> Optional[CompoundInfo]:
        """
        Normalize compound via PubChem.
        
        Args:
            name: Compound name to normalize
            
        Returns:
            CompoundInfo with standardized identifiers or None if not found
        """
        cache_key = f"compound:{name.lower()}"
        cached = self._get_cached(cache_key)
        if cached:
            return CompoundInfo(**cached) if cached.get("name") else None

        try:
            self._rate_limit("pubchem")
            
            # Search PubChem by name
            url = f"{PUBCHEM_API}/compound/name/{requests.utils.quote(name)}/JSON"
            
            response = requests.get(url, timeout=REQUEST_TIMEOUT)
            if response.status_code == 404:
                self._set_cached(cache_key, {})
                return None
            response.raise_for_status()
            
            data = response.json()
            compounds = data.get("PC_Compounds", [])
            if not compounds:
                self._set_cached(cache_key, {})
                return None
            
            compound = compounds[0]
            cid = compound.get("id", {}).get("id", {}).get("cid")
            
            # Get properties
            props = {p.get("urn", {}).get("label"): p.get("value", {}) 
                     for p in compound.get("props", [])}
            
            info = CompoundInfo(
                name=name,
                cid=cid,
                inchi_key=props.get("InChIKey", {}).get("sval"),
                smiles=props.get("SMILES", {}).get("sval"),
                molecular_formula=props.get("Molecular Formula", {}).get("sval"),
                molecular_weight=props.get("Molecular Weight", {}).get("fval")
            )
            
            self._set_cached(cache_key, info.__dict__)
            logger.info(f"Normalized compound '{name}' -> CID:{info.cid}")
            return info
            
        except Exception as e:
            logger.warning(f"Failed to normalize compound '{name}': {e}")
            self._set_cached(cache_key, {})
            return None

    def normalize_disease(self, name: str) -> Optional[DiseaseInfo]:
        """
        Normalize disease name.
        
        Note: Uses simple string matching. For production, integrate with
        MeSH or Disease Ontology APIs.
        
        Args:
            name: Disease name to normalize
            
        Returns:
            DiseaseInfo with standardized name or None
        """
        cache_key = f"disease:{name.lower()}"
        cached = self._get_cached(cache_key)
        if cached:
            return DiseaseInfo(**cached) if cached.get("name") else None

        # Simple normalization (case standardization)
        # TODO: Integrate with MeSH API for proper normalization
        normalized_name = name.strip().title()
        
        info = DiseaseInfo(
            name=normalized_name,
            mesh_id=None,  # Would need MeSH API
            doid=None,  # Would need Disease Ontology API
            synonyms=[]
        )
        
        self._set_cached(cache_key, info.__dict__)
        logger.info(f"Normalized disease '{name}' -> '{info.name}'")
        return info

    def normalize_batch(
        self, 
        genes: List[str] = None,
        compounds: List[str] = None,
        diseases: List[str] = None
    ) -> Dict[str, List]:
        """
        Normalize multiple entities in batch.
        
        Args:
            genes: List of gene names to normalize
            compounds: List of compound names to normalize
            diseases: List of disease names to normalize
            
        Returns:
            Dictionary with normalized results per entity type
        """
        results = {
            "genes": [],
            "compounds": [],
            "diseases": []
        }
        
        for gene in (genes or []):
            info = self.normalize_gene(gene)
            if info:
                results["genes"].append(info.__dict__)
        
        for compound in (compounds or []):
            info = self.normalize_compound(compound)
            if info:
                results["compounds"].append(info.__dict__)
        
        for disease in (diseases or []):
            info = self.normalize_disease(disease)
            if info:
                results["diseases"].append(info.__dict__)
        
        return results


def clear_cache() -> None:
    """Clear the normalization cache."""
    _CACHE.clear()
    logger.info("Entity normalization cache cleared")
