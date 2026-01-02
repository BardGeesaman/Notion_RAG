"""External API integrations for target enrichment."""
import requests
import logging
from typing import Optional
from functools import lru_cache

logger = logging.getLogger(__name__)

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"


@lru_cache(maxsize=100)
def fetch_uniprot_data(uniprot_id: str) -> Optional[dict]:
    """Fetch protein data from UniProt.
    
    Returns:
        dict with: gene_names, protein_name, function, organism, sequence_length,
        subcellular_location, disease_associations
    """
    try:
        url = f"{UNIPROT_API}/{uniprot_id}.json"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            return {
                "uniprot_id": uniprot_id,
                "protein_name": data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
                "gene_names": [g.get("geneName", {}).get("value") for g in data.get("genes", [])],
                "organism": data.get("organism", {}).get("scientificName"),
                "sequence_length": data.get("sequence", {}).get("length"),
                "function": _extract_comment(data, "FUNCTION"),
                "subcellular_location": _extract_comment(data, "SUBCELLULAR LOCATION"),
            }
        logger.warning(f"UniProt API returned {resp.status_code} for {uniprot_id}")
        return None
    except Exception as e:
        logger.error(f"UniProt API error for {uniprot_id}: {e}")
        return None


def _extract_comment(data: dict, comment_type: str) -> Optional[str]:
    """Extract specific comment type from UniProt data."""
    for comment in data.get("comments", []):
        if comment.get("commentType") == comment_type:
            texts = comment.get("texts", [])
            if texts:
                return texts[0].get("value")
    return None


@lru_cache(maxsize=100)
def fetch_chembl_target(chembl_id: str) -> Optional[dict]:
    """Fetch target data from ChEMBL.
    
    Returns:
        dict with: target_type, organism, pref_name, target_components
    """
    try:
        url = f"{CHEMBL_API}/target/{chembl_id}.json"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            return {
                "chembl_id": chembl_id,
                "pref_name": data.get("pref_name"),
                "target_type": data.get("target_type"),
                "organism": data.get("organism"),
                "target_components": data.get("target_components", []),
            }
        return None
    except Exception as e:
        logger.error(f"ChEMBL API error for {chembl_id}: {e}")
        return None


@lru_cache(maxsize=100)
def fetch_chembl_drugs_for_target(chembl_id: str, limit: int = 20) -> list[dict]:
    """Fetch approved drugs and clinical compounds for a target.
    
    Returns:
        List of dicts with: molecule_chembl_id, pref_name, max_phase, 
        molecule_type, indication_class
    """
    try:
        url = f"{CHEMBL_API}/mechanism.json"
        params = {"target_chembl_id": chembl_id, "limit": limit}
        resp = requests.get(url, params=params, timeout=15)
        if resp.status_code == 200:
            mechanisms = resp.json().get("mechanisms", [])
            drugs = []
            for mech in mechanisms:
                drugs.append({
                    "molecule_chembl_id": mech.get("molecule_chembl_id"),
                    "mechanism_of_action": mech.get("mechanism_of_action"),
                    "action_type": mech.get("action_type"),
                    "max_phase": mech.get("max_phase"),
                })
            return drugs
        return []
    except Exception as e:
        logger.error(f"ChEMBL drugs API error for {chembl_id}: {e}")
        return []


def search_uniprot_by_gene(gene_symbol: str, organism: str = "Homo sapiens") -> list[dict]:
    """Search UniProt for entries matching gene symbol."""
    try:
        query = f"gene:{gene_symbol} AND organism_name:{organism}"
        url = f"{UNIPROT_API}/search"
        params = {"query": query, "format": "json", "size": 5}
        resp = requests.get(url, params=params, timeout=10)
        if resp.status_code == 200:
            results = resp.json().get("results", [])
            return [{"uniprot_id": r.get("primaryAccession"), 
                     "name": r.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")}
                    for r in results]
        return []
    except Exception as e:
        logger.error(f"UniProt search error for {gene_symbol}: {e}")
        return []


def search_chembl_by_gene(gene_symbol: str, organism: str = "Homo sapiens") -> list[dict]:
    """Search ChEMBL for targets matching gene symbol."""
    try:
        url = f"{CHEMBL_API}/target/search.json"
        params = {"q": gene_symbol, "limit": 10}
        resp = requests.get(url, params=params, timeout=10)
        if resp.status_code == 200:
            targets = resp.json().get("targets", [])
            human_targets = []
            for target in targets:
                if target.get("organism", "").lower() == organism.lower():
                    human_targets.append({
                        "chembl_id": target.get("target_chembl_id"),
                        "pref_name": target.get("pref_name"),
                        "target_type": target.get("target_type")
                    })
            return human_targets
        return []
    except Exception as e:
        logger.error(f"ChEMBL search error for {gene_symbol}: {e}")
        return []
