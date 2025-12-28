"""AI-powered dataset finder for cross-repository search."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Set

import httpx
from difflib import SequenceMatcher

import logging

logger = logging.getLogger(__name__)


@dataclass
class DatasetResult:
    """Single dataset search result."""
    
    accession: str
    title: str
    description: str
    source: str  # "geo", "arrayexpress", "metabolomics_workbench"
    species: Optional[str] = None
    tissue: Optional[str] = None
    disease: Optional[str] = None
    assay_type: Optional[str] = None
    sample_count: Optional[int] = None
    url: Optional[str] = None
    score: float = 0.0  # Relevance score


@dataclass
class DatasetFinderResult:
    """Complete dataset finder result."""
    
    query: str
    extracted_terms: Dict[str, List[str]]
    results: List[DatasetResult]
    total_found: int
    sources_searched: List[str]
    sources_failed: List[str]


def extract_search_terms(query: str) -> Dict[str, List[str]]:
    """Extract structured search terms from natural language query using simple heuristics."""
    query_lower = query.lower()
    
    # Disease terms - common patterns
    diseases = []
    disease_patterns = [
        r'\b(cancer|carcinoma|tumor|tumour)\b',
        r'\b(diabetes|diabetic)\b',
        r'\b(alzheimer|parkinson|huntington)\b',
        r'\b(heart disease|cardiovascular)\b',
        r'\b(stroke|ischemic|hemorrhagic)\b',
        r'\b(depression|anxiety|bipolar)\b',
        r'\b(arthritis|osteoarthritis|rheumatoid)\b',
        r'\b(asthma|copd|lung disease)\b',
        r'\b(covid|sars|coronavirus)\b',
        r'\b(inflammatory|inflammation)\b',
        r'\b(infection|infectious|bacterial|viral)\b',
    ]
    
    for pattern in disease_patterns:
        matches = re.findall(pattern, query_lower)
        diseases.extend(matches)
    
    # Tissue/organ terms
    tissues = []
    tissue_patterns = [
        r'\b(brain|cerebral|neural|neuronal)\b',
        r'\b(heart|cardiac|myocardial)\b',
        r'\b(liver|hepatic)\b',
        r'\b(kidney|renal)\b',
        r'\b(lung|pulmonary|respiratory)\b',
        r'\b(muscle|skeletal muscle|smooth muscle)\b',
        r'\b(blood|plasma|serum)\b',
        r'\b(skin|dermal|epidermal)\b',
        r'\b(bone|skeletal|osteoblast)\b',
        r'\b(adipose|fat|lipid)\b',
        r'\b(breast|mammary)\b',
        r'\b(prostate|testicular|ovarian)\b',
        r'\b(colon|intestinal|gastric)\b',
    ]
    
    for pattern in tissue_patterns:
        matches = re.findall(pattern, query_lower)
        tissues.extend(matches)
    
    # Species terms
    species = []
    species_patterns = [
        r'\b(human|homo sapiens)\b',
        r'\b(mouse|mice|mus musculus)\b',
        r'\b(rat|rattus norvegicus)\b',
        r'\b(zebrafish|danio rerio)\b',
        r'\b(drosophila|fruit fly)\b',
        r'\b(yeast|saccharomyces)\b',
        r'\b(arabidopsis|plant)\b',
    ]
    
    for pattern in species_patterns:
        matches = re.findall(pattern, query_lower)
        species.extend(matches)
    
    # Assay type terms
    assay_types = []
    assay_patterns = [
        r'\b(rna-seq|rnaseq|transcriptome)\b',
        r'\b(microarray|gene expression)\b',
        r'\b(chip-seq|chipseq|chromatin)\b',
        r'\b(single cell|scRNA|sc-RNA)\b',
        r'\b(metabolomics|metabolite)\b',
        r'\b(proteomics|protein)\b',
        r'\b(methylation|bisulfite)\b',
        r'\b(atac-seq|atacseq|accessibility)\b',
    ]
    
    for pattern in assay_patterns:
        matches = re.findall(pattern, query_lower)
        assay_types.extend(matches)
    
    return {
        "disease": list(set(diseases)),
        "tissue": list(set(tissues)),
        "species": list(set(species)),
        "assay_type": list(set(assay_types)),
    }


def search_geo(terms: Dict[str, List[str]], max_results: int = 50) -> List[DatasetResult]:
    """Search GEO database using Entrez API."""
    results = []
    
    try:
        # Build search query for GEO
        search_parts = []
        
        for term_type, term_list in terms.items():
            for term in term_list:
                if term_type == "species":
                    search_parts.append(f'"{term}"[Organism]')
                elif term_type == "disease":
                    search_parts.append(f'"{term}"[All Fields]')
                elif term_type == "tissue":
                    search_parts.append(f'"{term}"[All Fields]')
                elif term_type == "assay_type":
                    search_parts.append(f'"{term}"[All Fields]')
        
        if not search_parts:
            return results
        
        search_query = " AND ".join(search_parts[:3])  # Limit to avoid too complex queries
        search_query += " AND gse[Entry Type]"  # Only GSE entries
        
        # Search Entrez
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        search_url = f"{base_url}/esearch.fcgi"
        
        search_params = {
            "db": "gds",
            "term": search_query,
            "retmax": max_results,
            "retmode": "json",
        }
        
        with httpx.Client(timeout=30) as client:
            search_response = client.get(search_url, params=search_params)
            search_response.raise_for_status()
            search_data = search_response.json()
        
        ids = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not ids:
            logger.info("No GEO datasets found for query")
            return results
        
        # Fetch details for found IDs
        summary_url = f"{base_url}/esummary.fcgi"
        summary_params = {
            "db": "gds",
            "id": ",".join(ids),
            "retmode": "json",
        }
        
        with httpx.Client(timeout=30) as client:
            summary_response = client.get(summary_url, params=summary_params)
            summary_response.raise_for_status()
            summary_data = summary_response.json()
        
        # Parse results
        for uid, data in summary_data.get("result", {}).items():
            if uid == "uids":
                continue
                
            accession = data.get("accession", "")
            if not accession.startswith("GSE"):
                continue
            
            title = data.get("title", "")
            summary = data.get("summary", "")
            organism = data.get("taxon", "")
            sample_count = data.get("n_samples", 0)
            
            # Try to extract assay type from platform info
            assay_type = None
            if "rna-seq" in summary.lower() or "rnaseq" in summary.lower():
                assay_type = "RNA-seq"
            elif "microarray" in summary.lower():
                assay_type = "Microarray"
            elif "chip-seq" in summary.lower():
                assay_type = "ChIP-seq"
            
            result = DatasetResult(
                accession=accession,
                title=title,
                description=summary,
                source="geo",
                species=organism,
                assay_type=assay_type,
                sample_count=int(sample_count) if sample_count else None,
                url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
            )
            results.append(result)
        
        logger.info(f"Found {len(results)} datasets from GEO")
        
    except Exception as e:
        logger.error(f"GEO search failed: {e}")
    
    return results


def search_arrayexpress(terms: Dict[str, List[str]], max_results: int = 50) -> List[DatasetResult]:
    """Search ArrayExpress database using REST API."""
    results = []
    
    try:
        base_url = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments"
        
        # Build query parameters
        params = {
            "format": "json",
            "sortby": "releasedate",
            "sortorder": "descending",
            "pagesize": min(max_results, 100),
        }
        
        # Add search terms
        keywords = []
        for term_type, term_list in terms.items():
            keywords.extend(term_list)
        
        if keywords:
            params["keywords"] = " ".join(keywords[:5])  # Limit keywords
        
        # Add species filter if available
        if "species" in terms and terms["species"]:
            species_term = terms["species"][0]
            if "human" in species_term.lower():
                params["species"] = "homo sapiens"
            elif "mouse" in species_term.lower():
                params["species"] = "mus musculus"
        
        with httpx.Client(timeout=30) as client:
            response = client.get(base_url, params=params)
            response.raise_for_status()
            data = response.json()
        
        experiments = data.get("experiments", {}).get("experiment", [])
        
        for exp in experiments:
            accession = exp.get("accession", "")
            title = exp.get("name", "")
            description = exp.get("description", [""])[0] if exp.get("description") else ""
            
            # Extract organism
            organisms = exp.get("organism", [])
            species = organisms[0] if organisms else None
            
            # Extract assay types
            experiment_types = exp.get("experimenttype", [])
            assay_type = experiment_types[0] if experiment_types else None
            
            # Extract sample count
            sample_count = None
            if "samples" in exp:
                sample_count = len(exp["samples"])
            
            result = DatasetResult(
                accession=accession,
                title=title,
                description=description,
                source="arrayexpress",
                species=species,
                assay_type=assay_type,
                sample_count=sample_count,
                url=f"https://www.ebi.ac.uk/arrayexpress/experiments/{accession}/",
            )
            results.append(result)
        
        logger.info(f"Found {len(results)} datasets from ArrayExpress")
        
    except Exception as e:
        logger.error(f"ArrayExpress search failed: {e}")
    
    return results


def search_metabolomics_workbench(terms: Dict[str, List[str]], max_results: int = 50) -> List[DatasetResult]:
    """Search Metabolomics Workbench using REST API."""
    results = []
    
    try:
        base_url = "https://www.metabolomicsworkbench.org/rest"
        
        # Search for studies
        search_url = f"{base_url}/study/study_id/ST/summary"
        
        with httpx.Client(timeout=30) as client:
            response = client.get(search_url)
            response.raise_for_status()
            
            # Parse tab-separated response
            lines = response.text.strip().split("\n")
            if len(lines) < 2:
                return results
            
            headers = lines[0].split("\t")
            
            for line in lines[1:max_results + 1]:
                fields = line.split("\t")
                if len(fields) < len(headers):
                    continue
                
                study_data = dict(zip(headers, fields))
                
                study_id = study_data.get("study_id", "")
                title = study_data.get("study_title", "")
                summary = study_data.get("study_summary", "")
                species = study_data.get("species", "")
                
                # Filter by terms if provided
                text_to_search = f"{title} {summary}".lower()
                matches_query = False
                
                for term_list in terms.values():
                    for term in term_list:
                        if term.lower() in text_to_search:
                            matches_query = True
                            break
                    if matches_query:
                        break
                
                # If no specific terms, include all metabolomics datasets
                if not any(terms.values()) or matches_query:
                    result = DatasetResult(
                        accession=study_id,
                        title=title,
                        description=summary,
                        source="metabolomics_workbench",
                        species=species,
                        assay_type="Metabolomics",
                        url=f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={study_id}",
                    )
                    results.append(result)
        
        logger.info(f"Found {len(results)} datasets from Metabolomics Workbench")
        
    except Exception as e:
        logger.error(f"Metabolomics Workbench search failed: {e}")
    
    return results


def deduplicate_results(results: List[DatasetResult]) -> List[DatasetResult]:
    """Remove duplicate datasets based on accession patterns and title similarity."""
    seen_accessions: Set[str] = set()
    deduplicated = []
    
    for result in results:
        # Check for known accession mapping patterns
        normalized_accession = result.accession
        
        # GSE <-> E-GEOD mapping
        if result.accession.startswith("GSE"):
            geo_id = result.accession[3:]  # Remove GSE prefix
            equivalent_ae = f"E-GEOD-{geo_id}"
            seen_accessions.add(equivalent_ae)
        elif result.accession.startswith("E-GEOD-"):
            ae_id = result.accession[7:]  # Remove E-GEOD- prefix
            equivalent_geo = f"GSE{ae_id}"
            seen_accessions.add(equivalent_geo)
        
        # Skip if we've seen this accession or its equivalent
        if normalized_accession in seen_accessions:
            continue
        
        # Check for title similarity with existing results
        is_duplicate = False
        for existing in deduplicated:
            similarity = SequenceMatcher(None, result.title.lower(), existing.title.lower()).ratio()
            if similarity > 0.9:  # Very high similarity threshold
                is_duplicate = True
                break
        
        if not is_duplicate:
            deduplicated.append(result)
            seen_accessions.add(normalized_accession)
    
    return deduplicated


def rank_results(results: List[DatasetResult], terms: Dict[str, List[str]]) -> List[DatasetResult]:
    """Rank results by relevance to search terms."""
    for result in results:
        score = 0.0
        text_to_score = f"{result.title} {result.description}".lower()
        
        # Score based on term matches
        for term_type, term_list in terms.items():
            for term in term_list:
                if term.lower() in text_to_score:
                    # Weight different term types
                    if term_type == "disease":
                        score += 3.0
                    elif term_type == "assay_type":
                        score += 2.0
                    elif term_type == "tissue":
                        score += 1.5
                    elif term_type == "species":
                        score += 1.0
        
        # Bonus for recent datasets (if we had date info)
        # Bonus for larger sample sizes
        if result.sample_count:
            score += min(result.sample_count / 100, 1.0)  # Up to 1 point for sample size
        
        result.score = score
    
    # Sort by score descending
    return sorted(results, key=lambda x: x.score, reverse=True)


def find_datasets_by_nl(
    query: str,
    repositories: Optional[List[str]] = None,
    max_results: int = 100,
) -> DatasetFinderResult:
    """
    Find datasets across repositories using natural language query.
    
    Args:
        query: Natural language search query
        repositories: List of repositories to search ("geo", "arrayexpress", "metabolomics_workbench")
        max_results: Maximum number of results to return
    
    Returns:
        DatasetFinderResult with search results and metadata
    """
    if repositories is None:
        repositories = ["geo", "arrayexpress", "metabolomics_workbench"]
    
    logger.info(f"Starting dataset search for query: '{query}'")
    
    # Extract search terms
    extracted_terms = extract_search_terms(query)
    logger.info(f"Extracted terms: {extracted_terms}")
    
    # Search each repository
    all_results = []
    sources_searched = []
    sources_failed = []
    
    per_source_limit = max(10, max_results // len(repositories))
    
    if "geo" in repositories:
        sources_searched.append("geo")
        try:
            geo_results = search_geo(extracted_terms, per_source_limit)
            all_results.extend(geo_results)
        except Exception as e:
            logger.error(f"GEO search failed: {e}")
            sources_failed.append("geo")
    
    if "arrayexpress" in repositories:
        sources_searched.append("arrayexpress")
        try:
            ae_results = search_arrayexpress(extracted_terms, per_source_limit)
            all_results.extend(ae_results)
        except Exception as e:
            logger.error(f"ArrayExpress search failed: {e}")
            sources_failed.append("arrayexpress")
    
    if "metabolomics_workbench" in repositories:
        sources_searched.append("metabolomics_workbench")
        try:
            mw_results = search_metabolomics_workbench(extracted_terms, per_source_limit)
            all_results.extend(mw_results)
        except Exception as e:
            logger.error(f"Metabolomics Workbench search failed: {e}")
            sources_failed.append("metabolomics_workbench")
    
    # Deduplicate and rank results
    deduplicated_results = deduplicate_results(all_results)
    ranked_results = rank_results(deduplicated_results, extracted_terms)
    
    # Limit final results
    final_results = ranked_results[:max_results]
    
    logger.info(f"Search complete: {len(final_results)} unique results from {len(sources_searched)} sources")
    
    return DatasetFinderResult(
        query=query,
        extracted_terms=extracted_terms,
        results=final_results,
        total_found=len(final_results),
        sources_searched=sources_searched,
        sources_failed=sources_failed,
    )


__all__ = [
    "DatasetResult",
    "DatasetFinderResult",
    "find_datasets_by_nl",
    "extract_search_terms",
]
