# Pathway Analysis ID Mapping Implementation Status

**Last Updated**: 2025-01-XX  
**Status**: ✅ Complete

## Overview

All ID mapping services for pathway analysis are fully implemented. These services map feature identifiers (genes, proteins, metabolites) to pathway database IDs (KEGG, Reactome, UniProt) for use in pathway enrichment analysis.

## Implementation Status

### ✅ Complete ID Mapping Functions

#### 1. UniProt Mapping

**Function**: `map_protein_to_uniprot(protein_id: str) -> Optional[str]`

**Location**: `amprenta_rag/analysis/id_mapping.py`

**Features**:
- ✅ Handles various input formats:
  - Gene symbols (e.g., "TP53")
  - UniProt IDs (already in correct format)
  - FASTA headers (e.g., "sp|P04637|TP53_HUMAN")
  - Ensembl protein IDs
- ✅ Uses UniProt REST API for mapping
- ✅ Caching to avoid repeated API calls
- ✅ Retry logic with exponential backoff

**API Used**: UniProt REST API (`https://rest.uniprot.org/uniprotkb/search`)

#### 2. KEGG Gene Mapping

**Function**: `map_gene_to_kegg(gene_symbol: str, organism: str = "hsa") -> Optional[str]`

**Location**: `amprenta_rag/analysis/id_mapping.py`

**Features**:
- ✅ Maps gene symbols to KEGG gene IDs (e.g., "TP53" → "hsa:7157")
- ✅ Uses MyGene.info API to get NCBI gene ID as intermediate step
- ✅ Falls back to direct KEGG search if NCBI mapping fails
- ✅ Supports organism codes (default: "hsa" for human)
- ✅ Caching and retry logic

**APIs Used**:
- MyGene.info API (`https://mygene.info/v3/query`)
- KEGG REST API (`https://rest.kegg.jp`)

#### 3. KEGG Metabolite Mapping

**Function**: `map_metabolite_to_kegg(metabolite_name: str) -> Optional[str]`

**Location**: `amprenta_rag/analysis/id_mapping.py`

**Features**:
- ✅ Maps metabolite names to KEGG compound IDs (e.g., "Glucose" → "C00031")
- ✅ Uses KEGG find API for compound search
- ✅ Case-insensitive matching
- ✅ Caching and retry logic

**API Used**: KEGG REST API (`https://rest.kegg.jp/find/compound/{name}`)

#### 4. KEGG Protein Mapping

**Function**: `map_protein_to_kegg(protein_id: str, organism: str = "hsa") -> Optional[str]`

**Location**: `amprenta_rag/analysis/id_mapping.py`

**Features**:
- ✅ Maps protein identifiers to KEGG gene IDs via UniProt
- ✅ Two-step process: protein → UniProt → KEGG
- ✅ Supports organism codes
- ✅ Caching and retry logic

**APIs Used**:
- UniProt REST API (via `map_protein_to_uniprot`)
- KEGG REST API (`https://rest.kegg.jp/conv/{organism}/uniprot:{id}`)

#### 5. Reactome Gene Mapping

**Function**: `map_gene_to_reactome(gene_symbol: str) -> Optional[str]`

**Location**: `amprenta_rag/analysis/id_mapping.py`

**Features**:
- ✅ Maps gene symbols to Reactome stable identifiers
- ✅ Uses Reactome identifier mapping API
- ✅ Caching and retry logic

**API Used**: Reactome Content Service (`https://reactome.org/ContentService/data/query/ids/map`)

#### 6. Reactome Protein Mapping

**Function**: `map_protein_to_reactome(protein_id: str) -> Optional[str]`

**Location**: `amprenta_rag/analysis/id_mapping.py`

**Features**:
- ✅ Maps protein identifiers to Reactome stable identifiers via UniProt
- ✅ Two-step process: protein → UniProt → Reactome
- ✅ Caching and retry logic

**APIs Used**:
- UniProt REST API (via `map_protein_to_uniprot`)
- Reactome Content Service (`https://reactome.org/ContentService/data/query/map/uniprot/{id}`)

### ✅ Batch Mapping Function

**Function**: `batch_map_features_to_pathway_ids(features_by_type: Dict[str, Set[str]], pathway_source: str = "KEGG", organism: str = "hsa") -> Dict[str, Dict[str, Optional[str]]]`

**Location**: `amprenta_rag/analysis/id_mapping.py`

**Features**:
- ✅ Batch mapping of multiple features by type
- ✅ Supports both KEGG and Reactome
- ✅ Rate limiting to be respectful to APIs
- ✅ Returns nested dictionary: `feature_type -> feature_name -> pathway_id`

## Integration with Pathway Analysis

All ID mapping functions are integrated into the pathway analysis module:

**Location**: `amprenta_rag/analysis/pathway_analysis.py`

- ✅ `map_features_to_kegg_pathways()` uses ID mapping functions
- ✅ `map_features_to_reactome_pathways()` uses ID mapping functions
- ✅ `perform_pathway_enrichment()` uses pathway mapping functions

## Features

### Caching

- ✅ In-memory cache for ID mappings (`_id_mapping_cache`)
- ✅ Avoids repeated API calls for the same identifiers
- ✅ Cache key includes mapping type and parameters

### Error Handling

- ✅ Graceful handling of API failures
- ✅ Returns `None` when mapping fails (non-blocking)
- ✅ Logging of errors and warnings
- ✅ Retry logic with exponential backoff

### Rate Limiting

- ✅ Respectful rate limiting for API calls
- ✅ KEGG: ~0.5 second delay between requests
- ✅ Reactome: ~0.1 second delay between requests
- ✅ Batch operations include rate limiting

### Retry Strategy

- ✅ HTTP session with retry adapter
- ✅ Retries on: 429 (rate limit), 500, 502, 503, 504
- ✅ Exponential backoff (factor=1, total=3 retries)

## API Resources

### UniProt

- **Base URL**: `https://rest.uniprot.org`
- **Search Endpoint**: `/uniprotkb/search`
- **Documentation**: https://www.uniprot.org/help/api
- **Rate Limits**: Reasonable use (no strict limits documented)

### KEGG

- **Base URL**: `https://rest.kegg.jp`
- **Endpoints Used**:
  - `/find/{db}/{query}` - Search by name
  - `/conv/{target}/{source}:{id}` - Convert IDs
  - `/link/pathway/{id}` - Get pathway links
  - `/get/{pathway_id}` - Get pathway info
- **Documentation**: https://www.kegg.jp/kegg/rest/keggapi.html
- **Rate Limits**: ~1 request/second recommended

### Reactome

- **Base URL**: `https://reactome.org/ContentService`
- **Endpoints Used**:
  - `/data/query/ids/map` - Map identifiers
  - `/data/query/map/uniprot/{id}` - Map UniProt IDs
- **Documentation**: https://reactome.org/ContentService/
- **Rate Limits**: ~10 requests/second

### MyGene.info

- **Base URL**: `https://mygene.info/v3`
- **Endpoint**: `/query`
- **Documentation**: https://docs.mygene.info/
- **Rate Limits**: Reasonable use

## Usage Examples

### Map Gene to KEGG

```python
from amprenta_rag.analysis.id_mapping import map_gene_to_kegg

kegg_id = map_gene_to_kegg("TP53", organism="hsa")
# Returns: "hsa:7157"
```

### Map Protein to UniProt

```python
from amprenta_rag.analysis.id_mapping import map_protein_to_uniprot

uniprot_id = map_protein_to_uniprot("TP53")
# Returns: "P04637"
```

### Map Metabolite to KEGG

```python
from amprenta_rag.analysis.id_mapping import map_metabolite_to_kegg

kegg_id = map_metabolite_to_kegg("Glucose")
# Returns: "C00031"
```

### Batch Mapping

```python
from amprenta_rag.analysis.id_mapping import batch_map_features_to_pathway_ids

features_by_type = {
    "gene": {"TP53", "BRCA1"},
    "protein": {"P04637", "P38398"},
    "metabolite": {"Glucose", "ATP"},
}

results = batch_map_features_to_pathway_ids(
    features_by_type,
    pathway_source="KEGG",
    organism="hsa"
)
# Returns: {
#     "gene": {"TP53": "hsa:7157", "BRCA1": "hsa:672", ...},
#     "protein": {"P04637": "hsa:7157", ...},
#     "metabolite": {"Glucose": "C00031", ...},
# }
```

## Testing

### Manual Testing

Test individual mapping functions:

```python
from amprenta_rag.analysis.id_mapping import (
    map_gene_to_kegg,
    map_protein_to_uniprot,
    map_metabolite_to_kegg,
)

# Test gene mapping
kegg_id = map_gene_to_kegg("TP53")
print(f"TP53 → KEGG: {kegg_id}")

# Test protein mapping
uniprot_id = map_protein_to_uniprot("TP53")
print(f"TP53 → UniProt: {uniprot_id}")

# Test metabolite mapping
kegg_id = map_metabolite_to_kegg("Glucose")
print(f"Glucose → KEGG: {kegg_id}")
```

### Integration Testing

Test with pathway analysis:

```python
from amprenta_rag.analysis.pathway_analysis import perform_pathway_enrichment

features = {"TP53", "BRCA1", "Glucose"}
feature_types = {"gene", "metabolite"}

results = perform_pathway_enrichment(
    input_features=features,
    input_feature_types=feature_types,
    pathway_sources=["KEGG"],
)
```

## Files

### Core Implementation

- `amprenta_rag/analysis/id_mapping.py` - All ID mapping functions
- `amprenta_rag/analysis/pathway_analysis.py` - Pathway mapping using ID functions
- `amprenta_rag/analysis/enrichment.py` - Enrichment analysis using pathway mapping

## Status: ✅ Production Ready

All ID mapping services are fully implemented and ready for use. The functions handle various input formats, include proper error handling and caching, and are integrated into the pathway analysis pipeline.

## Next Steps

### Recommended Enhancements

1. **Persistent Caching**: Add file-based or database caching for ID mappings
2. **More Organisms**: Support more KEGG organism codes beyond human
3. **Alternative APIs**: Add fallback APIs for better coverage
4. **Mapping Validation**: Add validation of mapped IDs
5. **Batch Optimization**: Optimize batch operations for better performance

