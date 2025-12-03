# Top 3 Feature Enhancements - Detailed Specifications

**Status**: Ready for implementation

This document provides detailed specifications for the top 3 highest-value feature enhancements, with full implementation details and Notion agent instructions where needed.

---

## 🥇 #1: Enhanced Dataset Feature Extraction & Caching

### **Strategic Value**

**Current State**: 
- Every signature scoring operation queries Notion APIs to extract dataset features
- For a dataset with 100 features, this means 100+ Notion API calls
- Scoring 10 datasets against 5 signatures = 50 scoring operations × multiple API calls = hundreds of API calls
- **Bottleneck**: Network latency and API rate limits

**Impact**: 
- **10-100x performance improvement** for batch scoring operations
- Reduces Notion API load significantly
- Enables real-time signature scoring

**Why This Matters**:
- Signature scoring is a core workflow
- Users will want to score datasets frequently
- Batch operations become feasible

---

### **Technical Implementation**

#### **1.1 Create Feature Cache Module**

**File**: `amprenta_rag/ingestion/dataset_feature_cache.py`

```python
from typing import Dict, Set, Optional
from datetime import datetime, timedelta
from functools import lru_cache
import logging

logger = logging.getLogger(__name__)

class DatasetFeatureCache:
    """
    In-memory cache for dataset feature sets to avoid repeated Notion API calls.
    
    Cache structure:
    {
        dataset_page_id: {
            "features_by_type": {
                "gene": set(...),
                "protein": set(...),
                "metabolite": set(...),
                "lipid": set(...)
            },
            "directions": {
                "gene": {"TP53": "↑", ...},
                ...
            },
            "timestamp": datetime,
            "omics_type": str
        }
    }
    """
    
    def __init__(self, ttl_seconds: int = 3600):
        """
        Args:
            ttl_seconds: Time-to-live for cache entries (default: 1 hour)
        """
        self.cache: Dict[str, Dict] = {}
        self.ttl = timedelta(seconds=ttl_seconds)
    
    def get_features(
        self, 
        dataset_page_id: str,
        force_refresh: bool = False
    ) -> Optional[Dict[str, Set[str]]]:
        """Get cached features, or None if not cached/stale."""
        if force_refresh or dataset_page_id not in self.cache:
            return None
        
        entry = self.cache[dataset_page_id]
        age = datetime.now() - entry["timestamp"]
        
        if age > self.ttl:
            logger.debug(
                "[CACHE] Cache entry for dataset %s is stale (age: %s), will refresh",
                dataset_page_id,
                age
            )
            return None
        
        return entry["features_by_type"]
    
    def set_features(
        self,
        dataset_page_id: str,
        features_by_type: Dict[str, Set[str]],
        directions: Optional[Dict[str, Dict[str, str]]] = None,
        omics_type: Optional[str] = None
    ):
        """Store features in cache."""
        self.cache[dataset_page_id] = {
            "features_by_type": features_by_type,
            "directions": directions or {},
            "timestamp": datetime.now(),
            "omics_type": omics_type
        }
        logger.debug(
            "[CACHE] Cached features for dataset %s (%d feature types)",
            dataset_page_id,
            len([k for k, v in features_by_type.items() if v])
        )
    
    def clear(self, dataset_page_id: Optional[str] = None):
        """Clear cache entry(ies)."""
        if dataset_page_id:
            if dataset_page_id in self.cache:
                del self.cache[dataset_page_id]
                logger.debug("[CACHE] Cleared cache for dataset %s", dataset_page_id)
        else:
            self.cache.clear()
            logger.debug("[CACHE] Cleared entire cache")
    
    def get_stats(self) -> Dict[str, int]:
        """Get cache statistics."""
        return {
            "cached_datasets": len(self.cache),
            "total_entries": sum(
                len(features) 
                for entry in self.cache.values() 
                for features in entry["features_by_type"].values()
            )
        }
```

#### **1.2 Enhance Feature Extraction to Support Files**

**Enhancement to**: `amprenta_rag/ingestion/multi_omics_scoring.py`

Add a new function:

```python
def extract_dataset_features_from_file(
    dataset_page_id: str,
    file_path: Optional[str] = None,
) -> Dict[str, Set[str]]:
    """
    Extract features directly from attached dataset files (CSV/TSV).
    
    This is faster and more accurate than querying Notion relations,
    especially for large datasets.
    
    Args:
        dataset_page_id: Notion page ID of dataset
        file_path: Optional path to dataset file (if not provided, fetch from Notion)
    
    Returns:
        Dictionary mapping feature_type → set of feature names
    """
    # 1. If file_path not provided, fetch file URL from Notion page
    # 2. Download file (CSV/TSV)
    # 3. Use existing normalization functions based on omics type
    # 4. Extract feature sets
    # 5. Return features_by_type dict
    pass
```

#### **1.3 Integrate Cache into Scoring Pipeline**

**Enhancement to**: `amprenta_rag/ingestion/signature_matching.py`

```python
# Add at module level
_feature_cache = DatasetFeatureCache(ttl_seconds=3600)

def find_matching_signatures_for_dataset(
    ...,
    use_cache: bool = True,
    force_refresh: bool = False,
) -> List[SignatureMatchResult]:
    """Enhanced with caching support."""
    
    # Check cache first
    if use_cache and not force_refresh:
        cached_features = _feature_cache.get_features(dataset_page_id)
        if cached_features:
            logger.debug(
                "[INGEST][SIGNATURE-MATCH] Using cached features for dataset %s",
                dataset_page_id
            )
            dataset_features_by_type = cached_features
        else:
            # Extract from Notion/file
            dataset_features_by_type = extract_dataset_features_by_type(...)
            # Cache for next time
            _feature_cache.set_features(
                dataset_page_id,
                dataset_features_by_type,
                directions=dataset_directions,
                omics_type=omics_type
            )
    else:
        # Fresh extraction
        dataset_features_by_type = extract_dataset_features_by_type(...)
        if use_cache:
            _feature_cache.set_features(...)
    
    # Continue with scoring...
```

#### **1.4 Add Batch Scoring Support**

**New function**: `amprenta_rag/ingestion/signature_matching.py`

```python
def score_datasets_against_signatures_batch(
    dataset_page_ids: List[str],
    signature_page_ids: Optional[List[str]] = None,
    overlap_threshold: float = 0.3,
    use_cache: bool = True,
) -> Dict[str, List[SignatureMatchResult]]:
    """
    Score multiple datasets against signatures efficiently using caching.
    
    Args:
        dataset_page_ids: List of dataset page IDs to score
        signature_page_ids: Optional list of signature page IDs (if None, uses all)
        overlap_threshold: Minimum overlap for matches
        use_cache: Whether to use feature cache
    
    Returns:
        Dictionary mapping dataset_page_id → list of matches
    """
    # 1. Pre-load all dataset features into cache
    # 2. Load all signatures once
    # 3. Score in batch
    # 4. Return results
    pass
```

---

### **Implementation Steps**

1. ✅ Create `dataset_feature_cache.py` module
2. ✅ Add caching to `extract_dataset_features_by_type`
3. ✅ Integrate cache into `find_matching_signatures_for_dataset`
4. ✅ Add batch scoring function
5. ✅ Add CLI option for batch scoring
6. ✅ Test with multiple datasets

### **Testing Plan**

```python
# Test cache hit/miss
# Test cache expiration
# Test batch scoring performance
# Compare before/after performance metrics
```

### **Notion Agent Instructions**

**None needed** - This is a pure code enhancement with no schema changes required.

---

## 🥈 #2: Automated Signature Discovery from Datasets

### **Strategic Value**

**Current State**:
- Signatures must be manually defined and curated
- Scientific discovery relies on manual pattern recognition
- Building a comprehensive signature library is slow

**Impact**:
- **Accelerate signature library growth** by 10x
- Discover patterns humans might miss
- Enable data-driven signature generation
- Support hypothesis generation

**Why This Matters**:
- Signatures are core to the platform's value
- More signatures = better dataset matching
- Enables discovery of novel biological patterns

---

### **Technical Implementation**

#### **2.1 Pattern Detection Engine**

**File**: `amprenta_rag/signatures/signature_discovery.py`

```python
def discover_signatures_from_datasets(
    dataset_page_ids: List[str],
    min_feature_frequency: float = 0.3,  # Must appear in 30% of datasets
    min_features_per_signature: int = 3,
    max_features_per_signature: int = 50,
    direction_consistency_threshold: float = 0.7,  # 70% consistent direction
    disease_filter: Optional[str] = None,
    omics_type_filter: Optional[str] = None,
) -> List[DiscoveredSignature]:
    """
    Discover recurring feature patterns across datasets.
    
    Algorithm:
    1. Extract all features from all datasets
    2. Count feature co-occurrence across datasets
    3. Cluster co-occurring features
    4. Check direction consistency
    5. Filter by frequency thresholds
    6. Generate candidate signatures
    
    Returns:
        List of DiscoveredSignature objects with:
        - features: List of (feature_name, feature_type, direction)
        - frequency: How many datasets contain this pattern
        - direction_consistency: Fraction of datasets with consistent direction
        - disease_contexts: Set of diseases where pattern appears
        - confidence_score: Statistical confidence
    """
    pass
```

#### **2.2 Co-Occurrence Analysis**

```python
def compute_feature_cooccurrence(
    dataset_features: Dict[str, Dict[str, Set[str]]],
) -> Dict[Tuple[str, str], float]:
    """
    Compute pairwise feature co-occurrence scores.
    
    Returns:
        Dictionary mapping (feature1, feature2) → co-occurrence score (0-1)
    """
    # Use Jaccard similarity or frequency-based scoring
    pass
```

#### **2.3 Candidate Signature Generation**

```python
def generate_candidate_signature_file(
    discovered_signature: DiscoveredSignature,
    output_path: str,
    signature_name: Optional[str] = None,
) -> str:
    """
    Generate a TSV file for a discovered signature.
    
    Format matches existing signature TSV format:
    feature_name    feature_type    direction    weight
    
    Returns:
        Path to generated file
    """
    pass
```

#### **2.4 Statistical Validation**

```python
def validate_discovered_signature(
    discovered_signature: DiscoveredSignature,
    all_datasets: List[str],
) -> ValidationMetrics:
    """
    Compute validation metrics for a discovered signature.
    
    Metrics:
    - Coverage: Fraction of datasets that match
    - Specificity: Fraction of non-matching datasets that don't match
    - Reproducibility: Consistency across datasets
    - Statistical significance: p-value
    """
    pass
```

---

### **Implementation Steps**

1. ✅ Create `signature_discovery.py` module
2. ✅ Implement pattern detection algorithms
3. ✅ Add co-occurrence analysis
4. ✅ Generate candidate signature files
5. ✅ Add validation metrics
6. ✅ Create CLI script: `scripts/discover_signatures.py`
7. ✅ Test with real datasets

### **CLI Script**

**File**: `scripts/discover_signatures.py`

```bash
python scripts/discover_signatures.py \
    --datasets dataset1_id dataset2_id dataset3_id \
    --min-frequency 0.3 \
    --output-dir ./discovered_signatures/ \
    --disease-filter "ALS" \
    --review
```

### **Notion Agent Instructions**

**None needed** - This generates candidate signature files for manual review before ingestion.

---

## 🥉 #3: Cross-Omics Pathway Analysis

### **Strategic Value**

**Current State**:
- System tracks individual features (genes, proteins, metabolites, lipids)
- No pathway-level context or enrichment
- Difficult to interpret multi-omics results biologically

**Impact**:
- **Biological interpretability** - pathways are more meaningful than individual features
- Enable pathway-level cross-omics reasoning
- Better scientific insights
- Connect features to known biological processes

**Why This Matters**:
- Scientists think in pathways, not individual features
- Pathways provide biological context
- Enables pathway-based signature scoring

---

### **Technical Implementation**

#### **3.1 Pathway Database Integration**

**File**: `amprenta_rag/analysis/pathway_mapping.py`

```python
def map_feature_to_pathways(
    feature_name: str,
    feature_type: str,
) -> List[PathwayMapping]:
    """
    Map a feature (gene/protein/metabolite) to known pathways.
    
    Uses:
    - KEGG API for metabolic pathways
    - Reactome API for signaling pathways
    - Local mappings for common features
    
    Returns:
        List of PathwayMapping objects:
        - pathway_id: str (e.g., "hsa00010" for KEGG)
        - pathway_name: str
        - pathway_type: str (metabolic, signaling, etc.)
        - confidence: float
    """
    pass
```

#### **3.2 Pathway Enrichment Analysis**

**File**: `amprenta_rag/analysis/pathway_enrichment.py`

```python
def compute_pathway_enrichment(
    dataset_features: Dict[str, Set[str]],
    background_features: Optional[Dict[str, Set[str]]] = None,
) -> List[EnrichedPathway]:
    """
    Compute pathway enrichment for a dataset.
    
    Uses hypergeometric test or Fisher's exact test.
    
    Returns:
        List of EnrichedPathway objects with:
        - pathway_id: str
        - pathway_name: str
        - p_value: float
        - enrichment_ratio: float
        - features_in_pathway: List[str]
    """
    pass
```

#### **3.3 Pathway-Based Cross-Omics Reasoning**

**Enhancement to**: `amprenta_rag/query/cross_omics_reasoning.py`

```python
def cross_omics_pathway_summary(
    program_page_id: str,
    top_k_pathways: int = 10,
) -> str:
    """
    Generate cross-omics summary focused on pathway-level insights.
    
    Algorithm:
    1. Extract all features from program's datasets
    2. Map features to pathways
    3. Identify enriched pathways across omics
    4. Find cross-omics pathway convergence
    5. Generate LLM summary focused on pathways
    """
    pass
```

---

### **Notion Agent Instructions**

**REQUIRED** - This enhancement requires Notion schema changes.

#### **Step 1: Create Pathways Database**

**Instructions for Notion Agent**:

```
Create a new database in Notion called "Pathways" with the following properties:

Properties:
- Name (Title) - Pathway name (e.g., "Glycolysis / Gluconeogenesis")
- Pathway ID (Text) - External pathway identifier (e.g., "hsa00010" for KEGG)
- Pathway Type (Select) - Options: "Metabolic", "Signaling", "Regulatory", "Other"
- Source (Select) - Options: "KEGG", "Reactome", "WikiPathways", "Custom"
- Description (Text) - Pathway description
- Related Genes (Relation) → Gene Features database
- Related Proteins (Relation) → Protein Features database  
- Related Metabolites (Relation) → Metabolite Features database
- Related Datasets (Relation) → Experimental Data Assets database
- Related Programs (Relation) → Programs database
```

#### **Step 2: Add Pathway Relations to Feature Databases**

**Instructions for Notion Agent**:

```
For each feature database (Gene Features, Protein Features, Metabolite Features):

Add property:
- Related Pathways (Relation) → Pathways database

This allows bidirectional linking between features and pathways.
```

#### **Step 3: Add Pathway Properties to Datasets**

**Instructions for Notion Agent**:

```
For Experimental Data Assets database:

Add properties:
- Enriched Pathways (Relation) → Pathways database
- Pathway Enrichment Summary (Rich Text) - Auto-populated summary of enriched pathways
```

#### **Step 4: Add Pathway Properties to Signatures**

**Instructions for Notion Agent**:

```
For Lipid Signatures database:

Add properties:
- Related Pathways (Relation) → Pathways database
- Pathway Context (Rich Text) - Pathway-level description of signature
```

---

### **Implementation Steps**

1. ✅ Create `pathway_mapping.py` module
2. ✅ Integrate KEGG API client
3. ✅ Integrate Reactome API client (optional)
4. ✅ Create pathway enrichment analysis module
5. ✅ Enhance cross-omics reasoning with pathway context
6. ✅ **Wait for Notion Agent** to create Pathways database
7. ✅ Implement pathway ingestion pipeline
8. ✅ Link pathways to features during ingestion
9. ✅ Test pathway-based queries

### **Dependencies**

- KEGG API access (free tier available)
- Reactome API (optional, for signaling pathways)
- Statistical libraries (scipy for enrichment tests)

---

## 📋 Summary Comparison

| Enhancement | Impact | Effort | Notion Agent | Dependencies |
|------------|--------|--------|--------------|--------------|
| **#1: Feature Caching** | ⭐⭐⭐⭐⭐ | Medium | None | None |
| **#2: Signature Discovery** | ⭐⭐⭐⭐⭐ | High | None | Statistical libs |
| **#3: Pathway Analysis** | ⭐⭐⭐⭐⭐ | High | **Yes** | Pathway APIs |

---

## 🎯 Recommendation Order

1. **Start with #1 (Feature Caching)** - Immediate performance win, no blockers
2. **Then #2 (Signature Discovery)** - Accelerates discovery, pure code
3. **Finally #3 (Pathway Analysis)** - Requires Notion schema changes, coordinate with agent

---

**Ready to implement any of these? I can start with #1 (Feature Caching) as it provides immediate value with no blockers!**

