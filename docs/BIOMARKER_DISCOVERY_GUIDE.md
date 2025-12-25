# Biomarker Discovery Guide

## Overview

The Biomarker Discovery module enables identification of differentially expressed features (genes, proteins, metabolites, lipids) that distinguish between experimental conditions or disease states. It combines multiple statistical and machine learning methods to produce robust, consensus-ranked candidate biomarkers.

**Use Cases:**
- Identify disease-associated lipids in Alzheimer's vs. control samples
- Find protein markers that differentiate drug responders from non-responders
- Discover gene signatures for patient stratification
- Validate existing signatures with independent datasets

---

## Methods

The module implements three complementary approaches for feature selection:

### 1. Statistical Tests

Classical hypothesis testing with multiple testing correction:

- **t-test** - Parametric comparison for normally distributed continuous data
- **Mann-Whitney U** - Non-parametric alternative for skewed distributions
- **ANOVA** - Multi-group comparisons (>2 conditions)
- **FDR correction** - Benjamini-Hochberg procedure to control false discovery rate

**Advantages:** Fast, interpretable, widely accepted in publications

**Limitations:** Univariate (doesn't capture interactions), sensitive to outliers

### 2. Stability Selection

Bootstrap-based feature selection with regularization:

- **Bootstrap resampling** - 100 iterations with 80% sample size
- **LassoCV** - L1-regularized logistic regression with cross-validation
- **Selection frequency** - Features selected in >50% of bootstrap iterations
- **Stability scores** - Proportion of iterations where feature was selected

**Advantages:** Handles collinearity, robust to noise, multivariate

**Limitations:** Requires sufficient samples (n > 50 recommended), computationally intensive

### 3. Cross-Validated Feature Importance

Ensemble-based importance ranking:

- **RandomForest** classifier with 100 trees
- **5-fold stratified cross-validation** - Averaged importance across folds
- **Mean decrease in impurity** - Gini importance metric
- **Normalized scores** - 0-100 scale for interpretability

**Advantages:** Captures non-linear relationships, robust to irrelevant features

**Limitations:** Biased toward high-cardinality features, requires balanced classes

---

## Consensus Ranking

The module combines results from all three methods using **average rank** across methods:

1. Each method produces a ranked list of features (best to worst)
2. For each feature, calculate rank from each method (1 = best)
3. Compute average rank: `(rank_ttest + rank_stability + rank_importance) / 3`
4. Sort features by average rank (lower = more consistent)

**Why consensus ranking?**
- Reduces method-specific biases
- Increases robustness to distributional assumptions
- Highlights features that perform well across diverse criteria

**Interpretation:**
- **Average rank < 10**: High-confidence biomarkers (top-ranked by all methods)
- **Average rank 10-50**: Strong candidates
- **Average rank > 50**: Weak or method-specific signals

---

## Dashboard Walkthrough

The Biomarker Discovery page has three tabs:

### Setup Tab

1. **Select Dataset** - Choose dataset with feature matrix
2. **Select Outcome Column** - Binary or multi-class grouping variable
3. **Configure Methods** - Enable/disable statistical tests, stability selection, importance
4. **Set Thresholds** - FDR cutoff (default: 0.05), stability frequency (default: 0.5)
5. **Run Discovery** - Executes analysis and generates results

### Methods Tab

Displays method-specific results in expandable sections:

- **Statistical Tests**: Feature, p-value, FDR, fold-change
- **Stability Selection**: Feature, selection frequency, stability score
- **CV Importance**: Feature, mean importance, std deviation

Each table is sortable and filterable for exploration.

### Results Tab

**Consensus Ranked Table:**
- Top 50 biomarkers sorted by average rank
- Columns: Feature name, average rank, ranks from each method
- Download CSV button for export

**Venn Diagram** (if >1 method enabled):
- Shows overlap between top features from each method
- Helps identify method-specific vs. consensus biomarkers

---

## API Reference

### GET /api/biomarker/methods

List available biomarker discovery methods.

**Response:**
```json
{
  "methods": [
    {"id": "ttest", "name": "T-Test", "type": "statistical"},
    {"id": "mannwhitney", "name": "Mann-Whitney U", "type": "statistical"},
    {"id": "anova", "name": "ANOVA", "type": "statistical"},
    {"id": "stability", "name": "Stability Selection", "type": "ml"},
    {"id": "importance", "name": "CV Feature Importance", "type": "ml"}
  ]
}
```

### POST /api/biomarker/discover

Run biomarker discovery analysis.

**Request:**
```json
{
  "dataset_id": "uuid-string",
  "outcome_column": "disease_status",
  "methods": ["ttest", "stability", "importance"],
  "fdr_threshold": 0.05,
  "stability_threshold": 0.5
}
```

**Response:**
```json
{
  "consensus_biomarkers": [
    {
      "feature_name": "Cer(d18:1/16:0)",
      "average_rank": 2.3,
      "ttest_rank": 1,
      "stability_rank": 3,
      "importance_rank": 3
    }
  ],
  "method_results": {
    "ttest": [...],
    "stability": [...],
    "importance": [...]
  }
}
```

---

## Best Practices

### Sample Size Requirements

- **Minimum**: 30 samples (15 per group for binary)
- **Recommended**: 50+ samples for stability selection
- **Optimal**: 100+ samples for robust feature importance

Underpowered analyses produce unstable results. Use **power analysis** to determine adequate sample size before running discovery.

### FDR Threshold Selection

- **Exploratory**: FDR < 0.10 (generate hypotheses)
- **Standard**: FDR < 0.05 (typical publication threshold)
- **Conservative**: FDR < 0.01 (high-confidence only)

Lower thresholds reduce false positives but may miss true biomarkers.

### Method Selection

**For small samples (n < 50):**
- Use Mann-Whitney U or t-test only
- Skip stability selection (unreliable with small n)

**For large samples (n > 100):**
- Use all three methods for consensus ranking
- Stability selection becomes more reliable

**For multi-group comparisons:**
- Use ANOVA instead of t-test
- Stability selection works for multi-class

### Interpretation Recommendations

1. **Validate top biomarkers** in independent cohorts
2. **Check biological plausibility** - Are top features known disease markers?
3. **Visualize distributions** - Box plots/violin plots for top features
4. **Pathway enrichment** - Do consensus biomarkers cluster in pathways?
5. **Literature search** - Cross-reference with existing biomarker literature

---

## Related Documentation

- [Advanced Analytics Features](ROADMAP.md#advanced-analytics-features) - See MOA Inference, Signature Scoring
- [Signature Management](USER_GUIDE.md) - Create signatures from discovered biomarkers
- [API Reference](API_REFERENCE.md) - Full API documentation

