# Statistical Analysis Guide

**Last Updated**: December 7, 2025

This document describes the built-in statistical analysis tools for comparing features across datasets.

---

## Overview

The platform provides **four statistical tests** for analyzing feature data:

1. **T-Test** - Compare means of two groups (parametric)
2. **ANOVA** - Compare means of multiple groups (parametric)
3. **Mann-Whitney U** - Compare distributions of two groups (non-parametric)
4. **Pearson Correlation** - Measure linear relationship between variables

All tests include:
- **Multiple testing correction** - FDR, Bonferroni, Holm, Sidak methods
- **Interactive dashboard** - Point-and-click statistical analysis
- **Programmatic access** - Python API for custom workflows
- **CSV export** - Download results for external analysis

---

## When to Use Each Test

### T-Test (Independent Samples)

**Use When**:
- Comparing **two groups** (e.g., disease vs control, treatment vs placebo)
- Data is **normally distributed** (or n > 30 per group)
- Comparing **means** of a continuous variable

**Examples**:
- Compare average log2FC of TP53 between ALS and PD datasets
- Compare C16 ceramide levels between CSF and plasma datasets
- Compare gene expression between responders and non-responders

**Assumptions**:
- Independence of observations
- Approximate normality (relaxed with large n)
- Continuous data

**Interpretation**:
- Small p-value (<0.05): Groups have different means
- Large p-value (≥0.05): No evidence of difference

---

### ANOVA (One-Way)

**Use When**:
- Comparing **three or more groups** simultaneously
- Data is **normally distributed**
- Testing if **any group** differs from others

**Examples**:
- Compare ceramide levels across multiple diseases (ALS, PD, AD, HD)
- Compare gene expression across different tissue types
- Compare metabolite abundance across different time points

**Assumptions**:
- Independence of observations
- Approximate normality
- Homogeneity of variance (equal variance across groups)

**Interpretation**:
- Small p-value (<0.05): At least one group differs
- Large p-value (≥0.05): No evidence of differences
- **Note**: ANOVA doesn't tell you *which* groups differ (needs post-hoc tests)

---

### Mann-Whitney U Test

**Use When**:
- Comparing **two groups** but data is **non-normal**
- Data is **ordinal** or has **outliers**
- Comparing **distributions** rather than means
- Small sample sizes (n < 30 per group)

**Examples**:
- Compare feature ranks between datasets with outliers
- Compare scores on ordinal scales
- Robust alternative to t-test for skewed data

**Advantages**:
- No normality assumption
- Robust to outliers
- Works with ordinal data

**Interpretation**:
- Small p-value (<0.05): Distributions differ
- Large p-value (≥0.05): No evidence of difference
- **Note**: Tests for stochastic dominance, not just location shift

---

### Pearson Correlation

**Use When**:
- Measuring **linear relationship** between two continuous variables
- Testing if features **co-vary** across datasets
- Identifying **correlated features**

**Examples**:
- Correlate C16 ceramide with C18 ceramide across datasets
- Correlate gene expression with protein abundance
- Correlate fold-changes with p-values (volcano plot validation)

**Assumptions**:
- Linear relationship
- Bivariate normality
- Continuous data

**Interpretation**:
- **r = +1**: Perfect positive correlation
- **r = 0**: No linear correlation
- **r = -1**: Perfect negative correlation
- **|r| > 0.7**: Strong correlation
- **|r| > 0.4**: Moderate correlation
- **|r| > 0.2**: Weak correlation
- **p-value < 0.05**: Correlation is statistically significant

---

## Using the Dashboard

### Access Statistical Analysis

```bash
# Start Streamlit dashboard
cd scripts/dashboard
streamlit run app.py

# Navigate to "Statistical Analysis" page in sidebar
# URL: http://localhost:8501
```

### Dashboard Workflow

#### Step 1: Select Test Type

**Options**:
- T-test (two groups)
- ANOVA (multiple groups)
- Mann-Whitney (two groups, non-parametric)
- Correlation (linear relationship)

**Tip**: Start with Mann-Whitney if unsure about normality

#### Step 2: Choose Multiple Testing Correction

**Options**:
1. **fdr_bh** (False Discovery Rate, Benjamini-Hochberg) - Recommended default
2. **bonferroni** - Most conservative, controls family-wise error rate
3. **holm** - Less conservative than Bonferroni
4. **sidak** - Similar to Bonferroni but less conservative
5. **none** - No correction (use only for single tests)

**When to Use**:
- **FDR (fdr_bh)**: Default for most cases, good balance of power and control
- **Bonferroni**: When false positives are very costly (e.g., clinical decisions)
- **None**: Only testing one feature (no multiple testing problem)

#### Step 3: Select Datasets

**Instructions**:
- Use multiselect dropdown to choose datasets
- Default: First 2 datasets pre-selected
- Can select as many as needed (but consider grouping)

**For Two-Group Tests** (T-test, Mann-Whitney):
- First N datasets = Group A
- Remaining datasets = Group B
- Set split point with number input

**For ANOVA**:
- All selected datasets used (2+ groups)

**For Correlation**:
- Values collected across selected datasets

#### Step 4: Enter Feature Name

**Format**: Exact feature name as stored in database

**Examples**:
- `TP53` (gene)
- `P04637` (protein UniProt ID)
- `Glutamate` (metabolite)
- `Cer(d18:1/16:0)` (lipid)

**Tips**:
- Feature names are case-sensitive
- Use autocomplete if available
- Check dataset for exact spelling

#### Step 5: Set Group Split (Two-Group Tests Only)

**Purpose**: Define which datasets belong to Group A vs Group B

**Example**:
- Selected 6 datasets
- Split point = 3
- Group A = Datasets 1-3 (first 3)
- Group B = Datasets 4-6 (remaining 3)

**Use Case**:
- Group A = ALS datasets
- Group B = Control datasets
- Compare feature between disease and control

#### Step 6: Run Test

Click **"Run test"** button to execute analysis

**Processing**:
1. Extracts feature values from selected datasets
2. Groups values according to split point
3. Runs statistical test
4. Applies multiple testing correction (if selected)
5. Displays results table

### Results Table

**Columns**:
- `feature`: Feature name tested
- `statistic`: Test statistic value
- `pvalue`: Raw p-value
- `interpretation`: Significance level (*, **, ***, or not significant)
- `pvalue_adj`: Adjusted p-value (if correction applied)
- `reject`: Boolean flag (True = significant after correction)

**Significance Levels**:
- `*** significant`: p < 0.001 (highly significant)
- `** significant`: p < 0.01 (very significant)
- `* significant`: p < 0.05 (significant)
- `not significant`: p ≥ 0.05 (no evidence of effect)

### Export Results

Click **"Export CSV"** button to download results table

**File format**: CSV with all result columns

**Use cases**:
- Import into R, Python, Excel for further analysis
- Create custom visualizations
- Combine with other datasets
- Archive results

---

## Interpretation Guide

### T-Test Results

**Example**:
```
feature: TP53
statistic: 3.45
pvalue: 0.002
interpretation: ** significant
pvalue_adj: 0.018 (FDR)
reject: True
```

**Interpretation**:
- TP53 differs significantly between groups (p=0.002)
- Effect remains significant after FDR correction (adj p=0.018)
- Positive statistic suggests Group A > Group B (check means to confirm)

**Next Steps**:
- Calculate mean ± SD for each group
- Create box plot or violin plot
- Check for outliers
- Validate with independent dataset

---

### ANOVA Results

**Example**:
```
feature: Cer(d18:1/16:0)
statistic: 8.92 (F-statistic)
pvalue: 0.0003
interpretation: *** significant
pvalue_adj: 0.003
reject: True
```

**Interpretation**:
- C16 ceramide differs across groups (p=0.0003)
- At least one group is different from others
- Effect remains significant after correction

**Limitations**:
- ANOVA doesn't specify *which* groups differ
- Need post-hoc tests (Tukey HSD, etc.) to identify pairs
- Significant ANOVA doesn't mean all groups differ

**Next Steps**:
- Run pairwise t-tests with correction
- Create box plot showing all groups
- Calculate effect sizes (eta-squared, omega-squared)

---

### Mann-Whitney Results

**Example**:
```
feature: APOE
statistic: 45.0 (U statistic)
pvalue: 0.023
interpretation: * significant
pvalue_adj: 0.092
reject: False
```

**Interpretation**:
- APOE distributions differ nominally (p=0.023)
- Effect NOT significant after FDR correction (adj p=0.092)
- Evidence of difference is weak (marginal significance)

**Key Point**:
- Mann-Whitney tests distributions, not just means
- Groups can differ in spread, median, or shape
- Non-significant after correction suggests possible false positive

**Next Steps**:
- Visualize distributions with histograms
- Check sample sizes (may lack power)
- Consider biological significance vs statistical significance

---

### Correlation Results

**Example**:
```
feature: Correlation analysis
statistic: 0.78 (Pearson r)
pvalue: 0.0001
interpretation: *** significant
```

**Interpretation**:
- Strong positive correlation (r=0.78)
- Highly significant (p=0.0001)
- 78% of variance in Y explained by linear relationship with X

**Correlation Strength**:
- |r| > 0.7: Strong
- |r| > 0.4: Moderate
- |r| > 0.2: Weak
- |r| ≤ 0.2: Negligible

**Cautions**:
- Correlation ≠ causation
- May be driven by outliers (check scatterplot)
- Only measures linear relationships

**Next Steps**:
- Create scatterplot with regression line
- Check for outliers (Cook's distance)
- Consider non-linear relationships (Spearman correlation)

---

## Multiple Testing Correction Explained

### The Problem

When testing multiple features, the probability of false positives increases:

**Example**:
- Test 100 features at α=0.05
- Expected false positives = 100 × 0.05 = **5 false discoveries**
- Without correction, ~5 features will appear significant by chance

### False Discovery Rate (FDR)

**Method**: Benjamini-Hochberg procedure

**How It Works**:
1. Rank p-values from smallest to largest
2. Find largest i where p(i) ≤ (i/m) × α
3. Reject all H0 with p ≤ p(i)

**Advantages**:
- Controls proportion of false discoveries (not family-wise error)
- More powerful than Bonferroni
- Appropriate for exploratory research

**When to Use**: Default for most genomics/omics research

**Example**:
```
Raw p-values: [0.001, 0.02, 0.04, 0.06, 0.10]
FDR adjusted:  [0.005, 0.05, 0.067, 0.075, 0.10]
Reject at α=0.05: First 2 features
```

---

### Bonferroni Correction

**Method**: Multiply each p-value by number of tests (or divide α by number of tests)

**Formula**: `p_adjusted = p_raw × m` (m = number of tests)

**Advantages**:
- Controls family-wise error rate (FWER)
- Simple to understand
- Guarantees no false positives (at chosen α level)

**Disadvantages**:
- Very conservative (low power)
- May miss true positives

**When to Use**: When false positives are very costly (e.g., clinical trials, confirmatory studies)

**Example**:
```
Testing 20 features, α=0.05
Bonferroni threshold: 0.05/20 = 0.0025
Only features with p<0.0025 considered significant
```

---

### Holm Correction

**Method**: Step-down procedure, less conservative than Bonferroni

**How It Works**:
1. Rank p-values
2. Compare smallest p to α/m
3. Compare next p to α/(m-1), etc.
4. Stop at first non-significant test

**Advantages**:
- More powerful than Bonferroni
- Still controls FWER
- Easy to understand

**When to Use**: When FWER control needed but Bonferroni too conservative

---

### When to Use No Correction

**Appropriate Cases**:
- Testing only **one feature** (no multiple testing)
- Purely **exploratory analysis** (hypothesis generation)
- Corrections applied **downstream** (in publication)

**Inappropriate Cases**:
- Testing many features and reporting all p-values
- Making strong claims based on uncorrected p-values

---

## Programmatic Access

### Run T-Test

```python
from amprenta_rag.analysis.statistical_tests import ttest_independent

# Example: Compare TP53 between two groups
group_a = [2.5, 3.1, 2.8, 3.4, 2.9]  # log2FC values from ALS datasets
group_b = [1.2, 1.5, 1.3, 1.8, 1.4]  # log2FC values from control datasets

result = ttest_independent(group_a, group_b)
print(result)
# Output: {'statistic': 10.45, 'pvalue': 0.0001, 'interpretation': '*** significant'}
```

---

### Run ANOVA

```python
from amprenta_rag.analysis.statistical_tests import anova_oneway

# Example: Compare across three diseases
als_values = [3.2, 3.5, 3.1, 3.4]
pd_values = [2.1, 2.3, 2.0, 2.4]
ad_values = [1.5, 1.7, 1.6, 1.8]

result = anova_oneway(als_values, pd_values, ad_values)
print(result)
# Output: {'statistic': 28.4, 'pvalue': 0.00001, 'interpretation': '*** significant'}
```

---

### Run Mann-Whitney

```python
from amprenta_rag.analysis.statistical_tests import mann_whitney

# Example: Non-parametric comparison with outliers
group_a = [2.1, 2.3, 2.2, 15.0]  # Note outlier (15.0)
group_b = [1.1, 1.2, 1.3, 1.4]

result = mann_whitney(group_a, group_b)
print(result)
# Output: {'statistic': 16.0, 'pvalue': 0.014, 'interpretation': '* significant'}
```

---

### Calculate Correlation

```python
from amprenta_rag.analysis.statistical_tests import pearson_corr

# Example: Correlate two features across datasets
c16_ceramide = [2.5, 3.0, 2.8, 3.2, 2.9]
c18_ceramide = [2.0, 2.4, 2.2, 2.6, 2.3]

result = pearson_corr(c16_ceramide, c18_ceramide)
print(result)
# Output: {'statistic': 0.96, 'pvalue': 0.01, 'interpretation': '** significant'}
```

---

### Apply Multiple Testing Correction

```python
from amprenta_rag.analysis.statistical_tests import adjust_pvalues

# Example: Correct 10 p-values using FDR
raw_pvalues = [0.001, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.25, 0.40, 0.60]

adjusted_pvals, reject_flags = adjust_pvalues(raw_pvalues, method="fdr_bh")
print("Adjusted p-values:", adjusted_pvals)
print("Reject null:", reject_flags)

# Output:
# Adjusted p-values: [0.010, 0.100, 0.100, 0.125, 0.160, ...]
# Reject null: [True, True, True, False, False, ...]
```

**Available Methods**:
- `"fdr_bh"` - Benjamini-Hochberg FDR
- `"bonferroni"` - Bonferroni correction
- `"holm"` - Holm step-down
- `"sidak"` - Sidak correction

---

## Real-World Examples

### Example 1: Compare Ceramide Levels (Disease vs Control)

**Scenario**: Test if C16 ceramide differs between ALS and control CSF

**Data**:
- ALS datasets (n=15): Mean log2FC = 2.8
- Control datasets (n=15): Mean log2FC = 0.5

**Workflow**:
1. Dashboard → Statistical Analysis
2. Select "T-test"
3. Choose FDR correction
4. Select 30 datasets (15 ALS + 15 control)
5. Set split point = 15
6. Enter feature: "Cer(d18:1/16:0)"
7. Run test

**Result**:
```
statistic: 8.45
pvalue: 0.000001
interpretation: *** significant
pvalue_adj: 0.000005
reject: True
```

**Conclusion**: C16 ceramide is significantly elevated in ALS CSF (p<0.000001, FDR-corrected)

---

### Example 2: Multi-Disease Comparison (ANOVA)

**Scenario**: Compare TP53 expression across 4 neurodegenerative diseases

**Data**:
- ALS (n=10): Mean = 3.2
- PD (n=10): Mean = 2.1
- AD (n=10): Mean = 2.8
- HD (n=10): Mean = 1.5

**Workflow**:
1. Dashboard → Statistical Analysis
2. Select "ANOVA"
3. Choose FDR correction
4. Select 40 datasets (10 per disease)
5. Enter feature: "TP53"
6. Run test

**Result**:
```
statistic: 15.6 (F-statistic)
pvalue: 0.00001
interpretation: *** significant
```

**Conclusion**: TP53 expression differs across diseases (p<0.00001). Need post-hoc tests to identify which pairs differ.

---

### Example 3: Feature Correlation

**Scenario**: Test if C16 and C18 ceramides correlate across datasets

**Data**: 20 datasets with both ceramide measurements

**Workflow**:
1. Dashboard → Statistical Analysis
2. Select "Correlation"
3. Select 20 datasets
4. Enter feature: "Cer(d18:1/16:0)" (or use programmatic approach for two features)
5. Run test

**Result** (programmatic):
```python
c16_values = [...]  # Extract from 20 datasets
c18_values = [...]  # Extract from same 20 datasets
result = pearson_corr(c16_values, c18_values)
# Output: {'statistic': 0.85, 'pvalue': 0.000001, 'interpretation': '*** significant'}
```

**Conclusion**: Strong positive correlation (r=0.85, p<0.000001). C16 and C18 ceramides tend to increase/decrease together.

---

## Best Practices

### 1. Check Assumptions

**Before t-test/ANOVA**:
- Visualize data with box plots or histograms
- Check for approximate normality (Q-Q plots)
- Check for equal variance (Levene's test)
- Use Mann-Whitney if assumptions violated

### 2. Consider Sample Size

**Rules of thumb**:
- T-test: n ≥ 30 per group (central limit theorem applies)
- ANOVA: n ≥ 20 per group
- Mann-Whitney: Works with small samples (n ≥ 5)
- Correlation: n ≥ 20 pairs recommended

**Low power warning**:
- Small samples may miss real effects (Type II error)
- Non-significant ≠ no difference (may lack power)

### 3. Report Effect Sizes

**P-values alone insufficient**:
- Report means ± SD for each group
- Report effect sizes (Cohen's d, eta-squared)
- Consider biological significance vs statistical significance

**Example**:
```
TP53 expression: ALS (3.2 ± 0.5) vs Control (3.1 ± 0.4)
T-test: p=0.04 (significant)
Cohen's d: 0.22 (small effect)
Conclusion: Statistically significant but small biological effect
```

### 4. Use Appropriate Correction

**Guidelines**:
- Testing 1-5 features: FDR or no correction
- Testing 10-100 features: FDR
- Testing 100-10,000 features: FDR
- Testing >10,000 features: Consider q-values (more sophisticated)
- Confirmatory study: Bonferroni

### 5. Validate Findings

**Cross-validation**:
- Split data into discovery and validation sets
- Test in independent cohort
- Replicate with different omics type
- Confirm with targeted experiments

---

## Troubleshooting

### Issue: "Insufficient Data"

**Problem**: Test returns NaN with "insufficient data" interpretation

**Causes**:
- Too few values in one or both groups
- Feature not present in selected datasets
- All values are NaN (missing data)

**Solutions**:
1. Check feature name spelling
2. Select more datasets
3. Verify feature exists in datasets
4. Check for missing values

**Diagnostic**:
```python
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset

db = next(get_db())
dataset = db.query(Dataset).filter(Dataset.name == "your-dataset").first()
feature = next((f for f in dataset.features if f.name == "TP53"), None)
if feature:
    print(f"Found feature, external_ids: {feature.external_ids}")
else:
    print("Feature not found in dataset")
```

---

### Issue: Non-Significant After Correction

**Problem**: p<0.05 becomes non-significant after FDR correction

**Explanation**: This is expected behavior—correction protects against false positives

**Options**:
1. **Accept result**: Effect may be false positive
2. **Increase sample size**: More power to detect real effects
3. **Use less stringent α**: α=0.10 instead of 0.05 (exploratory only)
4. **Report both**: Raw and adjusted p-values for transparency

---

### Issue: All Tests Non-Significant

**Problem**: No features show significant differences

**Possible Causes**:
1. **No real biological difference**: Groups may be similar
2. **Low power**: Sample sizes too small
3. **High variability**: Large within-group variance obscures effects
4. **Wrong test**: May need different statistical approach

**Solutions**:
- Increase sample sizes
- Use more sensitive test (Mann-Whitney instead of t-test)
- Check data quality (outliers, measurement errors)
- Consider effect size regardless of p-value

---

## See Also

- [Visualization Guide](VISUALIZATIONS.md) - Volcano plots, heatmaps for statistical results
- [Quality Checks Guide](QUALITY_CHECKS.md) - Ensure data quality before analysis
- [Usage Examples](USAGE_EXAMPLES.md) - Data preparation workflows
- [API Reference](API_REFERENCE.md) - Statistical functions documentation

---

**Questions or Issues?**

See [Troubleshooting Guide](TROUBLESHOOTING.md) or file an issue on GitHub.

