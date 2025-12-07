# Signature Discovery v1 - Validation Report

**Generated:** 2025-12-05T16:13:26.178538

## Executive Summary

- **Unit Tests:** ✅ PASSED
- **Datasets Analyzed:** 12
- **Discovery Tests Run:** 3

## Test Results

### 1. Unit Tests

✅ All unit tests passed

### 2. Dataset Coverage Analysis

#### LIPIDOMICS
- Total datasets: 3
- Datasets with features: 1
- Sample datasets:
  - Lipid Alterations in ASAH1-Deficient Cells: Insights into Ce: 740 features

#### METABOLOMICS
- Total datasets: 3
- Datasets with features: 0

#### PROTEOMICS
- Total datasets: 3
- Datasets with features: 1
- Sample datasets:
  - Validation of Solvent Proteome Profiling for Antimalarial Dr: 897 features

#### TRANSCRIPTOMICS
- Total datasets: 3
- Datasets with features: 1
- Sample datasets:
  - MYC-Driven Gliosis Impairs Neuron-Glia Communication in Amyo: 50344 features

### 3. Discovery Test Results

#### LIPIDOMICS
- Datasets processed: 1
- Signatures discovered: 1
- Sample signatures:
  - **AUTO_lipidomics_a60a2a28**: support=1, 740 components
    Features: Pi 20:4/o-18:0, Tg 16:0/22:5/p-18:0, Che 26:4, Pe 22:4/p-18:0, Sm d41:2

#### PROTEOMICS
- Datasets processed: 1
- Signatures discovered: 1
- Sample signatures:
  - **AUTO_proteomics_fd5a2043**: support=1, 897 components
    Features: RAB5B_HUMAN, DJC13_HUMAN, SNAG_HUMAN, TGM2_HUMAN, P5CR3_HUMAN

#### TRANSCRIPTOMICS
- Datasets processed: 1
- Signatures discovered: 1
- Sample signatures:
  - **AUTO_transcriptomics_df91840e**: support=1, 50344 components
    Features: GM43635, GM45709, GM5084, GM23217, 4930581F22Rik

## Key Findings

✅ Some datasets have features linked correctly.

✅ Discovery algorithm successfully found 3 signatures.

## Issues & Recommendations

- ✅ No critical issues found!

## Next Steps

1. If feature linking is broken: Fix ingestion pipeline to properly link features
2. If discovery finds 0 signatures: Test with lower thresholds or more homogeneous datasets
3. If discovery works: Test on larger dataset slices (10-20 datasets)
4. Validate discovered signatures manually for biological plausibility
