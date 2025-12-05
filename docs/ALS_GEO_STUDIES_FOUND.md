# ALS-Related GEO Studies Found

## Search Date
2025-12-04

## Search Results

### Total Studies Found
**29 unique ALS-related GEO studies**

### Search Terms Used
1. "ALS"
2. "amyotrophic lateral sclerosis" 
3. "motor neuron"

## Top Studies Found

### 1. GSE275841
- **Title**: MYC-Driven Gliosis Impairs Neuron-Glia Communication in Amyotrophic Lateral Sclerosis
- **Organism**: Mus musculus
- **Focus**: Glial cell dysfunction in ALS

### 2. GSE252892
- **Title**: TDP-43 nuclear loss in FTD-ALS causes widespread alternative polyadenylation
- **Organism**: Homo sapiens
- **Focus**: TDP-43 protein and RNA processing

### 3. GSE307044
- **Title**: Optineurin Shapes Basal and LPS-Induced Transcriptomes in BV2 Microglia
- **Organism**: Mus musculus
- **Focus**: OPTN gene in neurodegenerative diseases

### 4. GSE309506
- **Title**: Pharmacological rescue of motor neuron health by riluzole in patient-derived motor neurons
- **Organism**: Homo sapiens
- **Focus**: Riluzole treatment in ALS

### 5. GSE307456
- **Title**: IsomiR Utility in ALS Prognostication [control]
- **Organism**: Homo sapiens
- **Focus**: miRNA biomarkers for ALS

## All Study IDs Found

1. GSE252892 - TDP-43 nuclear loss in FTD-ALS
2. GSE265944
3. GSE275841 - MYC-Driven Gliosis in ALS
4. GSE278666
5. GSE288856
6. GSE289958
7. GSE290979
8. GSE290980
9. GSE298532
10. GSE300606
11. GSE301963
12. GSE302051
13. GSE302208
14. GSE302295
15. GSE302319
16. GSE302320
17. GSE302449
18. GSE305831
19. GSE306766
20. GSE307044 - Optineurin in Microglia
21. GSE307375
22. GSE307456 - IsomiR in ALS Prognostication
23. GSE308008
24. GSE308151
25. GSE309401
26. GSE309506 - Riluzole in Motor Neurons
27. GSE309664
28. GSE309740
29. GSE310834

## Import Instructions

To import any of these studies:

```bash
python scripts/harvest_repository_study.py \
    --study-id GSE<ID> \
    --repository GEO \
    --create-postgres \
    --ingest
```

Example:
```bash
python scripts/harvest_repository_study.py \
    --study-id GSE275841 \
    --repository GEO \
    --create-postgres \
    --ingest
```

