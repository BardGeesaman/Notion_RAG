# ALS Study Omics Data Inspection

## Study Imported
- **Study ID**: GSE275841
- **Title**: MYC-Driven Gliosis Impairs Neuron-Glia Communication in Amyotrophic Lateral Sclerosis
- **Dataset ID**: 69f55ca5-eef4-4195-966f-5bcb1f44136f

## ✅ Omics Data Confirmed!

### Data Type
**TRANSCRIPTOMICS (RNA-seq)**

- Expression profiling by high throughput sequencing
- RNA-seq experiment
- Quantitative gene expression measurements

### Available Data

#### 1. Expression Matrix (Series Matrix File)
- **Format**: Tab-delimited text file (.txt.gz)
- **Structure**: Gene IDs × Sample IDs matrix
- **Content**: Gene expression values (counts or normalized)
- **Location**: GEO FTP server
- **Size**: Available for download

#### 2. Sample Metadata
- **Cell Types**: 
  - Astrocytes (ACSA+)
  - Oligodendrocytes
  - Microglia
- **Genotypes**:
  - TDP-43_Q331K
  - SOD1
  - C9orf72
- **Time Points**: 3 months, etc.
- **Organism**: Mus musculus (mouse)

#### 3. Platform Information
- **Platform ID**: GPL24247
- **Type**: RNA-seq platform
- **Organism**: Mouse (taxid: 10090)

### Expression Data Structure

The Series Matrix file contains:
- **Header Section**: Metadata about the study
- **Data Table**: 
  - First column: Gene/transcript IDs
  - Subsequent columns: Sample IDs (GSM numbers)
  - Values: Expression measurements per gene per sample

### Data Access

**Series Matrix File URL**:
```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275841/matrix/GSE275841_series_matrix.txt.gz
```

### Sample Count
- Multiple samples across different conditions
- Cell type-specific expression profiles
- Genotype comparisons

## Summary

✅ **This study contains the omics data you're looking for!**

- Quantitative gene expression measurements
- Multiple experimental conditions
- Cell type-specific data
- Ready for:
  - Differential expression analysis
  - Feature extraction
  - Signature matching
  - Pathway analysis

## Next Steps

1. **Extract Expression Features**
   - Parse Series Matrix file
   - Extract gene IDs and expression values
   - Link to feature database

2. **Feature Linking**
   - Map genes to your feature database
   - Create feature-dataset associations

3. **Analysis**
   - Differential expression
   - Signature matching
   - Pathway enrichment

