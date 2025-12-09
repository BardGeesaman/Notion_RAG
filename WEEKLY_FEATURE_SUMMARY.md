# Amprenta RAG Platform: Functional Capabilities Summary

**Last Updated**: December 7, 2025  
**For**: Scientific Team  
**Update Frequency**: Weekly (Fridays)

---

## Current Functional Capabilities

### 1. Data Ingestion

#### 1.1 Multi-Omics Data Types Supported
- **Lipidomics**: CSV/TSV files with lipid species data
- **Metabolomics**: CSV/TSV files with metabolite measurements
- **Proteomics**: CSV/TSV files with protein quantification
- **Transcriptomics**: CSV/TSV files with gene expression data

#### 1.2 Automated Processing
- **Automatic Type Detection**: System identifies omics type from filename and file content
- **Batch Processing**: Ingest multiple files simultaneously
- **Automatic Linking**: Datasets automatically connected to relevant programs and experiments based on metadata

#### 1.3 Data Normalization
- **Feature Name Standardization**: Converts vendor-specific names to canonical formats
  - Example: `"CER 16:0"` becomes `"Cer(d18:1/16:0)"`
- **Cross-Platform Compatibility**: Harmonizes data from different instruments and vendors
- **Database Linking**: Automatic connection to HMDB, KEGG, UniProt identifiers

---

### 2. Signature Analysis

#### 2.1 Signature Definition
- **Multi-Omics Signatures**: Combine genes, proteins, metabolites, and lipids in one signature
- **Directional Information**: Specify up-regulation (↑) or down-regulation (↓) for each component
- **Component Weighting**: Assign importance scores to individual features

#### 2.2 Signature Discovery
- **Automatic Pattern Detection**: System identifies co-occurring features across datasets
- **Statistical Validation**: Ensures patterns are statistically significant
- **Candidate Ranking**: Prioritizes discovered signatures by confidence

#### 2.3 Signature Scoring
- **Dataset Matching**: Compare any dataset against your signature library
- **Overlap Analysis**: Calculate how many signature components are present
- **Direction Consistency**: Verify if up/down regulation matches expectations
- **Match Reporting**: Ranked list of matching datasets with detailed metrics

---

### 3. Semantic Search & Knowledge Retrieval

#### 3.1 Natural Language Queries
- **Ask Questions in Plain English**: 
  - "What do we know about ceramide dysregulation in ALS?"
  - "Which studies measured C16 ceramide in CSF?"
  - "Show me transcriptomics data related to sphingolipid metabolism"

#### 3.2 Search Across All Sources
- Literature pages (Zotero-imported papers)
- Experiment notebooks
- Dataset descriptions
- Signature definitions
- Email archives

#### 3.3 Advanced Filtering
- Filter by omics type (lipidomics, metabolomics, etc.)
- Filter by disease
- Filter by sample type (CSF, plasma, serum, tissue)
- Filter by model system (human, mouse, cell line)
- Combine multiple filters

---

### 4. AI-Powered Analysis

#### 4.1 Cross-Omics Reasoning
- **Program Summaries**: AI synthesizes all omics data for a research program
- **Signature Analysis**: AI explains how signatures appear across datasets
- **Feature Reports**: AI compiles all evidence for specific molecules
- **Dataset Context**: AI provides multi-omics context for individual datasets

#### 4.2 Intelligence Features
- **Convergent Findings**: Identifies patterns that appear across multiple omics types
- **Divergent Findings**: Highlights contradictory results that need investigation
- **Knowledge Gaps**: Points out what's missing or unclear
- **Mechanistic Insights**: Suggests connections between omics layers

---

### 5. Chemistry & Screening

#### 5.1 Compound Tracking
- **Compound Database**: Store chemical structures (SMILES, InChI)
- **Property Calculations**: Molecular weight, logP, H-bond donors/acceptors, etc.
- **Program Associations**: Link compounds to research programs

#### 5.2 High-Throughput Screening
- **Campaign Management**: Track screening campaigns with metadata
- **Hit Recording**: Store assay results with activity values and Z-scores
- **Hit Retrieval**: Query hits by campaign, activity threshold, or category
- **Compound-to-Hit Linking**: Connect screening hits to compound structures

---

### 6. Data Access Methods

#### 6.1 Command-Line Scripts
- **60+ Pre-Built Scripts**: Ready-to-use tools for common tasks
- **Examples**:
  - `ingest_lipidomics.py`: Ingest a lipidomics file
  - `batch_ingest_omics.py`: Ingest multiple files at once
  - `discover_signatures.py`: Find patterns across datasets
  - `rag_query.py`: Search your data with natural language
  - `generate_evidence_report.py`: Create cross-omics summaries

#### 6.2 REST API
- **Programmatic Access**: Build custom tools that interact with the platform
- **Available Endpoints**:
  - Compounds (list, get details, get linked programs)
  - Screening campaigns (list, get details, get hits)
  - Programs, Experiments, Datasets, Features, Signatures
- **Interactive Documentation**: Automatic API docs at `/docs` endpoint

#### 6.3 Python Library
- **Import and Use**: `from amprenta_rag.ingestion import ingest_lipidomics_file`
- **Programmatic Workflows**: Build custom analysis pipelines
- **Integration**: Embed platform capabilities in Jupyter notebooks or scripts

---

### 7. User Features

#### 7.1 Automatic Organization
- **Smart Linking**: Datasets automatically connected to programs/experiments
- **Confidence Scoring**: System only makes links when confident (≥80% certainty)
- **Ambiguity Handling**: Avoids incorrect links when metadata is unclear

#### 7.2 Quality Control
- **Error Reporting**: Clear messages when files have issues
- **Validation**: Automatic checks for required fields and data formats
- **Recovery**: Failed files don't stop batch processing

#### 7.3 Configuration
- **Environment-Based Settings**: Configure via `.env` file (no code changes)
- **Adjustable Thresholds**: Tune signature matching, auto-linking confidence, etc.
- **Enable/Disable Features**: Turn on/off auto-linking, caching, etc.

---

## Planned Features (Roadmap)

### Next Sprint (December 2025)

#### Visualization Dashboards
- **Interactive Plots**: Volcano plots, heatmaps, PCA plots for datasets
- **Signature Visualization**: Network diagrams showing signature components
- **Program Overview Dashboards**: Visual summary of all program data

#### Data Quality Checks
- **Automated Validation**: Flag datasets with missing data, outliers, or inconsistencies
- **Quality Scores**: Assign quality ratings to datasets
- **Quality Reports**: Generate summaries of data quality issues

#### Pathway Analysis
- **KEGG Integration**: Map features to KEGG pathways
- **Reactome Integration**: Map features to Reactome pathways
- **Enrichment Analysis**: Identify over-represented pathways in datasets
- **Pathway Visualization**: Visual pathway maps with highlighted features

---

### Next Month (January 2026)

#### Public Repository Integration
- **Metabolomics Workbench**: Auto-import datasets by study ID
- **GEO (Gene Expression Omnibus)**: Auto-import transcriptomics datasets
- **PRIDE**: Auto-import proteomics datasets
- **MetaboLights**: Auto-import metabolomics datasets

#### Advanced Statistical Analysis
- **Differential Expression**: Built-in DE analysis for transcriptomics
- **Statistical Tests**: T-tests, ANOVA, Mann-Whitney for omics data
- **Multiple Testing Correction**: FDR, Bonferroni corrections
- **Correlation Analysis**: Find correlated features across datasets

#### Experimental Design Tracking
- **Study Structure Capture**: Automatically extract case/control, time course, intervention groups
- **Design Types Supported**: Case vs control, time course, intervention groups, dose response, multi-factorial
- **Auto-Extraction**: Parse design from repository metadata (GEO, Metabolomics Workbench)
- **Design-Aware Statistics**: Statistical tests that respect experimental design
- **Time Course Analysis**: Repeated measures, trend analysis over time
- **Group Comparisons**: Multi-group ANOVA for intervention studies

#### Data Export
- **ISA-Tab Format**: Export to standardized metadata format
- **MAGE-TAB Format**: Export microarray/transcriptomics data
- **CSV/Excel Export**: Download processed data in common formats
- **Report Generation**: PDF reports with figures and tables

---

### Next Quarter (Q1 2026)

#### Multi-User Support
- **User Accounts**: Individual logins with passwords
- **Permissions System**: Control who can view/edit different data
- **Team Workspaces**: Separate spaces for different research groups
- **Audit Logs**: Track who did what and when

#### Collaboration Features
- **Comments & Annotations**: Add notes to datasets, signatures, experiments
- **Sharing**: Share specific datasets or analyses with colleagues
- **Notifications**: Get alerts when new data matches your interests
- **Version History**: Track changes to signatures and datasets over time

#### Advanced Queries
- **Similarity Search**: Find datasets similar to a given dataset
- **Temporal Analysis**: Track how signatures change over time
- **Cohort Comparison**: Compare groups of datasets (e.g., responders vs non-responders)
- **Custom Filters**: Save and reuse complex filter combinations

---

### Future Considerations (Q2 2026+)

#### Machine Learning
- **Predictive Models**: Train models to predict outcomes from omics data
- **Biomarker Discovery**: Automated identification of discriminative features
- **Clustering**: Unsupervised grouping of similar datasets
- **Dimensionality Reduction**: t-SNE, UMAP for visualization

#### Integration Enhancements
- **LIMS Integration**: Connect to laboratory information management systems
- **ELN Integration**: Deeper integration with electronic lab notebooks
- **Cloud Deployment**: AWS/Azure/GCP hosting options
- **Mobile Access**: Optimized interface for phones and tablets

#### Advanced Signature Features
- **Temporal Signatures**: Signatures that capture time-series patterns
- **Conditional Signatures**: "If X is up, then Y should be down" logic
- **Signature Evolution**: Track how signatures change as you add more data
- **Signature Prediction**: Predict which signatures a dataset will match before testing

---

## How to Use This Document

### For Scientists
- **What can I do now?**: See "Current Functional Capabilities" (sections 1-7)
- **What's coming soon?**: See "Planned Features - Next Sprint"
- **What about later?**: See "Next Month" and "Next Quarter" sections

### For Planning
- **This Week**: Features marked "Next Sprint" are in active development
- **This Month**: Features marked "Next Month" are being designed/scoped
- **This Quarter**: Features marked "Next Quarter" are in planning phase

### For Questions
- **"Can the system do X?"**: Check sections 1-7
- **"When will the system do X?"**: Check roadmap sections
- **"How do I do X?"**: See documentation in `docs/` directory

---

## Feature Request Process

Have an idea for a new feature? Here's how to suggest it:

1. **Check this document**: Is it already planned?
2. **Submit a request**: Create a GitHub issue or message the team channel
3. **Provide context**: Explain the scientific use case
4. **Expected timeline**: Team will respond with priority/timeline within 1 week

---

## Weekly Update Log

### Week of December 7, 2025
**New Capabilities**:
- ✅ Automatic program/experiment linking with confidence scoring
- ✅ REST API for compounds and screening data
- ✅ Batch ingestion with automatic omics type detection

**Documentation Added**:
- Auto-Linking Guide (how confidence scoring works)
- API Reference updates (new endpoints)
- Feature Caching Guide (optimization features)

**Next Week Focus**: Visualization dashboard prototypes

---

### Week of December 1, 2025
**New Capabilities**:
- ✅ Multi-omics signature scoring engine
- ✅ Chemistry compound database integration
- ✅ HTS campaign tracking

**Documentation Added**:
- Signature scoring documentation
- Chemistry database setup guide

**Next Week Focus**: Automatic linking implementation

---

### Week of November 24, 2025
**New Capabilities**:
- ✅ Cross-omics AI reasoning (GPT-4 integration)
- ✅ Signature discovery automation
- ✅ Batch processing framework

**Documentation Added**:
- Cross-omics reasoning examples
- Batch ingestion guide

**Next Week Focus**: Production hardening and error handling

---

## Document Maintenance

**Updated By**: Development team  
**Review Schedule**: Every Friday afternoon  
**Next Review**: December 14, 2025  
**Feedback**: Send questions/suggestions to team Slack channel

---

## Quick Reference

**Can I...?**
- ✅ Ingest lipidomics/metabolomics/proteomics/transcriptomics data
- ✅ Define multi-omics signatures
- ✅ Search data with natural language questions
- ✅ Get AI summaries of cross-omics evidence
- ✅ Track chemistry compounds and screening hits
- ✅ Access data via REST API
- ✅ Run batch operations on multiple files
- ⏳ Visualize data (coming December 2025)
- ⏳ Import from public repositories (coming January 2026)
- ⏳ Export to ISA-Tab (coming January 2026)
- ⏳ Multi-user access (coming Q1 2026)

**Legend**: ✅ Available Now | ⏳ Planned

---

**Questions?** See `docs/` directory for detailed guides or ask in team channel.
