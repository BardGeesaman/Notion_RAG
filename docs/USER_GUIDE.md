# Amprenta RAG User Guide

**Complete Guide to Multi-Omics Research Platform**

Version 1.0 | December 2025

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Getting Started](#2-getting-started)
3. [Data Management](#3-data-management)
4. [Signatures & Pattern Discovery](#4-signatures--pattern-discovery)
5. [Analysis & Insights](#5-analysis--insights)
6. [Chemistry & High-Throughput Screening](#6-chemistry--high-throughput-screening)
7. [Jupyter Integration](#7-jupyter-integration)
8. [Advanced Features](#8-advanced-features)
9. [Troubleshooting & FAQs](#9-troubleshooting--faqs)

---

## 1. Introduction

### What is Amprenta RAG?

Amprenta RAG is a comprehensive knowledge management platform designed for multi-omics research, with a particular focus on ceramide and sphingolipid neurodegeneration studies. The platform combines traditional database management with cutting-edge AI capabilities to help researchers discover patterns, generate insights, and accelerate scientific discovery.

### What Makes Amprenta Different?

**Multi-Omics Integration**: Unlike single-modality tools, Amprenta seamlessly integrates transcriptomics, proteomics, metabolomics, and lipidomics data in one unified system. You can ask questions that span multiple data types and get coherent answers.

**AI-Powered Search**: Using Retrieval-Augmented Generation (RAG), the system doesn't just search for keywords—it understands the semantic meaning of your questions and finds relevant information even when exact terms don't match.

**Signature-Based Discovery**: The platform helps you discover and validate multi-omics signatures—patterns of features that consistently appear together across datasets. These signatures can reveal biological mechanisms and predict outcomes.

**Chemistry Integration**: Built-in chemistry tools let you manage compound libraries, analyze structure-activity relationships (SAR), and connect molecular structures to biological outcomes.

**Interactive Analysis**: Through integrated Jupyter notebooks and Voila dashboards, you can explore data interactively, create custom visualizations, and generate publication-ready reports.

### Key Capabilities

- **Unified Data Repository**: Store and organize experiments, datasets, and results from multiple omics platforms
- **Intelligent Search**: Ask questions in natural language and get contextual answers with citations
- **Pattern Discovery**: Automatically detect signatures and score new datasets against known patterns
- **Cross-Omics Analysis**: Connect genes, proteins, metabolites, and lipids through pathway enrichment
- **Chemistry Tools**: Structure search, SAR analysis, compound registration with property calculation
- **Collaborative Science**: Share findings, track protocols, maintain audit trails
- **Reproducible Research**: Jupyter integration with version control and parameterized reports

### Who Should Use This Guide?

This guide is written for:

- **Computational Biologists**: Analyzing multi-omics datasets
- **Chemists**: Managing compound libraries and SAR studies
- **Data Scientists**: Building custom analysis pipelines
- **Principal Investigators**: Overseeing research programs
- **Lab Scientists**: Uploading experimental results and retrieving insights

No programming knowledge is required for most features, though Jupyter integration enables advanced customization for those who want it.

---

## 2. Getting Started

### System Requirements

**Supported Browsers**:
- Chrome 90+ (recommended)
- Firefox 88+
- Safari 14+
- Edge 90+

**Network**: Stable internet connection (the platform is web-based)

**Recommended Screen Resolution**: 1920x1080 or higher for optimal dashboard viewing

### First Login

When you first access Amprenta RAG, you'll see the login screen. Use the credentials provided by your system administrator.

**Default URL**: `http://localhost:8501` (local deployment) or your organization's custom domain.

After logging in successfully, you'll land on the **Home Dashboard**—your central hub for all activities.

### Understanding the Home Dashboard

The Home Dashboard provides a bird's-eye view of your research environment:

**Top Navigation Bar**:
- **Home**: Return to the main dashboard
- **Data**: Access programs, experiments, and datasets
- **Analysis**: Signatures, pathways, RAG search
- **Chemistry**: Compounds, HTS campaigns, SAR
- **Jupyter**: Open notebooks and dashboards
- **Admin**: System settings (if you have permissions)

**Quick Stats Panel** (top of page):
- Total programs
- Active experiments
- Datasets ingested
- Signatures discovered
- Recent activity feed

**Recent Activity**:
See the latest actions across your workspace—new datasets uploaded, signatures created, analyses run. This helps you track what your team is working on.

**Favorite Shortcuts**:
Pin frequently used pages for quick access.

### Navigation Patterns

**Breadcrumb Navigation**: At the top of each page, breadcrumbs show your current location. Click any breadcrumb to navigate back.

Example: `Home > Programs > Neurodegeneration Study > Experiment NRS-001`

**Search Bar**: The global search bar (top right) lets you quickly find programs, experiments, datasets, compounds, or signatures by name or ID.

**Sidebar**: Many pages have a collapsible sidebar with filters, options, or contextual help.

### Understanding Data Hierarchy

Amprenta organizes research in a logical hierarchy:

```
Program (e.g., "Ceramide Neurodegeneration Project")
  └─ Experiment (e.g., "Mouse Model Lipidomics")
      └─ Dataset (e.g., "Brain Tissue LC-MS Run 1")
          └─ Features (e.g., individual lipids, genes, proteins)
```

**Programs** represent large research initiatives. They contain multiple experiments.

**Experiments** are specific studies within a program (e.g., a particular cohort, time point, or condition).

**Datasets** are the actual data files—one per analytical run or assay.

**Features** are the measured entities—genes, proteins, metabolites, lipids, etc.

This hierarchy helps you organize data logically and maintain context when analyzing results.

### Your First Task: Create a Program

Let's walk through creating your first program step by step.

1. **Navigate to Programs**:
   Click "Data" in the top menu, then select "Programs" from the dropdown.

2. **Click "New Program"**:
   You'll see a green "New Program" button in the top-right corner.

3. **Fill in Program Details**:
   - **Name**: Give your program a clear, descriptive name (e.g., "ALS Biomarker Discovery")
   - **Description**: Explain the program's goals, scope, and expected outcomes (markdown supported)
   - **Principal Investigator**: Your name or the lead researcher
   - **Status**: Active, Planned, or Completed
   - **Start Date**: When the program began
   - **Tags**: Add keywords for easy searching (e.g., "ALS", "metabolomics", "biomarkers")

4. **Click "Create Program"**:
   The system creates your program and takes you to the Program Detail page.

**Congratulations!** You've created your first program. Now you're ready to add experiments and data.

---

## 3. Data Management

### Understanding Multi-Omics Data Types

Amprenta supports four major omics modalities, each with specialized handling:

#### Transcriptomics (Gene Expression)

**What it is**: Measurement of RNA levels (mRNA, microRNA, etc.) to understand which genes are active.

**Common platforms**:
- RNA-seq (Illumina, PacBio)
- Microarray (Affymetrix, Agilent)

**What Amprenta captures**:
- Gene symbols (HGNC standard)
- Ensembl IDs
- Log2 fold changes
- P-values and adjusted p-values
- FPKM/TPM values

**File formats accepted**: CSV, TSV, Excel, GEO SOFT files

#### Proteomics (Protein Abundance)

**What it is**: Measurement of protein levels and post-translational modifications.

**Common platforms**:
- LC-MS/MS
- TMT labeling
- LFQ (label-free quantification)

**What Amprenta captures**:
- UniProt IDs
- Protein names
- Abundance ratios
- Peptide counts
- Phosphorylation sites

**File formats accepted**: CSV, TSV, MaxQuant output, Proteome Discoverer files

#### Metabolomics (Small Molecules)

**What it is**: Measurement of small molecule metabolites (sugars, amino acids, etc.).

**Common platforms**:
- GC-MS
- LC-MS
- NMR

**What Amprenta captures**:
- HMDB IDs
- KEGG compound IDs
- m/z values
- Retention times
- Chemical formulas

**File formats accepted**: CSV, TSV, mzTab, XCMS output, mwTab files

#### Lipidomics (Lipid Species)

**What it is**: Specialized metabolomics focused on lipids.

**Common platforms**:
- LC-MS (lipidomics mode)
- Shotgun lipidomics

**What Amprenta captures**:
- Lipid Maps IDs
- Systematic lipid names (e.g., "PC(16:0/18:1)")
- Sum compositions
- Lipid classes
- Abundance values

**File formats accepted**: CSV, TSV, LipidSearch output

### Creating an Experiment

Experiments are where the work happens. Each experiment represents a specific study within your program.

**Step 1: Navigate to Your Program**

From the Programs list, click on the program you created earlier.

**Step 2: Click "New Experiment"**

On the program detail page, click the "New Experiment" button.

**Step 3: Fill in Experiment Details**

- **Name**: Descriptive name (e.g., "ALS Patient Plasma Metabolomics - Cohort 1")
- **Description**: Study design, sample information, analytical methods
- **Experiment Type**: Choose the primary omics type
- **Sample Type**: Tissue type (brain, plasma, cells, etc.)
- **Species**: Human, mouse, rat, etc.
- **Disease/Condition**: Relevant disease or experimental condition
- **Technology**: Platform used (Illumina NextSeq, Agilent 6550, etc.)
- **Date**: When the experiment was conducted

**Step 4: Add Protocol Information (Optional but Recommended)**

Good science requires documentation. Add:
- Sample preparation steps
- Instrument settings
- Data processing methods
- Quality control measures

This information becomes searchable and helps with reproducibility.

**Step 5: Save the Experiment**

Click "Create Experiment". You now have an experiment ready to receive datasets.

### Uploading Datasets

This is where your data enters the system. Amprenta makes uploading and processing data as painless as possible.

#### Preparing Your Data File

**Accepted Formats**: CSV, TSV (tab-separated), Excel (.xlsx)

**Required Columns**:

For **Transcriptomics**:
- `gene_symbol` or `gene_id` (required)
- `log2_fc` (log2 fold change)
- `p_value` or `adj_p_value`
- `expression_value` (optional: FPKM, TPM, raw counts)

For **Proteomics**:
- `protein_id` or `uniprot_id` (required)
- `protein_name`
- `log2_ratio` or `fold_change`
- `p_value` or `adj_p_value`
- `peptide_count` (optional)

For **Metabolomics**:
- `compound_name` or `hmdb_id` (required)
- `mz` (mass-to-charge ratio)
- `rt` (retention time)
- `fold_change`
- `p_value`

For **Lipidomics**:
- `lipid_name` (systematic nomenclature preferred)
- `lipid_class` (e.g., PC, PE, SM, Cer)
- `fold_change`
- `p_value`

**Example CSV** (Metabolomics):

```csv
compound_name,hmdb_id,mz,rt,fold_change,p_value
L-Glutamate,HMDB0000148,148.0604,2.45,1.52,0.001
Sphingosine,HMDB0000252,300.2897,15.2,2.34,0.0001
Ceramide(d18:1/16:0),HMDB0004949,538.5128,18.7,1.87,0.002
```

**Tips**:
- Use standard IDs when possible (HGNC gene symbols, UniProt IDs, HMDB IDs)
- Include p-values for statistical significance filtering
- Avoid special characters in column names
- UTF-8 encoding recommended

#### Upload Process

**Step 1: Navigate to Your Experiment**

From the experiment detail page, click "Upload Dataset".

**Step 2: Select Your File**

Click "Choose File" and select your prepared CSV/TSV/Excel file.

**Step 3: Configure Import Settings**

The system shows a preview of your data and asks you to map columns:

- **Feature ID Column**: Which column contains feature identifiers?
- **Feature Name Column**: Human-readable names (optional)
- **Value Column**: Fold changes, expression values, etc.
- **P-Value Column**: Statistical significance (optional)
- **Additional Metadata**: Any extra columns to preserve

The system automatically detects omics type and suggests mappings. Review and confirm.

**Step 4: Set Dataset Metadata**

- **Dataset Name**: Auto-generated from filename, but you can edit
- **Description**: Brief note about this specific dataset
- **Date Acquired**: When the data was generated
- **Quality Status**: Draft, Reviewed, or Approved

**Step 5: Click "Import"**

The system processes your file:
1. Validates format
2. Normalizes feature identifiers (converts to standard IDs)
3. Extracts features
4. Links to external databases (KEGG, Reactome, etc.)
5. Calculates summary statistics
6. Indexes for search

You'll see a progress indicator. Large files (10,000+ features) may take 1-2 minutes.

**Step 6: Review Import Summary**

After import completes, you see:
- Number of features imported
- Number of features normalized (successfully mapped to standard IDs)
- Warnings (if any IDs couldn't be mapped)
- Quality metrics (missing values, outliers)

Click "View Dataset" to explore your data.

### Viewing and Exploring Datasets

The **Dataset Detail Page** is your window into the data:

**Summary Panel** (top):
- Feature count
- Feature type breakdown (genes, proteins, metabolites)
- Date imported
- Status

**Feature Table**:
A sortable, filterable table showing all features in the dataset:

- **Feature ID**: Standard identifier (clickable to see details)
- **Feature Name**: Human-readable name
- **Value**: Fold change or abundance
- **P-Value**: Statistical significance
- **Direction**: Up, Down, or Neutral (based on thresholds)

**Filters**:
- **Significance**: Show only features with p < 0.05
- **Fold Change**: Filter by magnitude (e.g., |FC| > 1.5)
- **Direction**: Up-regulated, down-regulated, or both
- **Feature Type**: Filter by genes, proteins, metabolites, etc.

**Export Options**:
- Download filtered features as CSV
- Export full dataset with metadata
- Generate feature list for external tools

**Visualizations**:
- **Volcano Plot**: Fold change vs. p-value (interactive)
- **Distribution**: Histogram of values
- **Feature Type Breakdown**: Pie chart showing composition

**Actions**:
- **Create Signature from Dataset**: Turn significant features into a signature
- **Compare to Other Datasets**: Side-by-side comparison
- **Analyze with RAG**: Ask questions about this dataset
- **Open in Jupyter**: Launch notebook with this dataset pre-loaded

### Importing from Public Repositories

Amprenta connects directly to major public data repositories, letting you import published datasets with one click.

#### Supported Repositories

**GEO (Gene Expression Omnibus)**:
- Transcriptomics (RNA-seq, microarray)
- Search by GEO accession (GSE, GSM)
- Metadata automatically extracted

**PRIDE (Proteomics Identifications Database)**:
- Proteomics data
- Search by PXD accession
- Downloads processed result files

**MetaboLights**:
- Metabolomics studies
- Search by MTBLS accession
- Extracts MAF (Metabolite Assignment Files)

**Metabolomics Workbench**:
- Comprehensive metabolomics repository
- Search by ST accession
- Parses mwTab format automatically

#### How to Import from a Repository

**Step 1: Navigate to Repository Import**

From the main menu, click "Data" → "Import from Repository".

**Step 2: Choose Repository**

Select the repository (GEO, PRIDE, MetaboLights, or Metabolomics Workbench).

**Step 3: Search or Enter Accession**

You can either:
- **Search**: Enter keywords (disease, tissue, author) to find relevant studies
- **Direct Import**: Paste an accession number if you already know what you want

Example: Search for "ALS motor neuron" in GEO.

**Step 4: Review Study Metadata**

The system fetches metadata and shows:
- Title
- Authors
- Publication date
- Abstract
- Sample count
- Platform used

Review to confirm this is the study you want.

**Step 5: Select Datasets**

Many studies contain multiple datasets (different samples, conditions, or time points). Select which ones to import.

**Step 6: Assign to Program and Experiment**

Choose:
- Which program this data belongs to
- Create a new experiment or add to existing one

**Step 7: Click "Import"**

The system:
1. Downloads raw data files from the repository
2. Parses format-specific files (SOFT, mzTab, mwTab)
3. Normalizes features to standard IDs
4. Creates datasets in your experiment
5. Links to original repository for provenance

Import may take several minutes for large studies.

**Step 8: Review Imported Datasets**

Once complete, navigate to your experiment to see the newly imported datasets. They're tagged with repository provenance (e.g., "Source: GEO GSE12345").

### Batch Import

For high-throughput labs generating many datasets, batch import saves time.

**Prepare a Batch File** (CSV format):

```csv
file_path,experiment_id,dataset_name,omics_type
data/sample1.csv,EXP-001,Sample 1 Metabolomics,metabolomics
data/sample2.csv,EXP-001,Sample 2 Metabolomics,metabolomics
data/sample3.csv,EXP-001,Sample 3 Metabolomics,metabolomics
```

**Upload the Batch File**:

From "Data" → "Batch Import", upload your CSV. The system processes all datasets automatically and reports success/failure for each.

### Data Quality and Validation

Amprenta performs automatic quality checks on imported data:

**Identifier Validation**:
- Are gene symbols valid HGNC symbols?
- Do HMDB IDs exist in HMDB database?
- Are UniProt IDs current?

**Statistical Validation**:
- Are p-values in valid range [0, 1]?
- Are fold changes realistic?
- Are there obvious outliers?

**Completeness**:
- How many missing values?
- Are required columns present?
- Is metadata complete?

**Quality Score**:
Each dataset receives a quality score (0-100):
- 90-100: Excellent (ready for analysis)
- 70-89: Good (minor issues)
- 50-69: Fair (review warnings)
- <50: Poor (needs attention)

View quality reports on the dataset detail page. Address any warnings before publishing results.

---

## 4. Signatures & Pattern Discovery

### What Are Signatures?

A **signature** is a list of features (genes, proteins, metabolites, lipids) that consistently co-occur in certain conditions. Signatures represent biological mechanisms, disease states, or treatment responses.

**Examples**:

- **Alzheimer's Disease Signature**: 50 genes consistently upregulated in AD patients
- **Ceramide Stress Response**: 20 lipid species that increase together under oxidative stress
- **Mitochondrial Dysfunction Signature**: Proteins and metabolites indicating impaired mitochondria

Signatures are powerful because:
1. They capture multi-feature patterns (not just individual markers)
2. They're validated across multiple datasets
3. They can predict outcomes in new data
4. They reveal mechanisms

### Creating a Signature Manually

Sometimes you know exactly which features define a signature (e.g., from literature or expert knowledge).

**Step 1: Navigate to Signatures**

Click "Analysis" → "Signatures" in the main menu.

**Step 2: Click "New Signature"**

You'll see the Signature Builder interface.

**Step 3: Basic Information**

- **Name**: Descriptive name (e.g., "Neuroinflammation Gene Panel")
- **Description**: What does this signature represent? (markdown supported)
- **Signature Type**: Gene, Protein, Metabolite, Lipid, or Multi-Omics
- **Source**: Literature, expert opinion, discovered, etc.
- **Tags**: Keywords for searching

**Step 4: Add Features**

There are three ways to add features:

**Method 1: Manual Entry**
- Click "Add Feature"
- Enter feature ID (gene symbol, HMDB ID, etc.)
- Set weight (importance, usually 1.0 for equal weighting)
- Set direction: Up, Down, or Neutral
- Repeat for all features

**Method 2: Paste List**
- Click "Paste List"
- Paste a list of IDs (one per line)
- System validates and adds all at once

**Method 3: Import from File**
- Click "Import CSV"
- Upload a CSV with columns: `feature_id`, `weight`, `direction`

**Step 5: Set Significance Thresholds**

- **Default P-Value Threshold**: 0.05 (features below this are "significant")
- **Default Fold Change Threshold**: 1.5 (magnitude cutoff)

These control how strictly matches are scored.

**Step 6: Review and Save**

Review the feature list, then click "Create Signature". Your signature is now ready to match against datasets.

### Creating a Signature from a Dataset

The most common way to create signatures: extract them from differential expression results.

**Step 1: Open a Dataset**

Navigate to any dataset detail page.

**Step 2: Click "Create Signature from This Dataset"**

A dialog appears with options.

**Step 3: Define Selection Criteria**

Choose which features to include:

- **P-Value Cutoff**: e.g., p < 0.01 (only significant features)
- **Fold Change Cutoff**: e.g., |FC| > 2.0 (only large changes)
- **Direction**: Up-regulated only, down-regulated only, or both
- **Top N Features**: Limit to top features by significance

**Example**: "Include genes with p < 0.01 and |log2 FC| > 1.5"

**Step 4: Preview**

The system shows how many features match your criteria. Adjust thresholds if needed.

**Step 5: Name and Save**

Give your signature a name (auto-filled based on dataset name) and click "Create".

The new signature is created with:
- Features meeting your criteria
- Weights based on fold change magnitude
- Directions based on sign of fold change
- Link to source dataset for provenance

### Scoring Datasets Against Signatures

This is where signatures become useful—you can see if a signature "matches" a new dataset.

**Step 1: Navigate to Signature Detail Page**

Click on any signature from the Signatures list.

**Step 2: Click "Score Against Datasets"**

**Step 3: Select Datasets to Score**

Choose one or more datasets to compare:
- All datasets in an experiment
- Specific datasets you select
- Datasets from a particular program

**Step 4: Click "Calculate Scores"**

The system compares each dataset to the signature using statistical methods:

**Scoring Methods**:

1. **Feature Overlap**: How many signature features appear in the dataset?
2. **Direction Concordance**: Do features change in the same direction?
3. **Weighted Score**: Combines overlap, direction, and feature importance
4. **Enrichment P-Value**: Statistical significance (hypergeometric test)

**Step 5: Review Results**

You see a results table:

| Dataset | Score | P-Value | Features Matched | Direction Concordance |
|---------|-------|---------|------------------|----------------------|
| Dataset A | 0.85 | 0.001 | 42/50 (84%) | 38/42 (90%) |
| Dataset B | 0.32 | 0.15 | 18/50 (36%) | 10/18 (56%) |

**High scores** (>0.7) indicate strong signature presence.
**Low p-values** (<0.05) indicate statistical significance.
**High direction concordance** (>80%) indicates biological relevance.

**Step 6: Explore Match Details**

Click on any scored dataset to see:
- Which features matched
- Which features didn't match
- Feature-by-feature comparison
- Contribution of each feature to the score

### Signature Match Explainability

Understanding *why* a signature matches (or doesn't match) is crucial. Amprenta provides detailed explainability tools.

**Feature Contribution View**:

For each matched feature, see:
- **Contribution Score**: How much this feature contributed to the overall match
- **Sign Match**: Does the direction agree? (e.g., both upregulated)
- **Magnitude**: Fold change in dataset vs. expected direction
- **Statistical Significance**: Is this feature significant in the dataset?

**Visual Explanation**:

**Lollipop Chart**: Shows top contributing features (positive contributions in green, negative in red)

**Direction Concordance Bar**: Visual percentage of features with matching direction

**Pathway Context**: If features belong to known pathways, see which pathways are enriched

**Example Interpretation**:

> "This dataset scores 0.78 against the Mitochondrial Dysfunction signature. The top contributors are:
> - COX4I1 (cytochrome c oxidase): strongly downregulated as expected (+0.15 contribution)
> - ATP5A1 (ATP synthase): downregulated as expected (+0.12 contribution)
> - VDAC1 (voltage-dependent anion channel): upregulated, but signature expects downregulation (-0.05 contribution)
>
> Overall, 85% of signature features show concordant direction changes, indicating likely mitochondrial impairment in this sample."

### Cross-Omics Signatures

One of Amprenta's unique strengths: signatures that span multiple omics types.

**Example**: A "Ceramide Stress" signature might include:
- **Genes**: SMPD1, CERS2, SGPP1 (ceramide metabolism enzymes)
- **Proteins**: SMPD1 protein abundance
- **Metabolites**: Sphingosine, sphingosine-1-phosphate
- **Lipids**: Ceramide(d18:1/16:0), Ceramide(d18:1/24:0)

**Creating a Cross-Omics Signature**:

Follow the same process as creating a signature manually, but:
1. Set "Signature Type" to "Multi-Omics"
2. Add features from different omics types
3. Assign weights reflecting relative importance across modalities

**Scoring Cross-Omics Signatures**:

When scoring, Amprenta:
1. Scores each omics type separately
2. Combines scores with weights
3. Reports overall match plus breakdown by omics type

This reveals whether all omics layers agree or if there's dysregulation at specific levels (e.g., gene expression changes but protein doesn't follow).

### Automatic Signature Discovery

Amprenta can discover signatures automatically by finding features that co-occur across multiple datasets.

**Step 1: Navigate to "Signature Discovery"**

From "Analysis" → "Discover Signatures".

**Step 2: Select Datasets for Mining**

Choose datasets to analyze:
- All datasets in a program
- Datasets with specific tags
- Datasets from a particular condition or tissue type

**Step 3: Set Discovery Parameters**

- **Minimum Datasets**: How many datasets must a pattern appear in? (e.g., 3)
- **Minimum Features**: How many features define a signature? (e.g., 10)
- **Correlation Threshold**: How tightly must features co-vary? (e.g., 0.6)
- **P-Value Threshold**: Statistical significance for inclusion (e.g., 0.01)

**Step 4: Run Discovery**

The algorithm:
1. Identifies frequently co-occurring features across datasets
2. Clusters features by correlation
3. Validates statistical significance
4. Proposes candidate signatures

This may take several minutes for large dataset collections.

**Step 5: Review Candidate Signatures**

The system presents candidate signatures with:
- Proposed name (based on enriched pathways or GO terms)
- Feature list
- Prevalence (in how many datasets)
- Enrichment analysis (which pathways/functions)

**Step 6: Accept or Refine**

Review each candidate:
- Accept as-is
- Edit feature list or name
- Reject if not biologically meaningful

Accepted signatures are added to your signature library.

### Signature Validation

Before using a signature for predictions, validate it on independent datasets.

**Validation Workflow**:

1. **Training Set**: Create signature from subset of datasets
2. **Test Set**: Score signature against held-out datasets
3. **Performance Metrics**:
   - Sensitivity: How many true positives did it catch?
   - Specificity: How many negatives did it correctly exclude?
   - AUC: Area under ROC curve

**Example**:

> "Created Alzheimer's Disease signature from 20 AD patient datasets. Tested on 10 independent AD datasets and 10 healthy controls. Achieved:
> - Sensitivity: 90% (9/10 AD datasets matched)
> - Specificity: 100% (0/10 controls matched)
> - AUC: 0.95
>
> Signature is validated for AD detection."

---

## 5. Analysis & Insights

### RAG (Retrieval-Augmented Generation) Search

RAG is Amprenta's most powerful feature—it lets you ask questions in natural language and get intelligent, contextualized answers based on your data.

#### How RAG Works (Simplified)

Traditional search looks for exact keyword matches. RAG:
1. Understands the *semantic meaning* of your question
2. Searches both structured data (datasets, signatures) and unstructured text (descriptions, protocols, literature)
3. Retrieves relevant context
4. Generates a coherent answer with citations

**Example**:

**Question**: "What ceramide species are elevated in our ALS patient samples?"

**Traditional keyword search** would find documents containing "ceramide", "ALS", and "patient" but wouldn't understand the relationship.

**RAG search**:
1. Understands you're asking about *lipid abundance* in a specific *disease context*
2. Retrieves relevant datasets tagged with ALS
3. Finds features classified as ceramides with elevated fold changes
4. Generates: "In ALS patient samples, the following ceramide species are significantly elevated: Ceramide(d18:1/16:0) (FC=2.1, p=0.003), Ceramide(d18:1/24:0) (FC=1.8, p=0.01). These findings are consistent across 3 datasets in the ALS Motor Neuron Study [link]. See Dataset Detail for full results."

#### Using RAG Search

**Step 1: Navigate to RAG Search**

Click "Analysis" → "Ask a Question" or use the search bar with RAG mode enabled.

**Step 2: Enter Your Question**

Type naturally, as if asking a colleague. Good questions are:
- **Specific**: "What genes are upregulated in both metabolomics and transcriptomics datasets from the AD study?" (better than "Tell me about AD")
- **Contextual**: Include program names, conditions, or omics types
- **Actionable**: Ask for patterns, connections, or explanations

**Example Questions**:

- "Which pathways are enriched in our sphingolipid stress signature?"
- "Do any of our compounds inhibit targets related to neuroinflammation?"
- "What's the correlation between SMPD1 expression and ceramide levels across experiments?"
- "Are there published studies with similar lipidomic profiles to our ALS datasets?"

**Step 3: Review the Answer**

RAG provides:
- **Direct Answer**: Clear, concise response to your question
- **Supporting Evidence**: Links to datasets, experiments, signatures that informed the answer
- **Citations**: References to specific data points
- **Confidence Indicator**: High, medium, or low (based on data availability)

**Step 4: Follow-Up Questions**

RAG maintains context. Ask follow-ups:
- "Show me the top 10 features"
- "Which experiments showed this pattern?"
- "Are there any conflicting results?"

**Step 5: Save Useful Queries**

Click "Save Query" to bookmark questions you ask frequently. You can re-run them later as new data is added.

### Cross-Omics Pathway Analysis

Pathways are networks of genes, proteins, and metabolites that work together to carry out biological functions. Cross-omics pathway analysis reveals which pathways are affected across multiple data types.

#### Running Pathway Enrichment

**Step 1: Select Datasets**

Navigate to "Analysis" → "Pathway Analysis". Select datasets to analyze:
- Single dataset (within-omics)
- Multiple datasets from same experiment (time series)
- Multiple datasets from different omics types (cross-omics integration)

**Step 2: Choose Pathway Databases**

Select which pathway collections to query:
- **KEGG**: Metabolic and signaling pathways
- **Reactome**: Detailed reaction networks
- **GO Biological Process**: Broad functional categories
- **WikiPathways**: Community-curated pathways

**Step 3: Set Parameters**

- **Feature Selection**: Use only significant features (p < 0.05) or all features
- **Statistical Test**: Hypergeometric (over-representation) or GSEA (enrichment)
- **FDR Correction**: Multiple testing correction (Benjamini-Hochberg recommended)

**Step 4: Run Analysis**

Click "Analyze". The system:
1. Maps features to pathways
2. Calculates enrichment for each pathway
3. Applies statistical corrections
4. Ranks pathways by significance

**Step 5: Review Results**

**Pathway Table**:

| Pathway Name | Pathway ID | P-Value | FDR | Features Matched | Direction |
|--------------|------------|---------|-----|------------------|-----------|
| Sphingolipid metabolism | KEGG:00600 | 1.2e-08 | 1.5e-06 | 12/45 | Mixed |
| Apoptosis | KEGG:04210 | 0.003 | 0.02 | 8/87 | Up |

**Visualizations**:

**Enrichment Bar Chart**: Top pathways by significance

**Network View** (if Cytoscape.js is enabled): Interactive graph showing:
- Pathways as nodes (size = number of matched features)
- Connections between pathways (shared features)
- Color-coded by direction (red = upregulated, blue = downregulated)

**Feature-Level Details**:

Click on any pathway to see:
- Which features from your data are in this pathway
- Their fold changes and p-values
- Diagram of pathway with your features highlighted

#### Cross-Omics Integration

When analyzing multiple omics types together, Amprenta performs **multi-layer pathway analysis**.

**Example**:

You have:
- Transcriptomics dataset (gene expression changes)
- Metabolomics dataset (metabolite abundance changes)

Run cross-omics pathway analysis:

**Step 1**: Select both datasets

**Step 2**: System maps:
- Genes → Enzymes/Proteins in pathways
- Metabolites → Substrates/Products in pathways

**Step 3**: Results show pathways affected at **both** levels:

> "Sphingolipid metabolism pathway is enriched:
> - Gene level: SMPD1 (up 2.3-fold), CERS2 (up 1.8-fold)
> - Metabolite level: Ceramide(d18:1/16:0) (up 2.1-fold), Sphingosine (up 1.5-fold)
>
> This indicates coordinated dysregulation at transcriptional and metabolic levels."

**Interpretation**: When both gene expression and metabolite levels change together, it's strong evidence that the pathway is truly perturbed (not just transcriptional noise without functional consequences).

### Mechanism of Action (MOA) Inference

For drug discovery projects, MOA inference predicts how compounds work based on multi-omics signatures.

#### How MOA Inference Works

When you treat cells with a compound, it triggers multi-omics changes. By comparing these changes to known signatures, Amprenta can predict:
- Which targets the compound affects
- Which pathways are perturbed
- What the mechanism of action likely is

#### Running MOA Inference

**Step 1: Prepare Data**

You need:
- **HTS data**: Compound bioactivity (IC50, EC50, binding affinity)
- **Multi-omics data**: Transcriptomics, proteomics, metabolomics from treated samples

**Step 2: Navigate to MOA Inference**

Click "Analysis" → "MOA Inference".

**Step 3: Select Input Data**

- **Compound**: Choose the compound to analyze
- **HTS Campaign**: Select bioactivity data
- **Omics Datasets**: Select transcriptomics/proteomics/metabolomics datasets from treated samples

**Step 4: Configure Analysis**

- **Candidate Targets**: Let system propose targets (automatic) or specify known targets
- **Signature Library**: Compare against your signature library and/or public databases
- **Evidence Features**:
  - Omics concordance (do signatures match?)
  - Bioactivity support (does compound bind relevant targets?)
  - Pathway enrichment (affected pathways make sense?)
  - Literature support (are there publications linking compound to targets?)

**Step 5: Run Inference**

The system builds a probabilistic model:
1. Generates candidate MOAs (target/pathway hypotheses)
2. Scores each candidate based on evidence features
3. Ranks by probability and confidence intervals

This takes 1-2 minutes.

**Step 6: Review MOA Predictions**

**Top Candidates Table**:

| Rank | MOA Candidate | Probability | Confidence Interval | Evidence Score |
|------|---------------|-------------|---------------------|----------------|
| 1 | SMPD1 inhibition | 0.82 | 0.68-0.91 | High |
| 2 | Ceramide synthase modulation | 0.64 | 0.48-0.79 | Medium |
| 3 | Sphingosine kinase activation | 0.23 | 0.10-0.42 | Low |

**Evidence Breakdown** (click on any candidate):

**For SMPD1 inhibition**:
- **Omics concordance**: 85% (43/50 signature features match)
- **Bioactivity support**: IC50 = 120 nM against SMPD1 enzyme
- **Pathway enrichment**: Sphingolipid metabolism pathway highly enriched (p=1e-08)
- **Network proximity**: Compound targets closely connected to SMPD1 in protein network
- **Literature score**: 12 publications link this chemotype to SMPD1

**Visual Evidence**:

**Contribution Chart**: Shows which evidence features contributed most to the probability

**Pathway Overlay**: Highlights affected pathways with evidence from each omics layer

**Confidence Assessment**: Based on evidence consistency:
- **High confidence**: Multiple independent evidence types agree
- **Medium confidence**: Some conflicting evidence, but main hypothesis supported
- **Low confidence**: Weak or contradictory evidence

#### Using MOA Predictions

**Hypothesis Generation**: Use top candidates to design follow-up experiments

**Target Validation**: Prioritize candidates with high bioactivity support for experimental validation

**Mechanism Understanding**: Even if primary target is known, MOA inference reveals off-target effects and downstream pathway perturbations

**Example Workflow**:

> "Compound X shows anti-neurodegenerative activity in phenotypic screen. MOA inference predicts SMPD1 inhibition (82% probability) with high confidence. Design biochemical assay to directly test SMPD1 inhibition. If confirmed, pursue SMPD1 as therapeutic target for neurodegeneration."

### One-Click Narrative Reports

Generate publication-ready reports automatically.

**Step 1: Select Content**

Navigate to an experiment, dataset, or signature. Click "Generate Report".

**Step 2: Choose Report Template**

- **Executive Summary**: High-level findings for PIs or stakeholders
- **Methods & Results**: Detailed technical report for papers
- **Quality Assessment**: Data quality and validation report
- **Comparative Analysis**: Side-by-side comparison of datasets or signatures

**Step 3: Configure Report**

- **Sections to Include**: Check boxes for what to include (methods, results, visualizations, statistics)
- **Visualizations**: Which plots to auto-generate
- **Detail Level**: Summary, standard, or comprehensive

**Step 4: Generate**

Click "Generate Report". The system:
1. Gathers data
2. Runs statistics
3. Creates plots
4. Compiles narrative
5. Formats as PDF or HTML

**Step 5: Download or Share**

Download the report or get a shareable link.

**Example Report Structure**:

```
Dataset Analysis Report: ALS Patient Plasma Metabolomics

Executive Summary
- 127 metabolites detected
- 42 significantly altered (p < 0.05)
- Top pathway: Sphingolipid metabolism (p=1.2e-08)

Methods
- Sample: Plasma from 15 ALS patients, 15 controls
- Platform: Agilent 6550 Q-TOF LC-MS
- Processing: XCMS with manual curation

Results
- [Volcano plot]
- [Top 20 metabolites table]
- [Pathway enrichment chart]

Key Findings
- Ceramides elevated 1.5-2.5 fold in ALS patients
- Sphingosine-1-phosphate reduced 40%
- Pattern matches known neuroinflammation signatures

Conclusions
- Strong evidence of ceramide accumulation in ALS
- Consistent with sphingolipid dysregulation hypothesis
```

### Data Quality Watcher

Automated monitoring ensures your data stays reliable.

#### What the Quality Watcher Monitors

**Data Integrity**:
- Missing values increasing?
- New datasets with unusual distributions?
- Features with suspicious p-values (e.g., p=0.000000)?

**Consistency**:
- Do replicate datasets correlate?
- Are batch effects present?
- Do related datasets show expected patterns?

**Completeness**:
- Are metadata fields filled in?
- Are protocols documented?
- Are external links valid?

#### Setting Up Quality Monitoring

**Step 1: Navigate to Quality Watcher**

Click "Admin" → "Data Quality".

**Step 2: Configure Rules**

Define quality rules:
- **Missing Value Threshold**: Warn if >20% missing values
- **Outlier Detection**: Flag features >3 SD from mean
- **Replicate Correlation**: Expect r > 0.85 between replicates
- **Metadata Completeness**: Require protocol documentation

**Step 3: Set Alert Preferences**

Choose how to be notified:
- **Email**: Immediate alerts for critical issues
- **Dashboard Badge**: Notification icon with count
- **Weekly Digest**: Summary of all quality issues

**Step 4: Activate Monitoring**

Click "Start Monitoring". The system checks data daily and alerts you to issues.

#### Responding to Quality Alerts

When an alert fires:

**Step 1: Review Alert Details**

Click the notification to see:
- Which dataset triggered the alert
- What rule was violated
- Severity (critical, warning, info)

**Step 2: Investigate**

Open the problematic dataset. Check:
- Data preview (are values reasonable?)
- Distribution plots (outliers? bimodal?)
- Metadata (missing info?)

**Step 3: Take Action**

- **Re-import**: If it was a parsing error
- **Flag for Review**: Mark dataset as "Needs Review"
- **Document Issue**: Add note explaining the anomaly
- **Fix Metadata**: Fill in missing fields
- **Dismiss**: If it's a false positive

**Step 4: Mark as Resolved**

Once addressed, mark the alert as resolved. This removes it from your active alerts list and logs the resolution for audit purposes.

### Protocol Version Control & Deviation Tracking

Reproducibility requires knowing exactly what protocol was used. Amprenta tracks protocol versions and flags deviations.

#### Creating a Protocol

**Step 1: Navigate to Protocols**

Click "Data" → "Protocols".

**Step 2: Click "New Protocol"**

**Step 3: Fill in Protocol Details**

- **Name**: Descriptive name (e.g., "LC-MS Lipidomics - Brain Tissue")
- **Version**: 1.0 (increment for updates)
- **Category**: Sample Prep, Instrument Method, Data Processing, etc.
- **Content**: Full protocol text (markdown supported)
  - Materials
  - Step-by-step instructions
  - Instrument settings
  - Quality control steps

**Step 4: Save Protocol**

**Step 5: Assign to Experiments**

When creating experiments, select the protocol used. This links experiments to specific protocol versions.

#### Tracking Protocol Changes

When you update a protocol:

**Step 1: Open Existing Protocol**

**Step 2: Click "Create New Version"**

**Step 3: Make Changes**

Edit the protocol text.

**Step 4: Document Changes**

In the "Change Log" field, explain what changed:
- "Increased extraction time from 30 min to 45 min"
- "Changed column from C18 to HILIC"

**Step 5: Save as New Version**

The system creates version 1.1 (or 2.0 for major changes). The old version is preserved.

#### Deviation Tracking

Sometimes you deviate from the standard protocol. Document these:

**Step 1: Open Experiment**

**Step 2: Click "Report Deviation"**

**Step 3: Describe Deviation**

- **Step Affected**: Which protocol step deviated?
- **Nature of Deviation**: What was different?
- **Reason**: Why did it happen? (equipment failure, sample constraints, etc.)
- **Impact**: Did it affect results?

**Step 4: Save**

The deviation is logged and flagged for review.

#### Protocol Diff Tool

Compare two protocol versions side-by-side to see exactly what changed.

**Step 1: Open a Protocol**

**Step 2: Click "Compare Versions"**

**Step 3: Select Two Versions**

Choose v1.0 and v1.1 (or any two versions).

**Step 4: View Diff**

The system highlights:
- **Added text** in green
- **Removed text** in red
- **Modified text** in yellow

This makes it easy to see how protocols evolved over time.

---

## 6. Chemistry & High-Throughput Screening

### Compound Registration

Before running assays, register compounds in Amprenta's compound library.

#### Registering a Single Compound

**Step 1: Navigate to Compounds**

Click "Chemistry" → "Compound Library".

**Step 2: Click "Register Compound"**

**Step 3: Enter Compound Information**

- **SMILES**: Structural representation (required)
  Example: `CCO` for ethanol
- **Compound Name**: Common name
- **Molecular Formula**: Auto-calculated from SMILES
- **Molecular Weight**: Auto-calculated
- **Corporate ID**: Your internal tracking number
- **Batch ID**: Specific batch if multiple syntheses
- **Purity**: % purity (if known)
- **Vendor**: Supplier name
- **Notes**: Synthesis details, storage conditions, etc.

**Step 4: Calculate Properties**

Click "Calculate Properties". The system computes:
- LogP (lipophilicity)
- TPSA (topological polar surface area)
- H-Bond donors/acceptors
- Rotatable bonds
- Lipinski's Rule of Five compliance

**Step 5: Structure Validation**

The system checks:
- Is SMILES valid?
- Are there obvious structural errors (valence violations)?
- Is this compound already registered? (de-duplication check)

**Step 6: Save**

Click "Register". The compound is added to your library with a unique internal ID (e.g., `CMP-00001`).

#### Batch Registration

For large libraries, use batch registration.

**Prepare CSV File**:

```csv
smiles,compound_name,corporate_id,purity,vendor
CCO,Ethanol,CORP-001,99.5,Sigma-Aldrich
CC(=O)O,Acetic acid,CORP-002,99.9,Fisher Scientific
c1ccccc1,Benzene,CORP-003,99.0,Sigma-Aldrich
```

**Upload File**:

From "Chemistry" → "Compound Library" → "Batch Register", upload your CSV. The system:
1. Validates all SMILES
2. Calculates properties
3. Checks for duplicates
4. Registers valid compounds
5. Reports errors (invalid SMILES, duplicates)

You can download an error report and fix issues before re-uploading.

### Structure Search

Find compounds by structural similarity or substructure.

#### Substructure Search

Find all compounds containing a specific substructure (e.g., all compounds with a benzene ring).

**Step 1: Navigate to Structure Search**

Click "Chemistry" → "Structure Search".

**Step 2: Draw or Enter Substructure**

- **Option 1**: Draw using the molecular editor
- **Option 2**: Enter SMARTS pattern (e.g., `c1ccccc1` for benzene)

**Step 3: Click "Search"**

The system finds all compounds containing that substructure.

**Step 4: Review Results**

Results show:
- Compound structure
- Compound name
- Corporate ID
- Properties
- Links to associated HTS data (if any)

#### Similarity Search

Find compounds structurally similar to a query compound.

**Step 1: Enter Query Compound**

Draw structure or paste SMILES.

**Step 2: Set Similarity Threshold**

- **Tanimoto Coefficient**: 0.0-1.0 (0.7 = 70% similar)
- Higher thresholds = more similar compounds (fewer results)

**Step 3: Search**

The system calculates fingerprints and compares to your library.

**Step 4: Review Results**

Results ranked by similarity:

| Compound | Similarity | Structure | Corporate ID |
|----------|------------|-----------|--------------|
| Compound A | 0.95 | [structure] | CORP-123 |
| Compound B | 0.88 | [structure] | CORP-456 |

Click on any compound to see full details and associated data.

### HTS Campaign Management

High-throughput screening generates massive amounts of data. Amprenta organizes it.

#### Creating an HTS Campaign

**Step 1: Navigate to HTS**

Click "Chemistry" → "HTS Campaigns".

**Step 2: Click "New Campaign"**

**Step 3: Campaign Details**

- **Campaign Name**: Descriptive (e.g., "SMPD1 Inhibitor Screen Q4 2025")
- **Target**: Biological target (protein, pathway, phenotype)
- **Assay Type**: Biochemical, cell-based, phenotypic
- **Screen Type**: Primary screen, dose-response, counterscreen
- **Date Started**: Campaign start date
- **Status**: Active, Completed, On Hold

**Step 4: Upload Plate Layout**

If using plate-based screening, upload plate maps:
- Which wells contain compounds
- Which wells are controls (positive, negative)
- Replicate structure

**Example CSV**:

```csv
well,compound_id,concentration,control_type
A1,CMP-00001,10uM,
A2,CMP-00002,10uM,
A3,DMSO,,negative_control
H12,REFERENCE,,positive_control
```

**Step 5: Save Campaign**

The campaign is created and ready to receive data.

#### Uploading HTS Results

**Step 1: Open Campaign**

Click on the campaign from the list.

**Step 2: Click "Upload Results"**

**Step 3: Select Results File**

Upload CSV with columns:
- `well` or `compound_id`
- `readout` (e.g., % inhibition, fluorescence, viability)
- `raw_value` (optional)
- `normalized_value` (optional)
- `z_score` (optional)

**Step 4: Map Columns**

System auto-detects columns but verify:
- Which column is the primary readout?
- Are controls labeled correctly?

**Step 5: Import**

The system:
1. Links compounds to results
2. Calculates statistics (Z', CV)
3. Flags hits (based on thresholds you set)
4. Generates plate heatmaps

**Step 6: Review Quality Metrics**

After import, view:
- **Z' Factor**: Measure of assay quality (-∞ to 1, >0.5 is good)
- **Hit Rate**: Percentage of compounds flagged as hits
- **Control Performance**: Are positive/negative controls behaving as expected?

### HTS QC & Triage Assistant

Automated analysis helps you quickly identify hits and prioritize follow-ups.

#### Running QC Analysis

**Step 1: Open HTS Campaign**

**Step 2: Click "Run QC Analysis"**

The system checks:
- **Plate Effects**: Are edge wells different? (common artifact)
- **Control Consistency**: Do controls have low variability?
- **Hit Distribution**: Are hits evenly distributed or clustered? (clustering suggests artifacts)
- **Outlier Detection**: Flagging wells with extreme values

**Step 3: Review QC Report**

**Plate Heatmap**: Visual representation of entire plate with color-coding (red=high activity, blue=low activity)

**Edge Effect Analysis**: Comparison of edge vs. center wells

**Control Statistics**:
- Positive control: Mean, SD, CV
- Negative control: Mean, SD, CV

**Z' Factor**: Overall assay quality metric

**Recommendations**: System suggests actions (e.g., "Edge effects detected; consider excluding edge wells")

#### Hit Triage

After identifying hits, triage them for follow-up.

**Step 1: Set Hit Criteria**

Define what constitutes a "hit":
- **Activity Threshold**: e.g., >50% inhibition
- **Statistical Threshold**: e.g., >3 SD from negative control
- **Reproducibility**: Must appear in 2/3 replicates

**Step 2: Apply Criteria**

System flags compounds meeting criteria.

**Step 3: Traffic Light Scoring**

Compounds are scored in multiple dimensions:
- **Potency**: Red (high), Yellow (medium), Green (low)
- **Drug-likeness**: Green (Lipinski-compliant), Yellow (1 violation), Red (2+ violations)
- **Novelty**: Green (new chemotype), Yellow (known scaffold), Red (frequent hitter)
- **Selectivity**: Green (selective), Yellow (moderate off-targets), Red (promiscuous)

**Example**:

| Compound | Potency | Drug-like | Novelty | Selectivity | Priority |
|----------|---------|-----------|---------|-------------|----------|
| CMP-123 | 🔴 High | 🟢 Yes | 🟢 Novel | 🟢 Selective | 1 |
| CMP-456 | 🔴 High | 🟡 1 viol. | 🟡 Known | 🟢 Selective | 2 |
| CMP-789 | 🟡 Med | 🟢 Yes | 🟢 Novel | 🔴 Promiscuous | 5 |

**Step 4: Prioritize Follow-Ups**

Sort by priority score. Focus on compounds with mostly green/yellow lights.

**Step 5: Export Hit List**

Download prioritized hit list for dose-response testing or SAR expansion.

### Structure-Activity Relationship (SAR) Analysis

Understanding how structure relates to activity guides optimization.

#### Viewing SAR Data

**Step 1: Navigate to SAR Analysis**

Click "Chemistry" → "SAR Analysis".

**Step 2: Select Target or Campaign**

Choose the assay you want to analyze.

**Step 3: View Activity Summary**

**Scatter Plot**: Structure similarity (x-axis) vs. activity (y-axis)

Clusters indicate:
- **Activity cliffs**: Structurally similar compounds with very different activity (interesting for SAR)
- **Flat regions**: Structural changes don't affect activity (not useful for optimization)

#### Activity Cliff Detection

Activity cliffs are pairs of structurally similar compounds with large activity differences.

**Step 1: Click "Find Activity Cliffs"**

**Step 2: Set Parameters**

- **Similarity Threshold**: How similar? (e.g., Tanimoto > 0.85)
- **Activity Difference**: Minimum fold-change (e.g., 5-fold)

**Step 3: Review Cliff Pairs**

System shows pairs:

**Pair 1**: CMP-123 vs. CMP-124
- Similarity: 0.92 (very similar)
- Activity: CMP-123 = 10 nM, CMP-124 = 200 nM (20-fold difference)
- Structural difference: Single methyl → ethyl substitution at R1

**Interpretation**: Small structural change causes big activity loss. This position is critical for binding.

#### R-Group Analysis

If you have a series of analogs, analyze which substituents (R-groups) are best.

**Step 1: Define Core Structure**

Draw or enter the common scaffold shared by your series.

**Step 2: System Decomposes Series**

Amprenta identifies R-groups automatically.

**Step 3: View R-Group Table**

| R1 | R2 | Activity (IC50) | Lipophilicity | TPSA |
|----|----|-----------------|--------------  |------|
| -CH3 | -H | 50 nM | 2.3 | 45 |
| -C2H5 | -H | 200 nM | 2.8 | 45 |
| -CH3 | -Cl | 10 nM | 3.1 | 45 |

**Insights**:
- R1 = methyl is better than ethyl (smaller is better)
- R2 = chloro improves potency 5-fold
- Optimal: R1 = methyl, R2 = chloro → 10 nM potency

**Step 4: Design Next Round**

Based on SAR, propose new analogs to test.

### SAR Delta Explorer (Voila Dashboard)

Interactive dashboard for exploring SAR in real-time.

**Step 1: Open SAR Delta Explorer**

From "Chemistry" → "SAR Analysis", click "Open in Jupyter" → "SAR Delta Explorer".

**Step 2: Select Campaign or Series**

Choose which HTS campaign or compound series to analyze.

**Step 3: Interactive Exploration**

**Features**:
- **3D Structure Viewer**: Rotate and inspect compound structures (py3Dmol)
- **Activity Plot**: Dynamically filter by activity range
- **Matched Molecular Pairs**: Automatically identified pairs differing by single transformation
- **R-Group Grid**: Matrix view of R-group combinations with activity color-coding

**Step 4: Export Insights**

Click "Generate Report" to create a PDF summary of SAR findings.

---

## 7. Jupyter Integration

### Why Jupyter?

While Amprenta's dashboard covers most common tasks, sometimes you need custom analysis, specialized visualizations, or integration with other tools. Jupyter notebooks provide a flexible environment for:

- Custom statistical analyses
- Specialized visualizations (3D plots, network graphs)
- Integration with R or other languages
- Reproducible analysis workflows
- Publication-quality figures

### Opening Jupyter from the Dashboard

The easiest way to start a notebook is directly from Amprenta's interface.

**Step 1: Navigate to Data**

Open any:
- Experiment
- Dataset
- Signature
- Compound
- HTS Campaign

**Step 2: Click "Open in Jupyter"**

A button labeled "Open in Jupyter" appears on most detail pages.

**Step 3: Authenticate**

If this is your first time, you'll authenticate via SSO (single sign-on). The system uses JWT tokens to securely connect Jupyter to your Amprenta session.

**Step 4: Notebook Launches**

A new tab opens with JupyterLab. Your selected data (experiment, dataset, etc.) is **pre-loaded** into the notebook environment.

**Example**:

If you clicked "Open in Jupyter" from a dataset detail page, the notebook automatically contains:

```python
# Pre-loaded by Amprenta
from amprenta_client import AmprentaClient

client = AmprentaClient()
dataset = client.datasets.get("DATASET-12345")

# Your data is ready to use
features = dataset.features
print(f"Loaded {len(features)} features")
```

### Using Notebook Templates

Amprenta provides pre-built templates for common analyses.

#### Available Templates

**Getting Started**:
- Introduction to data access
- Basic plotting
- Filtering and statistics

**Molecule Analysis**:
- Structure visualization
- Property calculation
- Similarity search

**Signature Explorer**:
- Creating signatures
- Scoring datasets
- Visualizing matches

**Pathway Enrichment**:
- Running KEGG/Reactome enrichment
- Plotting pathway networks

**Dose-Response Fitting**:
- IC50/EC50 calculation
- Hill slope fitting
- Confidence intervals

**Compound Similarity Clustering**:
- Chemical space visualization
- Clustering by fingerprint
- Diversity analysis

#### Opening a Template

**Step 1: Navigate to Jupyter Hub**

Click "Jupyter" in the main menu → "Notebook Templates".

**Step 2: Select a Template**

Browse templates by category. Click on one to preview.

**Step 3: Create from Template**

Click "Use This Template". A new notebook is created with:
- Pre-filled code cells
- Comments explaining each step
- Example outputs (so you see what to expect)

**Step 4: Customize**

Edit the notebook:
- Replace example data with your data
- Adjust parameters
- Add your own analysis steps

**Step 5: Run**

Execute cells sequentially (Shift+Enter). The notebook runs your custom analysis.

### Voila Dashboards

Voila turns Jupyter notebooks into interactive web dashboards (no code visible, just widgets and outputs).

#### Available Dashboards

**HTS Plate Viewer**:
- Upload plate data
- Interactive heatmap
- Z' factor calculation
- Control visualization
- Export PNG/CSV

**Compound Triage Dashboard**:
- Import HTS results
- Traffic-light scoring (potency, drug-likeness, novelty)
- Pareto front optimization (balance multiple objectives)
- Hit list export

**SAR Delta Explorer**:
- Matched molecular pair analysis
- R-group table
- 3D structure viewer
- Activity cliff detection

**Signature Validation Console**:
- Score signature against test datasets
- Per-feature contribution
- Approve/reject workflow

**Pathway Impact Explorer**:
- Multi-omics pathway enrichment
- Cytoscape graph overlay
- Cross-omics toggle (show/hide layers)

#### Opening a Dashboard

**Step 1: Navigate to Voila Dashboards**

Click "Jupyter" → "Dashboards".

**Step 2: Select Dashboard**

Click on a dashboard name (e.g., "HTS Plate Viewer").

**Step 3: Dashboard Loads**

A clean interface appears with:
- Input controls (file upload, parameter settings)
- Visualization area
- Export buttons

**Step 4: Interact**

Use widgets to:
- Upload data
- Adjust parameters
- Switch views
- Export results

No coding required—just point-and-click.

### Custom Analysis Workflows

For advanced users, build custom workflows combining multiple tools.

#### Example: End-to-End Lipidomics Analysis

**Goal**: Analyze lipidomics dataset, find significant lipids, run pathway enrichment, generate report.

**Step 1: Create New Notebook**

From JupyterLab, click "New Notebook" → Python 3.

**Step 2: Load Data**

```python
from amprenta_client import AmprentaClient
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

client = AmprentaClient()
dataset = client.datasets.get("your-dataset-id")
features = dataset.features_as_dataframe()
```

**Step 3: Filter Significant Lipids**

```python
significant = features[(features['p_value'] < 0.05) & 
                       (abs(features['log2_fc']) > 1.0)]
print(f"{len(significant)} significant lipids")
```

**Step 4: Visualize**

```python
# Volcano plot
plt.figure(figsize=(10, 6))
plt.scatter(features['log2_fc'], -np.log10(features['p_value']), 
            c='gray', alpha=0.5, label='Not significant')
plt.scatter(significant['log2_fc'], -np.log10(significant['p_value']),
            c='red', alpha=0.8, label='Significant')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-Value')
plt.legend()
plt.title('Volcano Plot: Lipidomics Dataset')
plt.show()
```

**Step 5: Pathway Enrichment**

```python
# Run enrichment via API
enrichment = client.pathways.enrich(
    features=significant['feature_id'].tolist(),
    databases=['KEGG', 'Reactome']
)

# Display top pathways
print(enrichment.head(10))
```

**Step 6: Create Signature**

```python
# Create signature from significant features
signature = client.signatures.create(
    name="ALS Plasma Lipidomics Signature",
    features=significant['feature_id'].tolist(),
    weights=significant['log2_fc'].tolist(),
    directions=significant['direction'].tolist()
)
print(f"Created signature: {signature.id}")
```

**Step 7: Generate Report**

```python
# Use Papermill to parameterize and export
# (Report generation automatically creates PDF)
client.reports.generate(
    notebook_path="current_notebook.ipynb",
    output_format="pdf",
    title="ALS Lipidomics Analysis Report"
)
```

This entire workflow runs in one notebook and can be re-run on new datasets by changing the dataset ID.

### Publishing and Scheduling

Make your notebooks available to others or run them automatically.

#### Publishing a Dashboard

**Step 1: Finish Your Notebook**

Ensure all cells run correctly and outputs look good.

**Step 2: Add Widgets (if interactive)**

Use `ipywidgets` to add controls:

```python
import ipywidgets as widgets

dataset_selector = widgets.Dropdown(
    options=['Dataset A', 'Dataset B', 'Dataset C'],
    description='Dataset:'
)

threshold_slider = widgets.FloatSlider(
    value=0.05,
    min=0.01,
    max=0.1,
    step=0.01,
    description='P-value:'
)
```

**Step 3: Test with Voila**

From the notebook, click "File" → "Open with Voila". This previews how it will look as a dashboard.

**Step 4: Pin to Dashboard Catalog**

Click "Publish to Catalog". Give it:
- Title
- Description
- Tags (for searching)
- Access permissions (private, team, public)

**Step 5: Share Link**

Users can now access your dashboard from "Jupyter" → "Dashboards" → [Your Dashboard Name].

#### Scheduling Automated Reports

Run notebooks on a schedule (e.g., weekly quality reports).

**Step 1: Navigate to Scheduled Jobs**

Click "Jupyter" → "Scheduled Notebooks".

**Step 2: Click "New Scheduled Job"**

**Step 3: Configure Job**

- **Notebook**: Select which notebook to run
- **Parameters**: If using Papermill, specify parameter values
- **Schedule**: Cron expression (e.g., `0 9 * * MON` = 9 AM every Monday)
- **Output Destination**: Where to save results (email, shared folder, dashboard)

**Step 4: Activate**

Click "Activate Schedule". The notebook runs automatically at the specified times.

**Example Use Case**:

> "Every Monday at 9 AM, run the 'Weekly Data Quality Report' notebook. It checks all datasets imported in the past week, flags any quality issues, and emails a summary to the PI."

### Jupyter Best Practices

**Organize Your Notebooks**:
- Create folders by project or analysis type
- Use descriptive names (not "Untitled1.ipynb")
- Add markdown cells explaining what each section does

**Document Your Code**:
- Add comments for complex logic
- Use markdown cells for narrative explanations
- Include references to methods or papers

**Make Notebooks Reproducible**:
- Set random seeds (`np.random.seed(42)`)
- Record versions of packages used
- Include data provenance (dataset IDs, dates)

**Use Version Control**:
- Commit notebooks to Git (if your organization uses it)
- Track major changes

**Separate Data from Code**:
- Don't hard-code file paths
- Use parameters for dataset IDs
- Load data via API (not local files)

---

## 8. Advanced Features

### Concurrent Editing Safety

Multiple users can work on the same data simultaneously without conflicts.

#### How It Works

Amprenta uses **optimistic locking**:
1. When you open a dataset/signature/experiment, the system records the current version
2. You make edits locally
3. When you save, the system checks: has anyone else saved changes since you opened it?
4. If no: your changes save normally
5. If yes: you see a conflict warning

#### Resolving Conflicts

When a conflict occurs:

**Step 1: Conflict Dialog Appears**

You see:
- **Your Changes**: What you were trying to save
- **Current Version**: What's now in the database (someone else's changes)
- **Diff View**: Side-by-side comparison

**Step 2: Choose Resolution Strategy**

- **Keep Your Changes**: Override their changes (use carefully)
- **Merge**: Combine both sets of changes (system assists)
- **Discard Your Changes**: Reload their version

**Step 3: Save**

After resolving, save again.

**Best Practice**: Communicate with your team about who's editing what. Use comments or status fields to signal "I'm working on this."

### Bulk Operations

Perform actions on many items at once.

#### Bulk Tagging

**Step 1: Select Multiple Items**

From a list view (datasets, experiments, etc.), check boxes next to items.

**Step 2: Click "Bulk Actions" → "Add Tags"**

**Step 3: Enter Tags**

Type tags (comma-separated).

**Step 4: Apply**

All selected items are tagged.

#### Bulk Export

Export multiple datasets to CSV in one zip file.

**Step 1: Select Datasets**

Check boxes next to datasets to export.

**Step 2: Click "Bulk Actions" → "Export"**

**Step 3: Choose Format**

- CSV (features only)
- JSON (full metadata)
- Excel (multi-sheet workbook)

**Step 4: Download**

Zip file downloads containing all selected datasets.

#### Bulk Status Changes

Change status of multiple experiments at once (e.g., mark as "Completed").

**Step 1: Select Experiments**

**Step 2: Click "Bulk Actions" → "Change Status"**

**Step 3: Select New Status**

**Step 4: Confirm**

All experiments updated.

### API Access

For programmatic access outside Jupyter, use Amprenta's REST API.

#### Getting an API Key

**Step 1: Navigate to Settings**

Click your username → "Settings" → "API Keys".

**Step 2: Click "Generate New Key"**

**Step 3: Name Your Key**

Give it a descriptive name (e.g., "Analysis Pipeline Integration").

**Step 4: Set Permissions**

Choose what this key can do:
- Read-only
- Read + write datasets
- Full access

**Step 5: Copy Key**

The key is shown **once**. Copy and store it securely (e.g., in a password manager).

#### Using the API

**Base URL**: `http://your-amprenta-instance.com/api/v1`

**Authentication**: Include API key in header:
```
Authorization: Bearer YOUR_API_KEY
```

**Example (Python)**:

```python
import requests

headers = {"Authorization": "Bearer YOUR_API_KEY"}
response = requests.get(
    "http://your-amprenta-instance.com/api/v1/datasets",
    headers=headers
)

datasets = response.json()
print(f"Found {len(datasets)} datasets")
```

**Common Endpoints**:

- `GET /programs` - List programs
- `GET /experiments/{id}` - Get experiment details
- `GET /datasets/{id}/features` - Get dataset features
- `POST /signatures` - Create signature
- `POST /datasets/{id}/score` - Score dataset against signature
- `GET /compounds` - List compounds
- `GET /campaigns/{id}/results` - Get HTS results

See full API documentation at `/api/docs` (interactive Swagger UI).

### Admin Features

If you have admin privileges, you can manage users, permissions, and system settings.

#### User Management

**Step 1: Navigate to Admin Panel**

Click "Admin" → "Users".

**Step 2: View User List**

See all users, their roles, and last login.

**Step 3: Add New User**

Click "Add User". Fill in:
- Email
- Name
- Role (Admin, Researcher, Viewer)
- Initial password (they change on first login)

**Step 4: Manage Roles**

Click on a user to edit:
- Change role
- Activate/deactivate account
- Reset password

#### System Settings

**Step 1: Navigate to System Settings**

Click "Admin" → "Settings".

**Step 2: Configure Options**

- **Authentication**: SSO, local auth, LDAP
- **Storage**: Local filesystem, S3, etc.
- **API Rate Limits**: Requests per minute
- **Email Notifications**: SMTP settings
- **Backup Schedule**: Automated backups

**Step 3: Save Changes**

Changes take effect immediately (some may require restart).

### Backup and Export

Protect your data with regular backups.

#### Manual Backup

**Step 1: Navigate to Admin Panel**

Click "Admin" → "Backup".

**Step 2: Click "Create Backup"**

**Step 3: Select What to Backup**

- Full system (database + files)
- Database only
- Specific programs

**Step 4: Click "Start Backup"**

The system creates a backup archive. This may take several minutes for large datasets.

**Step 5: Download**

Once complete, download the backup file. Store it securely offsite.

#### Automatic Backups

**Step 1: Configure Backup Schedule**

In "Admin" → "Settings" → "Backup", set:
- Frequency (daily, weekly)
- Time (e.g., 2 AM)
- Retention (keep last N backups)

**Step 2: Set Destination**

Choose where backups are stored:
- Local directory
- S3 bucket
- Network share

**Step 3: Activate**

Click "Enable Automatic Backups". System handles it from then on.

### Audit Logs

Track all actions for compliance and troubleshooting.

**Step 1: Navigate to Audit Logs**

Click "Admin" → "Audit Logs".

**Step 2: Filter Logs**

Filter by:
- User
- Action type (create, update, delete)
- Resource type (dataset, signature, etc.)
- Date range

**Step 3: View Details**

Click on any log entry to see:
- Who performed the action
- What was changed
- When it happened
- IP address
- Before/after values (for edits)

**Step 4: Export**

Export audit logs for external analysis or archiving.

---

## 9. Troubleshooting & FAQs

### Common Issues

#### Issue: Cannot Upload Dataset

**Symptoms**: Upload fails with "Invalid file format" error.

**Solutions**:
1. **Check File Format**: Only CSV, TSV, and Excel (.xlsx) are supported. Convert other formats first.
2. **Check Encoding**: Files must be UTF-8 encoded. Open in a text editor and save as UTF-8.
3. **Check Column Names**: Avoid special characters (*, /, \, etc.) in column headers.
4. **Check File Size**: Files over 500 MB may time out. Split into smaller files.

#### Issue: Features Not Mapping to Standard IDs

**Symptoms**: After import, many features show "Unknown" or aren't normalized.

**Solutions**:
1. **Use Standard Nomenclature**: For genes, use HGNC symbols (not Ensembl IDs). For metabolites, use HMDB IDs.
2. **Check Spelling**: Typos prevent matching (e.g., "SPMD1" instead of "SMPD1").
3. **Update ID Mappings**: If using an older database, ID mappings may be outdated. Contact admin to refresh.

#### Issue: RAG Search Returns "No Results"

**Symptoms**: You ask a question but get no answer.

**Solutions**:
1. **Rephrase Question**: Try different wording. Instead of "Show me ALS data", try "Which datasets are from ALS patients?"
2. **Check Data Availability**: If you haven't imported relevant data yet, RAG can't find it.
3. **Check Tags**: Are your datasets tagged correctly? RAG uses tags for context.

#### Issue: Signature Score is 0 or Very Low

**Symptoms**: You score a dataset against a signature and get a near-zero score.

**Solutions**:
1. **Check Omics Type Match**: Don't score a proteomics signature against a metabolomics dataset (unless it's a multi-omics signature).
2. **Check Feature Overlap**: View details to see how many signature features appear in the dataset. If very few, the dataset may not be relevant.
3. **Adjust Thresholds**: Lower the p-value or fold change thresholds to allow more features to match.

#### Issue: Jupyter Notebook Won't Open

**Symptoms**: Clicking "Open in Jupyter" does nothing or shows an error.

**Solutions**:
1. **Check Authentication**: You may need to re-authenticate. Log out and back in.
2. **Check JupyterHub Status**: Ask admin if JupyterHub is running.
3. **Clear Browser Cache**: Sometimes cached credentials cause issues. Try incognito mode.

#### Issue: HTS Import Fails

**Symptoms**: HTS results upload fails or imports incomplete data.

**Solutions**:
1. **Check Column Mapping**: Ensure the `compound_id` or `well` column matches your compound library.
2. **Check Data Types**: Activity values should be numeric. Remove any text annotations.
3. **Check Controls**: Control labels (positive_control, negative_control) must match exactly (case-sensitive).

### Frequently Asked Questions

#### Q: Can I analyze data from multiple species?

**A**: Yes. When creating experiments, specify the species. Amprenta handles human, mouse, rat, and other model organisms. Cross-species comparisons use ortholog mapping.

#### Q: How do I share data with external collaborators?

**A**: Two options:
1. **Export**: Download datasets as CSV/Excel and share via email or file transfer
2. **Guest Access** (if enabled by admin): Admin can create guest accounts with limited permissions for specific programs

#### Q: What happens if I delete a dataset by mistake?

**A**: Datasets are **soft-deleted**—they're hidden but not permanently removed. Contact your admin within 30 days to restore. After 30 days, they're permanently deleted.

#### Q: Can I import data from other databases (e.g., ChEMBL, PubChem)?

**A**: Not directly via the interface, but:
1. Download data from those databases
2. Format as CSV
3. Upload to Amprenta as a dataset

For programmatic access, use the API to integrate with external databases.

#### Q: How do I cite Amprenta in a publication?

**A**: Use this citation (update with your version):

> "Data analysis was performed using Amprenta RAG (version 1.0, 2025), a multi-omics knowledge management platform with RAG capabilities."

Include a link to your organization's Amprenta instance or documentation.

#### Q: Can I use Amprenta offline?

**A**: The main system requires internet access. However:
- You can export data and analyze offline
- Jupyter notebooks can run locally (after exporting from Amprenta)
- For fully offline deployments, contact your admin about local hosting

#### Q: How often should I run data quality checks?

**A**: Recommended:
- **After every import**: Review quality score
- **Weekly**: Run Quality Watcher summary
- **Before publishing**: Full quality audit of datasets used in results

#### Q: What's the maximum dataset size?

**A**: Current limits (may vary by deployment):
- **Single dataset**: 100,000 features
- **File size**: 500 MB
- **Total program data**: No hard limit (depends on storage)

For larger datasets, contact admin about optimizations or batch processing.

#### Q: Can I use Amprenta for clinical data?

**A**: Amprenta is designed for research data. For clinical use:
1. Ensure your deployment meets regulatory requirements (HIPAA, GDPR, etc.)
2. Enable audit logging
3. Restrict access with role-based permissions
4. Consult with your compliance team

#### Q: How do I report a bug or request a feature?

**A**: Contact your system administrator or use the feedback form in the app (Help menu → "Send Feedback").

#### Q: Is my data private?

**A**: Yes. Data visibility is controlled by:
1. **Program-level permissions**: Only users assigned to a program can see its data
2. **Role-based access**: Viewers can't edit; only admins can delete
3. **Audit logs**: All accesses are logged

For multi-tenant deployments, data is fully segregated by organization.

---

## Conclusion

Congratulations on completing the Amprenta RAG User Guide! You now have the knowledge to:

- Organize research data with programs, experiments, and datasets
- Import multi-omics data and public repository datasets
- Create and score signatures for pattern discovery
- Use RAG to ask intelligent questions and get contextualized answers
- Analyze cross-omics pathways and infer mechanisms of action
- Manage compound libraries and HTS campaigns
- Perform SAR analysis
- Create custom analyses with Jupyter notebooks
- Build interactive dashboards with Voila
- Use advanced features for collaboration and reproducibility

### Next Steps

**For New Users**:
1. Create your first program
2. Upload a test dataset
3. Create a signature from that dataset
4. Ask a question with RAG
5. Explore a Voila dashboard

**For Advanced Users**:
1. Set up automated reports with scheduled notebooks
2. Build custom analysis pipelines with the API
3. Create organization-wide signature libraries
4. Integrate Amprenta with your existing tools

**For Administrators**:
1. Configure user roles and permissions
2. Set up automated backups
3. Enable SSO for your organization
4. Monitor system health and audit logs

### Getting Help

- **In-App Help**: Click the "?" icon in the top-right corner of any page
- **Documentation**: Visit the docs portal (link in Help menu)
- **Support**: Contact your system administrator
- **Community**: Join user forums (if available at your organization)

### Stay Updated

Amprenta is continuously improving. Check the **changelog** (Help → "What's New") to see new features and enhancements.

---

**Thank you for using Amprenta RAG!**

We're excited to see the discoveries you'll make with this platform. Science is better when knowledge is organized, accessible, and actionable—Amprenta makes that possible.

Happy researching! 🧬🔬
