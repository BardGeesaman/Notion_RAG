# Amprenta Platform - User Guide

Welcome to the Amprenta Multi-Omics RAG System! This comprehensive guide will help you navigate the platform, manage your data, and leverage advanced analytics to drive scientific discovery.

---

## 🚀 1. Getting Started

### Access & Login
*   **URL**: Access the platform via your deployed URL or `localhost:8501` for local instances.
*   **Login**: Use your assigned credentials. If you don't have an account, contact your administrator (or use the "Register" tab if enabled).
*   **Permissions**: Your access to specific features (e.g., Admin tools, Chemistry) is determined by your user role.
    *   **Admin**: Full access to all settings, user management, and logs.
    *   **Researcher**: Can create experiments, register compounds, and run analyses.
    *   **Viewer**: Read-only access to published data and dashboards.

### Dashboard Overview
The **Overview** page is your command center.
*   **Metrics**: See high-level stats (Total Experiments, Active Compounds, Datasets).
*   **Activity Feed**: A real-time timeline of actions taken by your team (e.g., "John registered compound AMP-123", "Sarah imported Dataset X").
*   **Widgets**: Customizable cards showing recent data or alerts (e.g., "3 New Draft Studies Found").
*   **Global Search**: Use the search bar (or `Cmd+K` / `Ctrl+K`) to jump to any entity (Gene, Compound, Experiment) instantly.

### Navigation
The sidebar is organized into functional groups:
*   **Data Management**: Core entities (Programs, Experiments, Datasets).
*   **Chemistry**: Compound registration, screening data, and SAR tools.
*   **Discovery**: Automated workflows, repository harvesting, and signature management.
*   **Analysis**: Statistical tools, pathway enrichment, and interactive visualizations.
*   **Knowledge**: RAG Query interface, Q&A Tracker, and Literature analysis.
*   **Collaboration**: Teams, Lab Notebook, and Scheduling tools.
*   **Admin**: System health monitoring, audit logs, and global settings.

---

## 🧪 2. Data Management (The Foundation)

### Programs & Experiments
Programs and Experiments organize your research.
*   **Programs**: Top-level research areas (e.g., "Neurodegeneration", "Oncology"). View all related experiments and datasets here.
    *   *Usage*: Create a Program first to act as a container for your work.
*   **Experiments**: Specific studies or campaigns.
    *   **Protocols**: Link standard operating procedures (SOPs) to experiments for reproducibility.
    *   **Metadata**: Capture critical context like `Cell Line`, `Treatment`, `Dose`, and `Timepoint`.

### Datasets & Samples
*   **Data Ingestion**: Use the **Import Data** page to upload raw files.
    *   **Supported Formats**: CSV, TSV, Excel.
    *   **Omics Types**: Transcriptomics (RNA-seq), Proteomics, Metabolomics, Lipidomics.
    *   **Validation**: The system checks for required columns (e.g., `gene_symbol`, `log2fc`, `p_value`) and data types automatically.
*   **Sample Inventory**: Track physical samples.
    *   **Location**: Specify Freezer, Rack, Box, and Position (e.g., "Freezer 1 > Shelf 2 > Box A > A1").
    *   **Lineage**: Track parent-child relationships (e.g., "Tissue Sample A" -> "RNA Aliquot A1").
*   **Datasets**: Browse imported data.
    *   **Feature View**: Click a dataset to see a table of all features (genes/metabolites) with search and sort capabilities.
    *   **Summary Stats**: View distributions of Log2FC and P-values to assess data quality.

---

## ⚗️ 3. Chemistry & Lead Optimization

### Compound Management
*   **Registration**: Use the **Chemistry > Compounds** page to register new molecules.
    *   **Input**: Enter a SMILES string directly or use the **Sketcher** to draw the structure.
    *   **Duplicate Check**: The system calculates a hash (InChIKey) to prevent registering the same molecule twice.
    *   **Properties**: Automatically calculates Molecular Weight (MW), LogP, TPSA, and Hydrogen Bond Donors/Acceptors.
    *   **ID**: Assigns a unique corporate ID (e.g., `AMP-00123`).

### Screening & SAR
*   **Assay Results**: View biochemical and cellular assay data linked to compounds.
    *   *Example*: IC50 values, % Inhibition, Cell Viability.
*   **SAR Analysis**:
    *   **R-Group Decomposition**: Select a core scaffold and see how different R-groups affect potency.
    *   **Activity Cliffs**: The system identifies pairs of similar molecules (high Tanimoto similarity) with drastically different activity (e.g., >10x difference in IC50), highlighting critical structural features.
    *   **Pharmacophore Search**: define 3D constraints (e.g., "Must have an aromatic ring here and a hydrogen donor there") to find matching compounds.

### Lead Optimization
*   **Candidate Selection**: Track compounds through the discovery funnel (Hit -> Lead -> Candidate).
*   **TPP Scoring**: Define a **Target Product Profile (TPP)** with criteria (e.g., "MW < 500", "IC50 < 10nM").
    *   **Traffic Light Scoring**: The system compares each compound against the TPP and assigns Green/Yellow/Red status to each property.
*   **Vendor Search**: Search external catalogs (MolPort, Mcule, Enamine) to find commercially available analogs for purchase.

---

## 🧬 4. Omics Analysis & Discovery

### Signatures & Patterns
*   **Signatures**: Collections of features (genes, proteins, metabolites) that change together in a biological context.
    *   **Discovery**: The system automatically identifies signatures using pattern mining algorithms (e.g., "Find genes that are always upregulated in ALS samples").
    *   **Scoring**: You can "score" a new dataset against established signatures to see if it matches known disease patterns.
*   **Cross-Omics Analysis**: Find convergence points.
    *   *Example*: If a transcriptomics dataset shows upregulation of the *SPTLC1* gene, and a lipidomics dataset shows increased *Ceramide* levels, the system highlights this convergence in the Sphingolipid Metabolism pathway.

### Analysis Tools
*   **Pathway Enrichment**:
    *   **Input**: A list of significant features (e.g., differentially expressed genes).
    *   **Databases**: Queries KEGG and Reactome APIs.
    *   **Output**: A ranked table of pathways with Fisher's Exact Test p-values and FDR correction.
*   **Statistical Analysis**:
    *   **PCA**: Visualize sample clustering to detect outliers or batch effects.
    *   **Volcano Plots**: Interactive plots showing statistical significance vs. magnitude of change. Click points to see gene names.
    *   **Clustering**: Hierarchical clustering dendrograms to group similar samples or features.

### Discovery Workflow
*   **Harvesting**: Automate the scanning of public repositories (GEO, Metabolomics Workbench).
    *   *Configure*: Set keywords (e.g., "Parkinson's") and frequency (e.g., "Weekly").
*   **Review**: The Harvester creates "Draft" entries. You can review the metadata and decide to **Import** (download and process) or **Discard** them.

---

## 🧠 5. RAG & Knowledge Retrieval

### RAG Query (Chat)
Ask natural language questions about your internal data and scientific literature.
*   **Models**: Choose between **GPT-4o** (best for complex reasoning) or **Claude 3.5 Sonnet** (great for writing and coding).
*   **Modes**:
    *   **Chat**: Standard fast conversation.
    *   **Deep Reasoning**: Activates a multi-step "Chain of Thought" process where the AI breaks down the problem, plans a search strategy, and synthesizes the answer.
*   **Citations**: Answers include inline citations (e.g., `[1]`) that link directly to the source evidence (an Experiment page, a Dataset row, or a PDF snippet).

### Scientific Q&A Tracker
*   **Ask & Save**: Ask a complex scientific question (e.g., "What is the effect of inhibition of Target X on Lipid Y?"). If the answer is valuable, save it to the registry.
*   **Version Control**: Scientific knowledge evolves. You can "Re-run" a saved question months later. The system will generate a new answer based on the latest data and compare it to the old one.
*   **Evidence**: Every saved answer preserves the specific chunks of text and data used to generate it, ensuring auditability.

### Literature Analysis
*   **Upload**: Drag and drop PDF publications.
*   **Critique**: The AI reads the paper and generates a structured critique:
    *   **Strengths**: What was done well?
    *   **Weaknesses**: "Low sample size", "Missing controls".
    *   **Unanswered Questions**: What future research is suggested?

---

## 🤝 6. Collaboration & Operations

*   **Teams**: Create teams (e.g., "Biology", "Chemistry") and Projects.
    *   **Access Control**: Assign users to specific projects. Users only see data for projects they are members of.
*   **Lab Notebook**: A digital notebook for daily observations.
    *   **Linking**: You can `@mention` Experiments or Datasets within a notebook entry to create a permanent link.
*   **Scheduling**:
    *   **Calendar**: Book shared equipment (e.g., Mass Spec, PCR machine).
    *   **Gantt**: Visualize experiment timelines and dependencies.
*   **Cost Tracking**: Monitor project budgets. Track expenses for reagents, labor, and outsourcing.

---

## ⚙️ 7. Administration

*   **Audit Logs**: (Admin only) A complete, immutable history of user actions.
    *   *Fields*: User, Action (Create/Update/Delete), Entity, Timestamp, IP Address.
*   **System Health**: Monitor the pulse of the platform.
    *   **Database**: Connection status and pool usage.
    *   **Queues**: Status of background jobs (Ingestion, Harvesting).
    *   **API**: Latency and error rates.
*   **Settings**:
    *   **Theme**: Toggle between Light and Dark mode.
    *   **API Keys**: Generate personal API keys for programmatic access (see [CLI Guide](CLI_GUIDE.md)).

---

## 💡 Pro Tips
*   **Keyboard Shortcuts**: Press `?` anywhere to see available shortcuts. Common ones:
    *   `Cmd+K`: Global Search
    *   `G then D`: Go to Datasets
    *   `G then E`: Go to Experiments
*   **Feedback**: Found a bug or have an idea? Use the **User Feedback** widget in the sidebar to send a report directly to the development team.
*   **Export**: Most data tables can be exported to CSV/Excel via the download button in the top right of the table.

### Advanced Usage (CLI)
For developers or bioinformaticians preferring command-line tools for batch ingestion and automation, please refer to the [CLI Guide](CLI_GUIDE.md).

---

## 🛠️ Appendix: Integrated Tools & Technologies

This platform integrates several powerful open-source and commercial tools to provide a seamless experience. Here's a brief explanation of the key technologies powering the system:

### Artificial Intelligence (AI) & RAG
*   **Retrieval-Augmented Generation (RAG)**: A technique that combines the vast knowledge of Large Language Models (LLMs) with your private data. When you ask a question, the system first retrieves relevant documents from your database and then uses the LLM to generate an answer based *only* on those facts.
*   **Vector Search (Pinecone)**: A specialized database that stores the "meaning" of text as mathematical vectors. This allows the system to find relevant information even if the keywords don't match exactly (e.g., matching "heart attack" with "myocardial infarction").
*   **Large Language Models (LLMs)**: The "brain" behind the reasoning. We support **GPT-4o** (OpenAI) and **Claude 3.5 Sonnet** (Anthropic) for high-quality scientific reasoning.

### Cheminformatics
*   **RDKit**: An open-source cheminformatics software. We use it to:
    *   Parse chemical structures (SMILES).
    *   Calculate molecular properties (Molecular Weight, LogP, TPSA).
    *   Perform substructure searches (finding a specific core in a library).
    *   Generate molecular fingerprints for similarity searching.
*   **Ketcher**: A web-based chemical structure editor (sketcher) that allows you to draw molecules directly in the browser.

### Bioinformatics & Analysis
*   **KEGG & Reactome**: Public databases of biological pathways. We query their APIs to perform enrichment analysis (identifying which pathways are overrepresented in your gene/metabolite lists).
*   **SciPy & NumPy**: Powerful Python libraries used for statistical calculations (PCA, t-tests, Clustering).
*   **UniProt**: The universal protein resource. We use it to map gene symbols to protein IDs and retrieve functional descriptions.

### Visualization
*   **Plotly**: An interactive charting library used for most graphs (Volcano plots, Heatmaps, Bar charts). You can zoom, pan, and hover over data points for details.
*   **Cytoscape.js**: A graph theory library used for visualizing biological networks (Protein-Protein Interactions) and knowledge graphs.
*   **Ag-Grid**: An advanced data grid component that provides Excel-like filtering, sorting, and pivoting capabilities for large datasets.

---

*Last Updated: December 2025*
