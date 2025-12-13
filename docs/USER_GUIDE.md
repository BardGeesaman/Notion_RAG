# Amprenta Platform - User Guide

Welcome to the Amprenta Multi-Omics RAG System! This comprehensive guide will help you navigate the platform, manage your data, and leverage advanced analytics to drive scientific discovery.

---

## 🚀 1. Getting Started

### Access & Login
*   **URL**: Access the platform via your deployed URL or `localhost:8501` for local instances.
*   **Login**: Use your assigned credentials. If you don't have an account, contact your administrator (or use the "Register" tab if enabled).
*   **Permissions**: Your access to specific features (e.g., Admin tools, Chemistry) is determined by your user role.

### Dashboard Overview
The **Overview** page is your command center.
*   **Metrics**: See high-level stats (Total Experiments, Active Compounds, etc.).
*   **Activity Feed**: Recent actions taken by your team.
*   **Widgets**: Customizable cards showing recent data or alerts.
*   **Global Search**: Use the search bar (or `Cmd+K` / `Ctrl+K`) to jump to any entity (Gene, Compound, Experiment) instantly.

### Navigation
The sidebar is organized into functional groups:
*   **Data Management**: Core entities (Programs, Experiments, Datasets).
*   **Chemistry**: Compound registration and analysis.
*   **Discovery**: Workflows, harvesting, and signatures.
*   **Analysis**: Statistical tools, pathway enrichment, and visualizations.
*   **Knowledge**: RAG Query, Q&A Tracker, Literature.
*   **Collaboration**: Teams, Lab Notebook, Scheduling.
*   **Admin**: System health, logs, and settings.

---

## 🧪 2. Data Management (The Foundation)

### Programs & Experiments
*   **Programs**: Top-level research areas (e.g., "Neurodegeneration"). View all related experiments and datasets here.
*   **Experiments**: Specific studies. You can create new experiments, assign protocols, and link datasets.
*   **Protocols**: Library of standard operating procedures (SOPs). Link these to experiments for reproducibility.

### Datasets & Samples
*   **Data Ingestion**: Use the **Import Data** page to upload raw files (CSV/Excel).
    *   Supports: Transcriptomics (RNA-seq), Proteomics, Metabolomics, Lipidomics.
    *   Validation: The system automatically checks for required columns and data types.
*   **Sample Inventory**: Track physical samples, their locations (Freezer/Box), and lineage (Parent -> Child samples).
*   **Datasets**: Browse imported data. Click a dataset to view its features (genes/metabolites), metadata, and summary statistics.

---

## ⚗️ 3. Chemistry & Lead Optimization

### Compound Management
*   **Registration**: Use the **Chemistry > Compounds** page to register new molecules.
    *   **Input**: Enter SMILES string or draw using the sketcher (if enabled).
    *   **Checks**: The system automatically checks for duplicates and calculates properties (MW, LogP).
    *   **ID**: Assigns a unique corporate ID (e.g., `AMP-00123`).

### Screening & SAR
*   **Assay Results**: View biochemical and cellular assay data associated with compounds.
*   **SAR Analysis**:
    *   **R-Group Decomposition**: Analyze how structural changes affect activity.
    *   **Activity Cliffs**: Visualize pairs of similar molecules with drastically different potency.
    *   **Pharmacophore Search**: Find compounds matching specific 3D feature constraints.

### Lead Optimization
*   **Candidate Selection**: Track compounds through the funnel (Hit -> Lead -> Candidate).
*   **TPP Scoring**: Compare compounds against your Target Product Profile (TPP) with "Traffic Light" scoring (Green/Yellow/Red).
*   **Vendor Search**: Search external catalogs (MolPort, Mcule) to procure analogs.

---

## 🧬 4. Omics Analysis & Discovery

### Signatures & Patterns
*   **Signatures**: Collections of features (genes, proteins) that change together.
    *   **Discovery**: The system automatically identifies signatures from datasets.
    *   **Scoring**: Score new datasets against established signatures to find biological connections.
*   **Cross-Omics Analysis**: Find convergence where transcriptomics, proteomics, and metabolomics point to the same pathway.

### Analysis Tools
*   **Pathway Enrichment**: Map your feature lists to KEGG and Reactome pathways.
    *   **Input**: List of genes/proteins/metabolites.
    *   **Output**: Ranked list of enriched pathways with FDR-corrected p-values.
*   **Statistical Analysis**: Run PCA, Volcano plots, and clustering directly in the browser.
*   **Visualizations**: Interactive heatmaps, network graphs (Cytoscape), and dendrograms.

### Discovery Workflow
*   **Harvesting**: Automate the scanning of public repositories (GEO, Metabolomics Workbench).
*   **Review**: Review "Draft" studies found by the harvester and choose to import them or discard them.

---

## 🧠 5. RAG & Knowledge Retrieval

### RAG Query (Chat)
Ask natural language questions about your internal data and scientific literature.
*   **Models**: Choose between GPT-4o, Claude 3.5, etc.
*   **Modes**:
    *   **Chat**: Standard conversation.
    *   **Deep Reasoning**: Multi-step "Chain of Thought" for complex problems.
*   **Citations**: Answers include inline citations `[1]` linked to the source data (Experiment X, Paper Y).

### Scientific Q&A Tracker
*   **Ask & Save**: Ask a complex scientific question once, and save the answer.
*   **Version Control**: Re-run the question later to see if new data changes the answer.
*   **Evidence**: Every saved answer preserves the retrieved context chunks for verification.

### Literature Analysis
*   **Upload**: Ingest PDFs or search Semantic Scholar/OpenAlex (if configured).
*   **Critique**: The AI automatically generates critical summaries, identifying "Strengths," "Weaknesses," and "Unanswered Questions."

---

## 🤝 6. Collaboration & Operations

*   **Teams**: Create teams and projects. Assign users to specific projects to control data access.
*   **Lab Notebook**: A digital notebook for daily observations. Entries can be linked to experiments.
*   **Scheduling**:
    *   **Calendar**: Book shared equipment.
    *   **Gantt**: View experiment timelines.
*   **Cost Tracking**: Monitor project budgets and experiment costs.

---

## ⚙️ 7. Administration

*   **Audit Logs**: (Admin only) View a complete history of user actions (login, data export, edits).
*   **System Health**: Monitor database status, API latency, and background job queues.
*   **Settings**:
    *   **Theme**: Toggle Light/Dark mode.
    *   **API Keys**: Manage your personal API keys (if enabled).

---

## 💡 Pro Tips
*   **Keyboard Shortcuts**: Press `?` anywhere to see available shortcuts.
*   **Feedback**: Found a bug or have an idea? Use the **User Feedback** widget in the sidebar.
*   **Export**: Most data tables can be exported to CSV/Excel via the download button in the top right of the table.

### Advanced Usage (CLI)
For developers or bioinformaticians preferring command-line tools for batch ingestion and automation, please refer to the [CLI Guide](CLI_GUIDE.md).

---

*Last Updated: December 2025*

