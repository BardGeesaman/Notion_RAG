"""Dashboard configuration and registry."""

from __future__ import annotations

import os

AUTH_DISABLED = os.environ.get("DISABLE_AUTH", "").lower() in ("1", "true", "yes")

# Legacy page lists - DEPRECATED, use PAGE_GROUPS instead
# Kept for backward compatibility during migration
# TODO: Remove after all imports are migrated to PAGE_GROUPS

# This will be defined after PAGE_REGISTRY and PAGE_GROUPS are defined below

# Mapping page name -> (module_path, function_name)
PAGE_REGISTRY = {
    "Overview": ("scripts.dashboard.pages.overview", "render_overview_page"),
    "Cockpit": ("scripts.dashboard.pages.cockpit", "render_cockpit_page"),
    "Workspaces": ("scripts.dashboard.pages.Workspaces", "render_workspaces_page"),
    "Getting Started": ("scripts.dashboard.pages.getting_started", "render_getting_started_page"),
    "Evaluation Wizard": ("scripts.dashboard.pages.evaluation_wizard", "render_evaluation_wizard"),
    "Chat": ("scripts.dashboard.pages.chat", "render_chat_page"),
    "Lab Notebook": ("scripts.dashboard.pages.lab_notebook", "render_lab_notebook_page"),
    "Sample Inventory": ("scripts.dashboard.pages.sample_inventory", "render_sample_inventory_page"),
    "Search": ("scripts.dashboard.pages.search", "render_search_page"),
    "Data Ingestion": ("scripts.dashboard.pages.ingestion", "render_ingestion_page"),
    "Repositories": ("scripts.dashboard.pages.repositories", "render_repositories_page"),
    "Discovery Workflow": ("scripts.dashboard.pages.discovery_workflow", "render_discovery_workflow_page"),
    "Variant Tracking": ("scripts.dashboard.pages.variants", "render_variants_page"),
    "Activity Feed": ("scripts.dashboard.pages.activity_feed", "render_activity_feed_page"),
    "Share Links": ("scripts.dashboard.pages.share_links", "render_share_links_page"),
    "Analysis Tools": ("scripts.dashboard.pages.analysis", "render_analysis_page"),
    "Experiment Planner": ("scripts.dashboard.pages.experiment_planner", "render_experiment_planner_page"),
    "Visualizations": ("scripts.dashboard.pages.visualizations", "render_visualizations_page"),
    "Cytoscape Demo": ("scripts.dashboard.pages.visualizations.cytoscape_network", "render"),
    "Entity Explorer": ("scripts.dashboard.pages.visualizations.entity_explorer", "render"),
    "Data Grid": ("scripts.dashboard.pages.visualizations.data_grid", "render"),
    "Genome Browser": ("scripts.dashboard.pages.visualizations.genome_browser", "render"),
    "Quality Checks": ("scripts.dashboard.pages.quality_checks", "render_quality_checks_page"),
    "Statistical Analysis": ("scripts.dashboard.pages.statistical_analysis", "render_statistical_analysis_page"),
    "Discovery": ("scripts.dashboard.pages.discovery", "render_discovery_page"),
    "Coverage Map": ("scripts.dashboard.pages.coverage", "render_coverage_page"),
    "Feature Recurrence": ("scripts.dashboard.pages.feature_recurrence", "render_feature_recurrence_page"),
    "Evidence Report": ("scripts.dashboard.pages.evidence_report", "render_evidence_report_page"),
    "Data Management": ("scripts.dashboard.pages.management", "render_management_page"),
    "System Health": ("scripts.dashboard.pages.system_health", "render_system_health_page"),
    "Audit Logs": ("scripts.dashboard.pages.audit_logs", "render_audit_logs_page"),
    "Audit Trail": ("scripts.dashboard.pages.audit_trail", "render_audit_trail_page"),
    "IP Portfolio": ("scripts.dashboard.pages.ip_portfolio", "render_ip_portfolio_page"),
    "Data Export": ("scripts.dashboard.pages.data_export", "render_data_export_page"),
    "Relationships": ("scripts.dashboard.pages.relationships", "render_relationships_page"),
    "Datasets": ("scripts.dashboard.pages.datasets", "render_datasets_page"),
    "Programs": ("scripts.dashboard.pages.programs", "render_programs_page"),
    "Experiments": ("scripts.dashboard.pages.experiments", "render_experiments_page"),
    "Protocols": ("scripts.dashboard.pages.protocols", "render_protocols_page"),
    "Q&A Tracker": ("scripts.dashboard.pages.qa_tracker", "render_qa_tracker_page"),
    "Teams & Projects": ("scripts.dashboard.pages.teams", "render_teams_page"),
    "Feedback": ("scripts.dashboard.pages.feedback", "render_feedback_page"),
    "Features": ("scripts.dashboard.pages.features", "render_features_page"),
    "Signatures": ("scripts.dashboard.pages.signatures", "render_signatures_page"),
    "Literature": ("scripts.dashboard.pages.literature", "render_literature_page"),
    "Emails": ("scripts.dashboard.pages.emails", "render_emails_page"),
    "RAG Chunks": ("scripts.dashboard.pages.rag_chunks", "render_rag_chunks_page"),
    "Chemistry": ("scripts.dashboard.pages.chemistry", "render_chemistry_page"),
    "Chemical Sketcher": ("scripts.dashboard.pages.chemical_sketcher", "render_chemical_sketcher_page"),
    "Generative Chemistry": ("scripts.dashboard.pages.generative_chemistry", "render_generative_chemistry_page"),
    "Image Analysis": ("scripts.dashboard.pages.image_analysis", "render_image_analysis_page"),
    "Compound Portfolio": ("scripts.dashboard.pages.compound_portfolio", "render_compound_portfolio_page"),
    "HTS QC": ("scripts.dashboard.pages.hts_qc", "render_hts_qc_page"),
    "MOA Inference": ("scripts.dashboard.pages.moa_inference", "render_moa_inference_page"),
    "SAR What-If": ("scripts.dashboard.pages.sar_whatif", "render_sar_whatif_page"),
    "ADMET Predictor": ("scripts.dashboard.pages.admet_predictor", "render_admet_predictor_page"),
    "Structural Alerts": ("scripts.dashboard.pages.structural_alerts", "render_structural_alerts_page"),
    "Target QSAR": ("scripts.dashboard.pages.target_qsar", "render_target_qsar_page"),
    "Biomarker Discovery": ("scripts.dashboard.pages.biomarker_discovery", "render_biomarker_discovery_page"),
    "Molecule Viewer": ("scripts.dashboard.pages.molecule_viewer", "render_molecule_viewer_page"),
    "Compound-Target Network": ("scripts.dashboard.pages.compound_target_network", "render_compound_target_network_page"),
    "Dose-Response Explorer": ("scripts.dashboard.pages.dose_response_explorer", "render_dose_response_explorer_page"),
    "Pathway Map Viewer": ("scripts.dashboard.pages.pathway_map_viewer", "render_pathway_map_viewer_page"),
    "Compound Ranking": ("scripts.dashboard.pages.compound_ranking", "render_compound_ranking_page"),
    "AI Extraction": ("scripts.dashboard.pages.ai_extraction", "render_ai_extraction_page"),
    "Sync Monitor": ("scripts.dashboard.pages.sync_monitor", "render_sync_monitor_page"),
    "Cross-Omics Pathways": ("scripts.dashboard.pages.cross_omics_pathways", "render_cross_omics_pathways_page"),
    "RAG Query": ("scripts.dashboard.pages.rag_query", "render_rag_query_page"),
    "Generic Assays": ("scripts.dashboard.pages.generic_assays", "render_generic_assays_page"),
    "Cross-Omics": ("scripts.dashboard.pages.cross_omics", "render_cross_omics_page"),
    "Import Data": ("scripts.dashboard.pages.import_data", "render_import_page"),
    "Compare": ("scripts.dashboard.pages.compare", "render_compare_page"),
    "Timeline": ("scripts.dashboard.pages.timeline", "render_timeline_page"),
    "Data Quality": ("scripts.dashboard.pages.data_quality", "render_data_quality_page"),
    "Workflows": ("scripts.dashboard.pages.workflows", "render_workflows_page"),
    "Literature Analysis": ("scripts.dashboard.pages.literature_analysis", "render_literature_analysis_page"),
    "Paper Search": ("scripts.dashboard.pages.paper_search", "render_paper_search_page"),
    "Screening": ("scripts.dashboard.pages.screening", "render_screening_page"),
    "Predictors": ("scripts.dashboard.pages.predictors", "render_predictors_page"),
    "Scoring": ("scripts.dashboard.pages.scoring", "render_scoring_page"),
    "Phenotypes": ("scripts.dashboard.pages.phenotypes", "render_phenotypes_page"),
    "Projector": ("scripts.dashboard.pages.projector", "render_projector_page"),
    "Candidate Selection": ("scripts.dashboard.pages.candidate_selection", "render_candidate_selection_page"),
    "Email Settings": ("scripts.dashboard.pages.email_settings", "render_email_settings_page"),
    "Data Lineage": ("scripts.dashboard.pages.data_lineage", "render_data_lineage_page"),
    "Feature Permissions": ("scripts.dashboard.pages.feature_permissions", "render_feature_permissions_page"),
    "Schedule": ("scripts.dashboard.pages.schedule", "render_schedule_page"),
    "Cost Tracking": ("scripts.dashboard.pages.cost_tracking", "render_cost_tracking_page"),
    "Data Retention": ("scripts.dashboard.pages.retention", "render_retention_page"),
    "Backup Admin": ("scripts.dashboard.pages.backup_admin", "render_backup_admin_page"),
    "Ontology Management": ("scripts.dashboard.pages.ontology", "render_ontology_page"),
    "Report History": ("scripts.dashboard.pages.report_history", "render_report_history_page"),
    "Notebook Co-Pilot": ("scripts.dashboard.pages.notebook_copilot", "render_notebook_copilot_page"),
    "Notebook Generator": ("scripts.dashboard.pages.notebook_generator", "render_notebook_generator_page"),
    "Notebook Gallery": ("scripts.dashboard.pages.notebook_gallery", "render_notebook_gallery_page"),
    "Notebooks": ("scripts.dashboard.pages.notebooks", "render_notebooks_page"),
    "Executive Digests": ("scripts.dashboard.pages.digest_manager", "render_digest_manager_page"),
    "Review Queue": ("scripts.dashboard.pages.review_queue", "render_review_queue_page"),
    "AutoML Templates": ("scripts.dashboard.pages.automl_launcher", "render_automl_launcher_page"),
    "Pipeline Runner": ("scripts.dashboard.pages.pipeline_runner", "render_pipeline_runner_page"),
    "Nextflow Orchestrator": ("scripts.dashboard.pages.nextflow_orchestrator", "render_nextflow_orchestrator_page"),
    "Experiment Optimizer": ("scripts.dashboard.pages.experiment_optimizer", "render_experiment_optimizer_page"),
    "Model Registry": ("scripts.dashboard.pages.model_registry", "render_model_registry_page"),
    "Model Monitoring": ("scripts.dashboard.pages.model_monitoring", "render_model_monitoring_page"),
    "Batch Correction": ("scripts.dashboard.pages.batch_correction", "render_batch_correction_page"),
    "Sphingolipid Imbalance": ("scripts.dashboard.pages.sphingolipid_imbalance", "render_sphingolipid_imbalance_page"),
    "Graph Explorer": ("scripts.dashboard.pages.graph_explorer", "render_graph_explorer_page"),
    "Protein Structures": ("scripts.dashboard.pages.protein_structures", "render_protein_structures_page"),
    "Binding Sites": ("scripts.dashboard.pages.binding_sites", "render_binding_sites_page"),
    "Docking Runs": ("scripts.dashboard.pages.docking_runs", "render_docking_runs_page"),
    "Docking Triage": ("scripts.dashboard.pages.docking_triage", "render_docking_triage_page"),
    "Connectivity Map": ("scripts.dashboard.pages.connectivity_map", "render_connectivity_map_page"),
    "Single-Cell Viewer": ("scripts.dashboard.pages.single_cell_viewer", "render_single_cell_viewer_page"),
    "Flow Cytometry": ("scripts.dashboard.pages.flow_cytometry", "render_flow_cytometry_page"),
    "Spectral Matching": ("scripts.dashboard.pages.spectral_matching", "render_spectral_matching_page"),
    "CRISPR Analysis": ("scripts.dashboard.pages.crispr_analysis", "render_crispr_analysis_page"),
    "Multi-Omics Integration": ("scripts.dashboard.pages.multi_omics_integration", "render_multi_omics_integration_page"),
    "Variant Analysis": ("scripts.dashboard.pages.variant_analysis", "render_variant_analysis_page"),
    "Company Settings": ("scripts.dashboard.pages.company_settings", "render_company_settings_page"),
    "Job Queue": ("scripts.dashboard.pages.job_queue", "render_job_queue_page"),
}

# Navigation groups for organized sidebar
PAGE_GROUPS = {
    "Home": ["Overview", "Cockpit", "Getting Started"],
    "Discovery": ["Paper Search", "Datasets", "Experiments", "Programs", "Literature", "Review Queue", "Discovery Workflow", "Discovery", "Repositories", "Search", "Features", "Signatures", "Workspaces"],
    "Chemistry": ["Chemistry", "Chemical Sketcher", "Generative Chemistry", "Compound Portfolio", "ADMET Predictor", "Structural Alerts", "Compound Ranking", "Compound-Target Network", "SAR What-If", "Candidate Selection"],
    "HTS": ["HTS QC", "Dose-Response Explorer", "Screening"],
    "Structure": ["Protein Structures", "Binding Sites", "Docking Runs", "Docking Triage", "Molecule Viewer"],
    "Omics": ["Single-Cell Viewer", "CRISPR Analysis", "Variant Analysis", "Multi-Omics Integration", "Connectivity Map", "Biomarker Discovery", "Cross-Omics", "Cross-Omics Pathways", "Sphingolipid Imbalance"],
    "ML/Predictive": ["Target QSAR", "MOA Inference", "Predictors", "Scoring", "Phenotypes", "Model Monitoring", "Model Registry", "AutoML Templates", "Experiment Optimizer"],
    "Visualization": ["Projector", "Graph Explorer", "Pathway Map Viewer", "Visualizations", "Cytoscape Demo", "Entity Explorer", "Data Grid", "Genome Browser"],
    "Notebooks": ["Notebook Co-Pilot", "Notebook Generator", "Notebook Gallery", "Notebooks"],
    "Pipelines": ["Nextflow Orchestrator", "Pipeline Runner", "Sync Monitor", "AI Extraction", "Batch Correction"],
    "Analysis": ["Analysis Tools", "Experiment Planner", "RAG Query", "Statistical Analysis", "Quality Checks", "Literature Analysis", "Compare", "Timeline", "Data Quality", "Coverage Map", "Feature Recurrence", "Evidence Report", "Image Analysis", "Flow Cytometry"],
    "Collaboration": ["Activity Feed", "Teams & Projects", "Share Links", "Variant Tracking"],
    "Reports": ["Executive Digests", "Spectral Matching", "Report History"],
    "Admin": ["Company Settings", "Job Queue", "Audit Logs", "Audit Trail", "IP Portfolio", "Data Export", "Import Data", "Data Management", "System Health", "Workflows", "Feature Permissions", "Data Retention", "Backup Admin", "Ontology Management", "Cost Tracking", "Schedule", "Email Settings", "Data Lineage", "Feedback", "Data Ingestion"],
    "Other": ["Lab Notebook", "Sample Inventory", "Chat", "Evaluation Wizard", "Emails", "RAG Chunks", "Protocols", "Q&A Tracker", "Relationships", "Generic Assays"],
}

GROUP_ORDER = ["Home", "Discovery", "Chemistry", "HTS", "Structure", "Omics", "ML/Predictive", "Visualization", "Analysis", "Notebooks", "Pipelines", "Collaboration", "Reports", "Admin", "Other"]

GROUP_ICONS = {
    "Home": "🏠", "Discovery": "🔍", "Chemistry": "⚗️", "HTS": "🧪",
    "Structure": "🔬", "Omics": "🧬", "ML/Predictive": "🤖", "Visualization": "📊",
    "Analysis": "📈", "Notebooks": "📓", "Pipelines": "🔄", "Collaboration": "👥",
    "Reports": "📋", "Admin": "⚙️", "Other": "📁"
}

# Derived page lists for backward compatibility
# These are derived from PAGE_REGISTRY and PAGE_GROUPS to maintain consistency
ALL_PAGES = list(PAGE_REGISTRY.keys())
ADMIN_PAGES = PAGE_GROUPS.get("Admin", [])


