"""
SQLAlchemy models for the multi-omics platform.

These models map domain models to Postgres tables. They maintain compatibility
with Notion IDs during migration by storing notion_page_id for dual-write support.
"""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from typing import List, Optional

from sqlalchemy import (
    JSON,
    BigInteger,
    Boolean,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Table,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import ARRAY, UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base
from amprenta_rag.models.chemistry import (
    compound_program,
    hts_campaign_program,
    biochemical_result_program,
)


# Helper function to generate UUID defaults
def generate_uuid() -> uuid.UUID:
    """
    Generate a new UUID for use as default.

    Returns:
        A new UUID4 instance
    """
    return uuid.uuid4()


# Association tables for many-to-many relationships

program_experiment_assoc = Table(
    "program_experiment",
    Base.metadata,
    Column("program_id", UUID(as_uuid=True), ForeignKey("programs.id"), primary_key=True),
    Column("experiment_id", UUID(as_uuid=True), ForeignKey("experiments.id"), primary_key=True),
)

program_dataset_assoc = Table(
    "program_dataset",
    Base.metadata,
    Column("program_id", UUID(as_uuid=True), ForeignKey("programs.id"), primary_key=True),
    Column("dataset_id", UUID(as_uuid=True), ForeignKey("datasets.id"), primary_key=True),
)

experiment_dataset_assoc = Table(
    "experiment_dataset",
    Base.metadata,
    Column("experiment_id", UUID(as_uuid=True), ForeignKey("experiments.id"), primary_key=True),
    Column("dataset_id", UUID(as_uuid=True), ForeignKey("datasets.id"), primary_key=True),
)

dataset_feature_assoc = Table(
    "dataset_feature",
    Base.metadata,
    Column("dataset_id", UUID(as_uuid=True), ForeignKey("datasets.id"), primary_key=True),
    Column("feature_id", UUID(as_uuid=True), ForeignKey("features.id"), primary_key=True),
)

dataset_signature_assoc = Table(
    "dataset_signature",
    Base.metadata,
    Column("dataset_id", UUID(as_uuid=True), ForeignKey("datasets.id"), primary_key=True),
    Column("signature_id", UUID(as_uuid=True), ForeignKey("signatures.id"), primary_key=True),
    Column("match_score", Float, nullable=True),
)

signature_feature_assoc = Table(
    "signature_feature",
    Base.metadata,
    Column("signature_id", UUID(as_uuid=True), ForeignKey("signatures.id"), primary_key=True),
    Column("feature_id", UUID(as_uuid=True), ForeignKey("features.id"), primary_key=True),
)

program_signature_assoc = Table(
    "program_signature",
    Base.metadata,
    Column("program_id", UUID(as_uuid=True), ForeignKey("programs.id"), primary_key=True),
    Column("signature_id", UUID(as_uuid=True), ForeignKey("signatures.id"), primary_key=True),
)


class Program(Base):
    """Research program/therapeutic area."""

    __tablename__ = "programs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(500), nullable=False, index=True)
    description = Column(Text, nullable=True)
    disease = Column(ARRAY(String), nullable=True)  # List of disease terms

    # Timestamps
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(
        DateTime,
        default=lambda: datetime.now(timezone.utc),
        onupdate=lambda: datetime.now(timezone.utc),
        nullable=False,
    )

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # External identifiers (stored as JSON for flexibility)
    external_ids = Column(JSON, nullable=True)

    # Audit
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_by: Mapped[Optional["User"]] = relationship(foreign_keys=[created_by_id])

    # Relationships
    experiments: Mapped[List["Experiment"]] = relationship(
        secondary=program_experiment_assoc, back_populates="programs"
    )
    datasets: Mapped[List["Dataset"]] = relationship(
        secondary=program_dataset_assoc, back_populates="programs"
    )
    signatures: Mapped[List["Signature"]] = relationship(
        secondary=program_signature_assoc, back_populates="programs"
    )
    compounds: Mapped[List["Compound"]] = relationship(
        "Compound", secondary=compound_program, back_populates="programs", viewonly=False
    )
    hts_campaigns: Mapped[List["HTSCampaign"]] = relationship(
        "HTSCampaign", secondary=hts_campaign_program, back_populates="programs", viewonly=False
    )
    biochemical_results: Mapped[List["BiochemicalResult"]] = relationship(
        "BiochemicalResult",
        secondary=biochemical_result_program,
        back_populates="programs",
        viewonly=False,
    )


class Experiment(Base):
    """Experimental study/assay."""

    __tablename__ = "experiments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(500), nullable=False, index=True)
    type = Column(String(100), nullable=True)  # in_vivo, in_vitro, patient, etc.
    description = Column(Text, nullable=True)
    disease = Column(ARRAY(String), nullable=True)
    matrix = Column(ARRAY(String), nullable=True)  # CSF, plasma, tissue, etc.
    model_systems = Column(ARRAY(String), nullable=True)

    # Additional experiment metadata
    targets = Column(ARRAY(String), nullable=True)  # Target molecules/proteins
    modality = Column(ARRAY(String), nullable=True)  # Treatment modality
    stage = Column(String(100), nullable=True)  # Disease stage
    biomarker_role = Column(ARRAY(String), nullable=True)  # Biomarker roles
    treatment_arms = Column(ARRAY(String), nullable=True)  # Treatment arms

    # Experimental design metadata
    design_type = Column(String(50), nullable=True)  # case_control, time_course, intervention, dose_response, multi_factorial, observational
    design_metadata = Column(JSON, nullable=True)  # Structured design details (factors, levels, etc.)
    sample_groups = Column(JSON, nullable=True)  # {"control": ["S1", "S2"], "case": ["S3", "S4"], "timepoints": {"T0": [...], "T1": [...]}}

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # Audit
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_by: Mapped[Optional["User"]] = relationship(foreign_keys=[created_by_id])

    # Project
    project_id = Column(UUID(as_uuid=True), ForeignKey("projects.id"), nullable=True)

    # Version for concurrent editing safety
    version = Column(Integer, default=1, nullable=False)

    # Archive fields for data retention
    is_archived = Column(Boolean, default=False, nullable=False)
    archived_at = Column(DateTime, nullable=True)

    # Relationships
    programs: Mapped[List["Program"]] = relationship(
        secondary=program_experiment_assoc, back_populates="experiments"
    )
    datasets: Mapped[List["Dataset"]] = relationship(
        secondary=experiment_dataset_assoc, back_populates="experiments"
    )
    project: Mapped[Optional["Project"]] = relationship(back_populates="experiments")
    generic_assay_results: Mapped[List["GenericAssayResult"]] = relationship(
        back_populates="experiment"
    )


class Feature(Base):
    """Biological feature (gene, protein, metabolite, lipid)."""

    __tablename__ = "features"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(500), nullable=False, index=True)
    feature_type = Column(String(50), nullable=False, index=True)  # gene, protein, metabolite, lipid
    normalized_name = Column(String(500), nullable=True, index=True)
    aliases = Column(ARRAY(String), nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # External identifiers
    external_ids = Column(JSON, nullable=True)

    # Relationships
    datasets: Mapped[List["Dataset"]] = relationship(
        secondary=dataset_feature_assoc, back_populates="features"
    )
    signatures: Mapped[List["Signature"]] = relationship(
        secondary=signature_feature_assoc, back_populates="features"
    )
    structures: Mapped[List["ProteinStructure"]] = relationship(
        "ProteinStructure", back_populates="feature"
    )


class Dataset(Base):
    """Experimental dataset (omics data)."""

    __tablename__ = "datasets"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(500), nullable=False, index=True)
    omics_type = Column(String(50), nullable=False, index=True)  # transcriptomics, proteomics, etc.
    description = Column(Text, nullable=True)
    file_paths = Column(ARRAY(String), nullable=True)
    file_urls = Column(ARRAY(String), nullable=True)
    organism = Column(ARRAY(String), nullable=True)
    sample_type = Column(ARRAY(String), nullable=True)
    disease = Column(ARRAY(String), nullable=True)

    # DVC versioning (optional)
    dvc_version = Column(String(32), nullable=True)
    dvc_metadata = Column(JSON, nullable=True)
    dvc_pushed = Column(Boolean, default=False, nullable=False)

    # Scientific metadata from mwTab
    methods = Column(Text, nullable=True)  # Extracted methods from mwTab
    summary = Column(Text, nullable=True)  # Study summary
    results = Column(Text, nullable=True)  # Results description
    conclusions = Column(Text, nullable=True)  # Conclusions
    dataset_source_type = Column(String(100), nullable=True)  # e.g., "Processed table"
    data_origin = Column(String(100), nullable=True)  # e.g., "External – Open Dataset"
    mwtab_json = Column(JSON, nullable=True)

    # Sample-level design metadata
    sample_group = Column(String(100), nullable=True)  # control, case, treated, untreated, placebo
    timepoint = Column(String(50), nullable=True)  # T0, T1, 24h, Week4, baseline
    intervention = Column(String(200), nullable=True)  # Drug name, treatment condition
    dose = Column(String(50), nullable=True)  # 10mg, 100nM, high, low
    replicate_id = Column(String(50), nullable=True)  # For biological/technical replicates

    # Signature matching
    signature_match_score = Column(Float, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # External identifiers (repository IDs, etc.)
    external_ids = Column(JSON, nullable=True)

    # Version for concurrent editing safety
    version = Column(Integer, default=1, nullable=False)

    # Archive fields for data retention
    is_archived = Column(Boolean, default=False, nullable=False)
    archived_at = Column(DateTime, nullable=True)

    # Ingestion status
    ingestion_status = Column(
        String(32), nullable=False, default="pending"
    )  # values: pending, in_progress, complete, failed

    # Quality metrics
    quality_score = Column(Float, nullable=True)
    quality_status = Column(String(32), nullable=True)
    quality_issues = Column(JSON, nullable=True)

    # Audit
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_by: Mapped[Optional["User"]] = relationship(foreign_keys=[created_by_id])

    # Relationships
    programs: Mapped[List["Program"]] = relationship(
        secondary=program_dataset_assoc, back_populates="datasets"
    )
    experiments: Mapped[List["Experiment"]] = relationship(
        secondary=experiment_dataset_assoc, back_populates="datasets"
    )
    features: Mapped[List["Feature"]] = relationship(
        secondary=dataset_feature_assoc, back_populates="datasets"
    )
    signatures: Mapped[List["Signature"]] = relationship(
        secondary=dataset_signature_assoc, back_populates="datasets"
    )


class SignatureComponent(Base):
    """A single component of a multi-omics signature."""

    __tablename__ = "signature_components"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    signature_id = Column(UUID(as_uuid=True), ForeignKey("signatures.id"), nullable=False, index=True)
    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=True, index=True)
    feature_name = Column(String(500), nullable=True)  # For features not yet in DB
    feature_type = Column(String(50), nullable=False)
    direction = Column(String(20), nullable=True)  # ↑, ↓, neutral, complex
    weight = Column(Float, nullable=True, default=1.0)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    signature: Mapped["Signature"] = relationship(back_populates="components")
    feature: Mapped[Optional["Feature"]] = relationship("Feature")


class Signature(Base):
    """A multi-omics signature definition."""

    __tablename__ = "signatures"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(500), nullable=False, unique=True, index=True)
    description = Column(Text, nullable=True)
    modalities = Column(ARRAY(String), nullable=True)  # List of feature types present

    # Signature metadata fields
    short_id = Column(String(100), nullable=True, index=True)  # Short identifier for the signature
    biomarker_role = Column(ARRAY(String), nullable=True)  # Biomarker roles
    phenotype_axes = Column(ARRAY(String), nullable=True)  # Phenotype axes
    data_ownership = Column(String(100), nullable=True)  # Data ownership classification

    # Validation workflow
    validation_status = Column(String(20), nullable=True, index=True)  # pending, approved, rejected

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # Audit
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_by: Mapped[Optional["User"]] = relationship(foreign_keys=[created_by_id])

    # Relationships
    components: Mapped[List["SignatureComponent"]] = relationship(
        back_populates="signature", cascade="all, delete-orphan"
    )
    features: Mapped[List["Feature"]] = relationship(
        secondary=signature_feature_assoc, back_populates="signatures"
    )
    datasets: Mapped[List["Dataset"]] = relationship(
        secondary=dataset_signature_assoc, back_populates="signatures"
    )
    programs: Mapped[List["Program"]] = relationship(
        secondary=program_signature_assoc, back_populates="signatures"
    )


class ReportArtifact(Base):
    """Generated report artifact metadata."""

    __tablename__ = "report_artifacts"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    entity_type = Column(String(50), nullable=False, index=True)  # program/dataset/signature/feature
    entity_id = Column(UUID(as_uuid=True), nullable=False, index=True)
    format = Column(String(10), nullable=False)  # html/pdf
    file_path = Column(String, nullable=False)
    params_hash = Column(String(128), nullable=False, index=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Audit
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_by: Mapped[Optional["User"]] = relationship(foreign_keys=[created_by_id])


class ReportSchedule(Base):
    """Scheduled generation of report artifacts."""

    __tablename__ = "report_schedules"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    entity_type = Column(String(50), nullable=False, index=True)  # program/dataset/signature
    entity_id = Column(UUID(as_uuid=True), nullable=True, index=True)  # allow null for "all of type"
    format = Column(String(10), nullable=False)  # html/pdf
    cron_expression = Column(String(100), nullable=False)
    enabled = Column(Boolean, default=True, nullable=False)
    last_run_at = Column(DateTime, nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    created_by: Mapped[Optional["User"]] = relationship(foreign_keys=[created_by_id])


class RepositorySubscription(Base):
    """Stored subscriptions for external repository searches."""

    __tablename__ = "repository_subscriptions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    name = Column(String(255), nullable=False)
    repository_source = Column(String(50), nullable=False, index=True)  # GEO, PRIDE, MetaboLights, MW, all
    query_params = Column(JSON, nullable=True)
    notify_email = Column(Boolean, default=False, nullable=False)
    notify_in_app = Column(Boolean, default=True, nullable=False)
    is_active = Column(Boolean, default=True, nullable=False)
    last_checked = Column(DateTime, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    user: Mapped[Optional["User"]] = relationship(foreign_keys=[user_id])


class Alert(Base):
    """Alert generated from repository subscription matches."""

    __tablename__ = "alerts"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    subscription_id = Column(UUID(as_uuid=True), ForeignKey("repository_subscriptions.id"), nullable=False)
    dataset_id = Column(UUID(as_uuid=True), ForeignKey("datasets.id"), nullable=False)
    is_read = Column(Boolean, default=False, nullable=False)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)

    subscription: Mapped["RepositorySubscription"] = relationship("RepositorySubscription")
    dataset: Mapped["Dataset"] = relationship("Dataset")


class GenomicsIndex(Base):
    """Reference transcriptome index for Salmon/Kallisto."""

    __tablename__ = "genomics_indices"

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    organism = Column(String(100), nullable=False)  # "Homo sapiens"
    tool = Column(String(20), nullable=False)  # "salmon" | "kallisto"
    version = Column(String(50), nullable=False)  # "GRCh38_v110"
    file_path = Column(String(500), nullable=False)
    file_size_bytes = Column(BigInteger)
    uploaded_by = Column(UUID(as_uuid=True), ForeignKey("users.id"))
    uploaded_at = Column(DateTime, default=datetime.utcnow)
    metadata_ = Column("metadata", JSON)

    __table_args__ = (UniqueConstraint("organism", "tool", "version", name="uix_index_org_tool_ver"),)


class PipelineJob(Base):
    """Track genomics pipeline job execution."""

    __tablename__ = "pipeline_jobs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    status = Column(String(20), default="pending")  # pending|running|complete|failed
    tool = Column(String(20), nullable=False)  # "salmon" | "kallisto"
    input_fastq_path = Column(String(500))
    index_id = Column(UUID(as_uuid=True), ForeignKey("genomics_indices.id"))
    output_dir = Column(String(500))
    result_file = Column(String(500))
    progress_percent = Column(Integer, default=0)
    error_message = Column(Text)
    started_at = Column(DateTime)
    completed_at = Column(DateTime)
    created_by = Column(UUID(as_uuid=True), ForeignKey("users.id"))
    created_at = Column(DateTime, default=datetime.utcnow)

    index = relationship("GenomicsIndex")


class NextflowJob(Base):
    """Track Nextflow pipeline executions."""

    __tablename__ = "nextflow_jobs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    pipeline_name = Column(String(100), nullable=False)  # "nf-core/rnaseq"
    pipeline_version = Column(String(20))  # "3.12.0"
    status = Column(String(20), default="pending")  # pending|running|complete|failed
    sample_sheet_path = Column(String(500))
    genome = Column(String(50))  # "GRCh38"
    output_dir = Column(String(500))
    work_dir = Column(String(500))  # Nextflow work directory
    nextflow_log = Column(String(500))  # Path to .nextflow.log
    multiqc_report = Column(String(500))  # Path to multiqc_report.html
    progress_percent = Column(Integer, default=0)
    error_message = Column(Text)
    started_at = Column(DateTime)
    completed_at = Column(DateTime)
    created_by = Column(UUID(as_uuid=True), ForeignKey("users.id"))
    created_at = Column(DateTime, default=datetime.utcnow)
    params = Column(JSON)  # Additional Nextflow params


class MLModel(Base):
    """Machine learning model registry entry."""

    __tablename__ = "ml_models"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    version = Column(String(50), nullable=False)
    model_type = Column(String(100), nullable=False)  # "admet_classification", "admet_regression"
    framework = Column(String(50), nullable=False)  # "xgboost", "sklearn"
    artifact_path = Column(String(500), nullable=False)

    # Training metadata
    training_dataset_id = Column(UUID(as_uuid=True), ForeignKey("datasets.id"), nullable=True)
    features = Column(JSON, nullable=True)  # List of feature names
    hyperparameters = Column(JSON, nullable=True)

    # Performance
    metrics = Column(JSON, nullable=True)  # {"auc": 0.92, "rmse": 0.45}

    # Lifecycle
    status = Column(String(50), default="active")  # "active", "archived", "training"
    description = Column(Text, nullable=True)

    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc))
    updated_at = Column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        onupdate=lambda: datetime.now(timezone.utc),
    )

    __table_args__ = (UniqueConstraint("name", "version", name="uq_ml_model_name_version"),)


class BatchCorrection(Base):
    """Batch effect correction record (e.g., ComBat)."""

    __tablename__ = "batch_corrections"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    method = Column(String(50), nullable=False)
    batch_map = Column(JSON, nullable=True)
    corrected_dataset_id = Column(UUID(as_uuid=True), ForeignKey("datasets.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc))


class PhenotypeGeneAssociation(Base):
    """Association between an HPO phenotype term and a gene symbol (HPOA genes_to_phenotype)."""

    __tablename__ = "phenotype_gene_associations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)

    # HPO term
    hpo_id = Column(String(20), nullable=False, index=True)  # e.g. "HP:0001250"
    hpo_name = Column(String(500), nullable=True)

    # Gene
    gene_symbol = Column(String(50), nullable=False, index=True)
    ncbi_gene_id = Column(String(30), nullable=True)

    # Optional context from HPOA
    disease_id = Column(String(50), nullable=True)
    frequency = Column(String(50), nullable=True)
    source = Column(String(100), nullable=True)

    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    __table_args__ = (
        UniqueConstraint("hpo_id", "gene_symbol", "disease_id", name="uq_hpo_gene_disease"),
    )


class GraphEdge(Base):
    """Generic evidence graph edge between two entities.

    This is a flexible table intended to represent relationships across the system
    (compound→target, target→pathway, dataset→signature, etc.).

    TODO: Consider adding a first-class `GraphNode` table once we need node-level
    metadata, canonicalization, and typed constraints beyond `provenance`.
    """

    __tablename__ = "graph_edges"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)

    source_entity_type = Column(String(50), nullable=False, index=True)
    source_entity_id = Column(UUID(as_uuid=True), nullable=False, index=True)

    target_entity_type = Column(String(50), nullable=False, index=True)
    target_entity_id = Column(UUID(as_uuid=True), nullable=False, index=True)

    relationship_type = Column(String(100), nullable=False, index=True)
    confidence = Column(Float, nullable=True)
    evidence_source = Column(String(100), nullable=True, index=True)
    provenance = Column(JSON, nullable=True)

    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        onupdate=lambda: datetime.now(timezone.utc),
        nullable=False,
    )

    __table_args__ = (
        UniqueConstraint(
            "source_entity_type",
            "source_entity_id",
            "target_entity_type",
            "target_entity_id",
            "relationship_type",
            "evidence_source",
            name="uq_graph_edge",
        ),
        Index(
            "ix_graph_edges_src_rel",
            "source_entity_type",
            "source_entity_id",
            "relationship_type",
        ),
        Index(
            "ix_graph_edges_tgt_rel",
            "target_entity_type",
            "target_entity_id",
            "relationship_type",
        ),
    )


class ProteinStructure(Base):
    """Protein structure record linked to a biological Feature (typically a gene/protein).

    Stores upstream identifiers (PDB, AlphaFold) and prep metadata; actual files are
    stored in `StructureFile`.
    """

    __tablename__ = "protein_structures"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)

    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=True, index=True)
    feature: Mapped[Optional["Feature"]] = relationship(
        "Feature", foreign_keys=[feature_id], back_populates="structures"
    )

    pdb_id = Column(String(10), nullable=True, index=True)
    alphafold_uniprot_id = Column(String(20), nullable=True, index=True)
    source = Column(String(50), nullable=False, default="unknown")  # pdb | alphafold | uploaded

    resolution = Column(Float, nullable=True)
    method = Column(String(100), nullable=True)
    chain_ids = Column(ARRAY(String), nullable=True)

    prep_status = Column(String(50), nullable=True, default="raw")
    prep_log = Column(Text, nullable=True)

    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        onupdate=lambda: datetime.now(timezone.utc),
        nullable=False,
    )

    files: Mapped[List["StructureFile"]] = relationship(
        "StructureFile", back_populates="structure", cascade="all, delete-orphan"
    )

    binding_sites: Mapped[List["BindingSite"]] = relationship(
        "BindingSite", back_populates="structure", cascade="all, delete-orphan"
    )


class StructureFile(Base):
    """File artifact associated with a ProteinStructure (PDB/mmCIF/etc.)."""

    __tablename__ = "structure_files"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    structure_id = Column(UUID(as_uuid=True), ForeignKey("protein_structures.id"), nullable=False, index=True)
    structure: Mapped["ProteinStructure"] = relationship("ProteinStructure", back_populates="files")

    file_type = Column(String(50), nullable=False)  # pdb, cif, fasta, etc.
    file_path = Column(String(500), nullable=False)
    file_size_bytes = Column(BigInteger, nullable=True)
    md5_hash = Column(String(32), nullable=True)

    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)


class BindingSite(Base):
    """Binding site / pocket detected on a ProteinStructure."""

    __tablename__ = "binding_sites"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    structure_id = Column(UUID(as_uuid=True), ForeignKey("protein_structures.id"), nullable=False, index=True)

    pocket_rank = Column(Integer, nullable=False)
    score = Column(Float, nullable=True)
    volume = Column(Float, nullable=True)

    center_x = Column(Float, nullable=True)
    center_y = Column(Float, nullable=True)
    center_z = Column(Float, nullable=True)

    residues = Column(JSON, nullable=True)  # list of "chain:res:num"
    pocket_pdb_path = Column(String(500), nullable=True)

    detection_method = Column(String(50), nullable=False, default="fpocket")
    detected_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    structure: Mapped["ProteinStructure"] = relationship(
        "ProteinStructure", foreign_keys=[structure_id], back_populates="binding_sites"
    )

    __table_args__ = (
        UniqueConstraint("structure_id", "detection_method", "pocket_rank", name="uq_binding_site_rank"),
    )


class DockingRun(Base):
    """A docking job run against a protein structure and optional binding site."""

    __tablename__ = "docking_runs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    structure_id = Column(UUID(as_uuid=True), ForeignKey("protein_structures.id"), nullable=False, index=True)
    binding_site_id = Column(UUID(as_uuid=True), ForeignKey("binding_sites.id"), nullable=True, index=True)

    center_x = Column(Float, nullable=True)
    center_y = Column(Float, nullable=True)
    center_z = Column(Float, nullable=True)

    size_x = Column(Float, nullable=True)
    size_y = Column(Float, nullable=True)
    size_z = Column(Float, nullable=True)

    status = Column(String(50), nullable=False, default="pending")  # pending|running|completed|failed
    compound_ids = Column(JSON, nullable=True)  # list of Compound.id UUID strings
    exhaustiveness = Column(Integer, nullable=False, default=8)
    total_compounds = Column(Integer, nullable=True)
    completed_compounds = Column(Integer, nullable=False, default=0)

    started_at = Column(DateTime(timezone=True), nullable=True)
    completed_at = Column(DateTime(timezone=True), nullable=True)
    error_log = Column(Text, nullable=True)

    structure: Mapped["ProteinStructure"] = relationship("ProteinStructure", foreign_keys=[structure_id])
    binding_site: Mapped[Optional["BindingSite"]] = relationship("BindingSite", foreign_keys=[binding_site_id])
    poses: Mapped[List["DockingPose"]] = relationship(
        "DockingPose", back_populates="docking_run", cascade="all, delete-orphan"
    )


class DockingPose(Base):
    """Docked pose for a compound in a DockingRun."""

    __tablename__ = "docking_poses"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    docking_run_id = Column(UUID(as_uuid=True), ForeignKey("docking_runs.id"), nullable=False, index=True)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=False, index=True)

    binding_affinity = Column(Float, nullable=True)  # kcal/mol
    pose_rank = Column(Integer, nullable=True)
    rmsd_lb = Column(Float, nullable=True)
    rmsd_ub = Column(Float, nullable=True)

    pose_pdbqt_path = Column(String(500), nullable=True)
    docked_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    docking_run: Mapped["DockingRun"] = relationship("DockingRun", back_populates="poses", foreign_keys=[docking_run_id])
    compound: Mapped["Compound"] = relationship("Compound", foreign_keys=[compound_id])
    quality: Mapped[Optional["PoseQuality"]] = relationship(
        "PoseQuality", back_populates="pose", cascade="all, delete-orphan", uselist=False
    )
    interactions: Mapped[List["PoseInteraction"]] = relationship(
        "PoseInteraction", back_populates="pose", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_docking_poses_run_rank", "docking_run_id", "pose_rank"),
    )


class PoseQuality(Base):
    """Per-pose quality/control metrics and interaction summary (e.g., PLIP)."""

    __tablename__ = "pose_qualities"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    pose_id = Column(UUID(as_uuid=True), ForeignKey("docking_poses.id"), nullable=False, index=True, unique=True)

    num_hbonds = Column(Integer, nullable=False, default=0)
    num_hydrophobic = Column(Integer, nullable=False, default=0)
    num_salt_bridges = Column(Integer, nullable=False, default=0)
    num_pi_stacking = Column(Integer, nullable=False, default=0)
    num_pi_cation = Column(Integer, nullable=False, default=0)
    num_halogen = Column(Integer, nullable=False, default=0)
    num_metal = Column(Integer, nullable=False, default=0)
    total_interactions = Column(Integer, nullable=False, default=0)

    has_clashes = Column(Boolean, nullable=False, default=False)
    ligand_efficiency = Column(Float, nullable=True)

    analyzed_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    pose: Mapped["DockingPose"] = relationship("DockingPose", foreign_keys=[pose_id], back_populates="quality")


class PoseInteraction(Base):
    """Detailed interactions for a pose (e.g., PLIP-derived)."""

    __tablename__ = "pose_interactions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    pose_id = Column(UUID(as_uuid=True), ForeignKey("docking_poses.id"), nullable=False, index=True)

    interaction_type = Column(String(50), nullable=False)  # hbond|hydrophobic|saltbridge|pistacking|pication|halogen|metal
    ligand_atom = Column(String(200), nullable=True)
    protein_residue = Column(String(200), nullable=True)
    distance = Column(Float, nullable=True)
    angle = Column(Float, nullable=True)

    pose: Mapped["DockingPose"] = relationship("DockingPose", foreign_keys=[pose_id], back_populates="interactions")

    __table_args__ = (Index("ix_pose_interactions_pose_type", "pose_id", "interaction_type"),)


gene_protein_map = Table(
    "gene_protein_map",
    Base.metadata,
    Column("gene_symbol", String, index=True, nullable=False),
    Column("uniprot_id", String, index=True, nullable=False),
    Column("meta", JSON, nullable=True),
)

feature_pathway_map = Table(
    "feature_pathway_map",
    Base.metadata,
    Column("feature_id", UUID(as_uuid=True), index=True, nullable=False),
    Column("pathway_id", UUID(as_uuid=True), index=True, nullable=False),
    Column("meta", JSON, nullable=True),
)

# Re-export models from subpackages for backward compatibility
from amprenta_rag.models.chemistry import (  # noqa: F401
    Compound,
    HTSCampaign,
    HTSResult,
    BiochemicalResult,
    ActivityResult,
    TargetProductProfile,
    CandidateNomination,
    GenericAssayResult,
    ADMEResult,  # noqa: F401
    PKStudy,  # noqa: F401
    ToxicologyResult,  # noqa: F401
)
from amprenta_rag.models.auth import (  # noqa: F401
    User,
    Team,
    TeamMember,
    Project,
    FeaturePermission,
    AuditLog,
)
from amprenta_rag.models.automation import (  # noqa: F401
    WorkflowRule,
    WorkflowExecution,
)
from amprenta_rag.models.content import (  # noqa: F401
    Literature,
    RAGChunk,
    LiteratureCritique,
)
from amprenta_rag.models.discovery import (  # noqa: F401
    DiscoveredStudy,
    HarvestSchedule,
    DiscoveryJob,
)
from amprenta_rag.models.user_prefs import (  # noqa: F401
    Note,
    Bookmark,
    Notification,
    UserFavorite,
    EmailSubscription,
    SavedFilter,
    Comment,
)
from amprenta_rag.models.content import Email  # noqa: F401
from amprenta_rag.models.sample import (  # noqa: F401
    Sample,
    StorageLocation,
    SampleTransfer,
)
from amprenta_rag.models.qa import (  # noqa: F401
    SavedQuestion,
    SavedAnswer,
)
from amprenta_rag.models.misc import (  # noqa: F401
    ScheduledEvent,
    RetentionPolicy,
    OntologyTerm,
    ExperimentTemplate,  # noqa: F401
    GeneticVariant,  # noqa: F401
    CostEntry,  # noqa: F401
)
from amprenta_rag.models.eln import (  # noqa: F401
    Protocol,
    ExperimentProtocol,
    LabNotebookEntry,
    LabNotebookEntryAssociation,
)
