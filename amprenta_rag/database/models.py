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
