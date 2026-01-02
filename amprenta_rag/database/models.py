"""
SQLAlchemy models for the multi-omics platform.

These models map domain models to Postgres tables. They maintain compatibility
with Notion IDs during migration by storing notion_page_id for dual-write support.
"""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from enum import Enum as PyEnum
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
    func,
)
from sqlalchemy.dialects.postgresql import ARRAY, JSONB, UUID
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


class ActivityEventType(str, PyEnum):
    """Types of activity events that can be tracked."""
    
    COMPOUND_ADDED = "compound_added"
    EXPERIMENT_CREATED = "experiment_created"
    MODEL_TRAINED = "model_trained"
    HIT_CONFIRMED = "hit_confirmed"
    STATUS_CHANGED = "status_changed"
    NOTEBOOK_REVIEWED = "notebook_reviewed"
    COMMENT_ADDED = "comment_added"
    MENTION_RECEIVED = "mention_received"
    ENTITY_SHARED = "entity_shared"
    ENTITY_UNSHARED = "entity_unshared"
    REVIEW_SUBMITTED = "review_submitted"
    REVIEW_ASSIGNED = "review_assigned"
    REVIEW_DECIDED = "review_decided"
    REVIEW_REMINDER = "review_reminder"
    REVIEW_ESCALATED = "review_escalated"


class LifecycleStatus(str, PyEnum):
    """Lifecycle status for data entities."""
    ACTIVE = "active"           # Normal, visible
    QUARANTINED = "quarantined" # Hidden but recoverable, under review
    INVALID = "invalid"         # Visible with warning, known issues
    ARCHIVED = "archived"       # Soft delete, hidden


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
    pinned_dashboards: Mapped[List["PinnedDashboard"]] = relationship(
        "PinnedDashboard", back_populates="program", cascade="all, delete-orphan"
    )


class PinnedDashboard(Base):
    """Notebook dashboards pinned/locked to a Program."""

    __tablename__ = "pinned_dashboards"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    notebook_path = Column(String(500), nullable=False, index=True)
    program_id = Column(UUID(as_uuid=True), ForeignKey("programs.id", ondelete="CASCADE"), nullable=False, index=True)
    pinned_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"), nullable=True, index=True)
    pinned_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    display_name = Column(String(500), nullable=True)
    config = Column(JSONB, nullable=True)

    program: Mapped["Program"] = relationship("Program", back_populates="pinned_dashboards")
    pinned_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[pinned_by_id])


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

    # Lifecycle status (replaces is_archived pattern)
    lifecycle_status = Column(String(20), default='active', nullable=False, index=True)
    retention_exempt = Column(Boolean, default=False, nullable=False)

    @property
    def is_archived_compat(self) -> bool:
        """Backward compatibility: Check if entity is hidden."""
        return self.lifecycle_status in ('archived', 'quarantined')

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

    # Lifecycle status (replaces is_archived pattern)
    lifecycle_status = Column(String(20), default='active', nullable=False, index=True)
    retention_exempt = Column(Boolean, default=False, nullable=False)

    @property
    def is_archived_compat(self) -> bool:
        """Backward compatibility: Check if entity is hidden."""
        return self.lifecycle_status in ('archived', 'quarantined')

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


class SingleCellDataset(Base):
    """Single-cell dataset metadata linked to a Dataset (AnnData / .h5ad)."""

    __tablename__ = "single_cell_datasets"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    dataset_id = Column(UUID(as_uuid=True), ForeignKey("datasets.id"), nullable=False, index=True)

    h5ad_path = Column(String(500), nullable=False)
    file_size_bytes = Column(BigInteger, nullable=True)
    n_cells = Column(Integer, nullable=True)
    n_genes = Column(Integer, nullable=True)

    normalization_method = Column(String(100), nullable=True)
    hvg_count = Column(Integer, nullable=True)
    clustering_method = Column(String(100), nullable=True)
    clustering_resolution = Column(Float, nullable=True)

    has_pca = Column(Boolean, nullable=False, default=False)
    has_umap = Column(Boolean, nullable=False, default=False)
    processing_status = Column(String(50), nullable=False, default="pending")
    processing_log = Column(Text, nullable=True)

    ingested_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    processed_at = Column(DateTime(timezone=True), nullable=True)

    dataset: Mapped["Dataset"] = relationship("Dataset", foreign_keys=[dataset_id])
    annotations: Mapped[List["CellAnnotation"]] = relationship(
        "CellAnnotation", back_populates="single_cell_dataset", cascade="all, delete-orphan"
    )
    clusters: Mapped[List["CellCluster"]] = relationship(
        "CellCluster", back_populates="single_cell_dataset", cascade="all, delete-orphan"
    )
    markers: Mapped[List["CellTypeMarker"]] = relationship(
        "CellTypeMarker", back_populates="single_cell_dataset", cascade="all, delete-orphan"
    )


class CellAnnotation(Base):
    """Per-cell annotations and embeddings for a SingleCellDataset."""

    __tablename__ = "cell_annotations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    single_cell_dataset_id = Column(
        UUID(as_uuid=True), ForeignKey("single_cell_datasets.id"), nullable=False, index=True
    )

    barcode = Column(String(200), nullable=False)
    cluster_id = Column(Integer, nullable=True)
    cell_type = Column(String(200), nullable=True)

    umap_1 = Column(Float, nullable=True)
    umap_2 = Column(Float, nullable=True)

    n_genes_detected = Column(Integer, nullable=True)
    total_counts = Column(Float, nullable=True)
    pct_mitochondrial = Column(Float, nullable=True)

    single_cell_dataset: Mapped["SingleCellDataset"] = relationship(
        "SingleCellDataset", foreign_keys=[single_cell_dataset_id], back_populates="annotations"
    )

    __table_args__ = (
        UniqueConstraint("single_cell_dataset_id", "barcode", name="uq_cell_annotation_barcode"),
        Index("ix_cell_annotations_dataset_cluster", "single_cell_dataset_id", "cluster_id"),
    )


class CellCluster(Base):
    """Cluster-level annotation for a SingleCellDataset."""

    __tablename__ = "cell_clusters"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    single_cell_dataset_id = Column(
        UUID(as_uuid=True), ForeignKey("single_cell_datasets.id"), nullable=False, index=True
    )
    cluster_id = Column(Integer, nullable=False)

    cell_type = Column(String(200), nullable=True)
    n_cells = Column(Integer, nullable=True)
    marker_feature_ids = Column(JSON, nullable=True)  # list[UUID] or list[str] depending on pipeline
    description = Column(Text, nullable=True)

    single_cell_dataset: Mapped["SingleCellDataset"] = relationship(
        "SingleCellDataset", foreign_keys=[single_cell_dataset_id], back_populates="clusters"
    )

    __table_args__ = (
        UniqueConstraint("single_cell_dataset_id", "cluster_id", name="uq_cell_cluster"),
    )


class CellTypeMarker(Base):
    """Differential expression marker for a cluster/cell type."""

    __tablename__ = "cell_type_markers"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    single_cell_dataset_id = Column(
        UUID(as_uuid=True), ForeignKey("single_cell_datasets.id"), nullable=False, index=True
    )
    cluster_id = Column(Integer, nullable=False)

    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=True, index=True)
    gene_symbol = Column(String(50), nullable=True, index=True)

    log2_fold_change = Column(Float, nullable=True)
    pval = Column(Float, nullable=True)
    pval_adj = Column(Float, nullable=True)
    pct_in_cluster = Column(Float, nullable=True)
    pct_out_cluster = Column(Float, nullable=True)

    single_cell_dataset: Mapped["SingleCellDataset"] = relationship(
        "SingleCellDataset", foreign_keys=[single_cell_dataset_id], back_populates="markers"
    )
    feature: Mapped[Optional["Feature"]] = relationship("Feature", foreign_keys=[feature_id])

    __table_args__ = (
        Index("ix_cell_type_markers_dataset_cluster", "single_cell_dataset_id", "cluster_id"),
    )


class SpectralLibrary(Base):
    """Reference spectral library (e.g., MGF) for lipidomics matching."""

    __tablename__ = "spectral_libraries"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(200), nullable=False, index=True)
    version = Column(String(100), nullable=True, index=True)
    source_url = Column(String(500), nullable=True)
    file_path = Column(String(500), nullable=True)
    n_spectra = Column(Integer, nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    references: Mapped[List["SpectralReference"]] = relationship(
        "SpectralReference", back_populates="library", cascade="all, delete-orphan"
    )

    __table_args__ = (
        UniqueConstraint("name", "version", name="uq_spectral_library_name_version"),
    )


class SpectralReference(Base):
    """A single reference spectrum entry from a spectral library."""

    __tablename__ = "spectral_references"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    library_id = Column(UUID(as_uuid=True), ForeignKey("spectral_libraries.id"), nullable=False, index=True)

    lipid_name = Column(String(500), nullable=False, index=True)
    lipid_class = Column(String(200), nullable=True, index=True)
    smiles = Column(Text, nullable=True)
    inchi_key = Column(String(50), nullable=True, index=True)

    precursor_mz = Column(Float, nullable=False, index=True)
    precursor_type = Column(String(50), nullable=True)
    collision_energy = Column(Float, nullable=True)

    spectrum = Column(JSON, nullable=False)  # {"mz": [...], "intensity": [...]}

    library: Mapped["SpectralLibrary"] = relationship("SpectralLibrary", back_populates="references")

    __table_args__ = (
        Index("ix_spectral_ref_library_precursor", "library_id", "precursor_mz"),
    )


class LipidAnnotation(Base):
    """Match between an internal Feature and a SpectralReference."""

    __tablename__ = "lipid_annotations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=False, index=True)
    spectral_reference_id = Column(UUID(as_uuid=True), ForeignKey("spectral_references.id"), nullable=False, index=True)

    spectral_score = Column(Float, nullable=False)
    mz_error_ppm = Column(Float, nullable=True)
    matched_peaks = Column(Integer, nullable=True)
    total_peaks = Column(Integer, nullable=True)

    is_confident = Column(Boolean, nullable=False, default=False)
    is_ambiguous = Column(Boolean, nullable=False, default=False)
    rank = Column(Integer, nullable=True)

    manually_reviewed = Column(Boolean, nullable=False, default=False)
    review_status = Column(String(50), nullable=True)  # confirmed|rejected
    reviewed_at = Column(DateTime(timezone=True), nullable=True)

    feature: Mapped["Feature"] = relationship("Feature", foreign_keys=[feature_id])
    spectral_reference: Mapped["SpectralReference"] = relationship("SpectralReference", foreign_keys=[spectral_reference_id])

    __table_args__ = (
        Index("ix_lipid_annotations_feature_rank", "feature_id", "rank"),
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

    # Lifecycle status
    lifecycle_status = Column(String(20), default='active', nullable=False, index=True)
    retention_exempt = Column(Boolean, default=False, nullable=False)

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


class LINCSGene(Base):
    """LINCS L1000 gene metadata."""

    __tablename__ = "lincs_genes"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    entrez_id = Column(Integer, nullable=False, unique=True, index=True)
    gene_symbol = Column(String(50), nullable=True, index=True)
    gene_title = Column(String(500), nullable=True)
    is_landmark = Column(Boolean, nullable=False, default=False)

    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=True, index=True)
    feature: Mapped[Optional["Feature"]] = relationship("Feature", foreign_keys=[feature_id])


class LINCSSignature(Base):
    """LINCS gene expression signature (typically z-scores for many genes)."""

    __tablename__ = "lincs_signatures"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    sig_id = Column(String(100), nullable=False, unique=True, index=True)

    pert_iname = Column(String(200), nullable=True, index=True)
    pert_id = Column(String(100), nullable=True, index=True)
    pert_type = Column(String(100), nullable=True)

    cell_id = Column(String(100), nullable=True, index=True)
    pert_time = Column(Float, nullable=True)
    pert_dose = Column(Float, nullable=True)

    gene_expression = Column(JSON, nullable=True)  # {entrez_id: z-score}

    tas = Column(Float, nullable=True)
    distil_id = Column(String(500), nullable=True)
    ingested_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)


class ConnectivityScore(Base):
    """Connectivity score between a query Signature and a LINCS signature."""

    __tablename__ = "connectivity_scores"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    query_signature_id = Column(UUID(as_uuid=True), ForeignKey("signatures.id"), nullable=False, index=True)
    lincs_signature_id = Column(UUID(as_uuid=True), ForeignKey("lincs_signatures.id"), nullable=False, index=True)

    score = Column(Float, nullable=False)
    score_type = Column(String(50), nullable=False, default="connectivity")
    p_value = Column(Float, nullable=True)
    computed_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    query_signature: Mapped["Signature"] = relationship("Signature", foreign_keys=[query_signature_id])
    lincs_signature: Mapped["LINCSSignature"] = relationship("LINCSSignature", foreign_keys=[lincs_signature_id])

    __table_args__ = (
        UniqueConstraint(
            "query_signature_id",
            "lincs_signature_id",
            "score_type",
            name="uq_connectivity_score",
        ),
        Index("ix_connectivity_scores_query_score_type", "query_signature_id", "score_type"),
    )


class CRISPRScreen(Base):
    """CRISPR screen run metadata."""

    __tablename__ = "crispr_screens"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    dataset_id = Column(UUID(as_uuid=True), ForeignKey("datasets.id"), nullable=False, index=True)

    name = Column(String(500), nullable=False, index=True)
    library_type = Column(String(100), nullable=True, index=True)
    cell_line = Column(String(200), nullable=True, index=True)
    treatment = Column(String(200), nullable=True, index=True)
    control_label = Column(String(200), nullable=True)
    treatment_label = Column(String(200), nullable=True)
    status = Column(String(50), nullable=True, index=True)  # pending|running|completed|failed

    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    dataset: Mapped["Dataset"] = relationship("Dataset", foreign_keys=[dataset_id])
    guides: Mapped[List["CRISPRGuide"]] = relationship(
        "CRISPRGuide", back_populates="screen", cascade="all, delete-orphan"
    )
    results: Mapped[List["CRISPRResult"]] = relationship(
        "CRISPRResult", back_populates="screen", cascade="all, delete-orphan"
    )


class CRISPRGuide(Base):
    """CRISPR guide information for a screen."""

    __tablename__ = "crispr_guides"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    screen_id = Column(UUID(as_uuid=True), ForeignKey("crispr_screens.id"), nullable=False, index=True)

    guide_seq = Column(String(200), nullable=False, index=True)
    gene_symbol = Column(String(100), nullable=True, index=True)
    gene_id = Column(String(100), nullable=True, index=True)

    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=True, index=True)
    feature: Mapped[Optional["Feature"]] = relationship("Feature", foreign_keys=[feature_id])

    screen: Mapped["CRISPRScreen"] = relationship("CRISPRScreen", back_populates="guides", foreign_keys=[screen_id])

    __table_args__ = (
        UniqueConstraint("screen_id", "guide_seq", name="uq_crispr_guide_screen_seq"),
        Index("ix_crispr_guides_screen_gene", "screen_id", "gene_symbol"),
    )


class CRISPRResult(Base):
    """Per-gene result for a CRISPR screen."""

    __tablename__ = "crispr_results"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    screen_id = Column(UUID(as_uuid=True), ForeignKey("crispr_screens.id"), nullable=False, index=True)

    gene_symbol = Column(String(100), nullable=True, index=True)
    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=True, index=True)

    beta_score = Column(Float, nullable=True)
    p_value = Column(Float, nullable=True)
    fdr = Column(Float, nullable=True)
    neg_lfc = Column(Float, nullable=True)
    pos_lfc = Column(Float, nullable=True)

    rank = Column(Integer, nullable=True, index=True)
    is_hit = Column(Boolean, nullable=False, default=False, index=True)

    screen: Mapped["CRISPRScreen"] = relationship("CRISPRScreen", back_populates="results", foreign_keys=[screen_id])
    feature: Mapped[Optional["Feature"]] = relationship("Feature", foreign_keys=[feature_id])

    __table_args__ = (
        Index("ix_crispr_results_screen_gene", "screen_id", "gene_symbol"),
    )


class MultiOmicsExperiment(Base):
    """Multi-omics experiment for latent factor analysis (e.g., MOFA-style)."""

    __tablename__ = "multi_omics_experiments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(500), nullable=False, index=True)
    description = Column(Text, nullable=True)

    dataset_ids = Column(JSON, nullable=True)  # [dataset_id, ...]
    sample_mapping = Column(JSON, nullable=True)  # mapping across datasets -> unified sample ids

    n_factors = Column(Integer, nullable=True)
    convergence_mode = Column(String(100), nullable=True)

    status = Column(String(50), nullable=True, index=True)  # pending|running|completed|failed
    processing_log = Column(Text, nullable=True)

    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    processed_at = Column(DateTime, nullable=True)

    factors: Mapped[List["LatentFactor"]] = relationship(
        "LatentFactor", back_populates="experiment", cascade="all, delete-orphan"
    )


class LatentFactor(Base):
    """A single latent factor within a MultiOmicsExperiment."""

    __tablename__ = "latent_factors"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(
        UUID(as_uuid=True), ForeignKey("multi_omics_experiments.id"), nullable=False, index=True
    )
    factor_index = Column(Integer, nullable=False, index=True)

    variance_explained = Column(JSON, nullable=True)  # {"omics_type": float, ...}
    description = Column(Text, nullable=True)

    experiment: Mapped["MultiOmicsExperiment"] = relationship(
        "MultiOmicsExperiment", back_populates="factors", foreign_keys=[experiment_id]
    )
    loadings: Mapped[List["FactorLoading"]] = relationship(
        "FactorLoading", back_populates="factor", cascade="all, delete-orphan"
    )
    scores: Mapped[List["FactorScore"]] = relationship(
        "FactorScore", back_populates="factor", cascade="all, delete-orphan"
    )

    __table_args__ = (
        UniqueConstraint("experiment_id", "factor_index", name="uq_latent_factor_experiment_index"),
    )


class FactorLoading(Base):
    """Feature loading for a specific latent factor."""

    __tablename__ = "factor_loadings"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    factor_id = Column(UUID(as_uuid=True), ForeignKey("latent_factors.id"), nullable=False, index=True)
    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=False, index=True)

    loading = Column(Float, nullable=False)
    abs_loading = Column(Float, nullable=False, index=True)
    omics_type = Column(String(50), nullable=True, index=True)

    factor: Mapped["LatentFactor"] = relationship("LatentFactor", back_populates="loadings", foreign_keys=[factor_id])
    feature: Mapped["Feature"] = relationship("Feature", foreign_keys=[feature_id])

    __table_args__ = (
        UniqueConstraint("factor_id", "feature_id", name="uq_factor_loading_factor_feature"),
        Index("ix_factor_loadings_factor_abs", "factor_id", "abs_loading"),
    )


class FactorScore(Base):
    """Per-sample score for a specific latent factor."""

    __tablename__ = "factor_scores"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    factor_id = Column(UUID(as_uuid=True), ForeignKey("latent_factors.id"), nullable=False, index=True)
    sample_id = Column(String(200), nullable=False, index=True)
    score = Column(Float, nullable=False)

    factor: Mapped["LatentFactor"] = relationship("LatentFactor", back_populates="scores", foreign_keys=[factor_id])

    __table_args__ = (
        UniqueConstraint("factor_id", "sample_id", name="uq_factor_score_factor_sample"),
        Index("ix_factor_scores_factor_sample", "factor_id", "sample_id"),
    )


class VariantSet(Base):
    """A set of genetic variants imported from a source file (VCF/TSV/etc.)."""

    __tablename__ = "variant_sets"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(500), nullable=False, index=True)
    description = Column(Text, nullable=True)

    source_file = Column(String(500), nullable=True)
    source_type = Column(String(100), nullable=True, index=True)  # vcf/tsv/clinvar/etc.

    n_variants = Column(Integer, nullable=True)
    n_genes = Column(Integer, nullable=True)

    status = Column(String(50), nullable=True, index=True)  # pending|processing|completed|failed
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    variants: Mapped[List["Variant"]] = relationship(
        "Variant", back_populates="variant_set", cascade="all, delete-orphan"
    )
    burdens: Mapped[List["GeneBurden"]] = relationship(
        "GeneBurden", back_populates="variant_set", cascade="all, delete-orphan"
    )


class Variant(Base):
    """A single variant record in a VariantSet."""

    __tablename__ = "variants"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    variant_set_id = Column(UUID(as_uuid=True), ForeignKey("variant_sets.id"), nullable=False, index=True)

    chromosome = Column(String(20), nullable=False, index=True)
    position = Column(Integer, nullable=False, index=True)
    ref_allele = Column(String(500), nullable=False)
    alt_allele = Column(String(500), nullable=False)

    rs_id = Column(String(50), nullable=True, index=True)
    hgvs_genomic = Column(String(500), nullable=True, index=True)
    hgvs_coding = Column(String(500), nullable=True, index=True)
    hgvs_protein = Column(String(500), nullable=True, index=True)

    gene_symbol = Column(String(100), nullable=True, index=True)
    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=True, index=True)

    consequence = Column(String(200), nullable=True, index=True)
    impact = Column(String(50), nullable=True, index=True)
    gnomad_af = Column(Float, nullable=True)

    variant_set: Mapped["VariantSet"] = relationship("VariantSet", back_populates="variants", foreign_keys=[variant_set_id])
    feature: Mapped[Optional["Feature"]] = relationship("Feature", foreign_keys=[feature_id])
    annotations: Mapped[List["VariantAnnotation"]] = relationship(
        "VariantAnnotation", back_populates="variant", cascade="all, delete-orphan"
    )

    __table_args__ = (
        UniqueConstraint(
            "variant_set_id",
            "chromosome",
            "position",
            "ref_allele",
            "alt_allele",
            name="uq_variant_set_locus",
        ),
        Index("ix_variants_set_gene", "variant_set_id", "gene_symbol"),
    )


class VariantAnnotation(Base):
    """Annotation for a variant from an external source (ClinVar/VEP/CADD/etc.)."""

    __tablename__ = "variant_annotations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    variant_id = Column(UUID(as_uuid=True), ForeignKey("variants.id"), nullable=False, index=True)

    annotation_source = Column(String(100), nullable=False, index=True)
    clinvar_id = Column(String(100), nullable=True, index=True)
    clinical_significance = Column(String(200), nullable=True, index=True)
    review_status = Column(String(200), nullable=True)
    condition = Column(String(500), nullable=True)

    cadd_phred = Column(Float, nullable=True)
    revel_score = Column(Float, nullable=True)
    sift_prediction = Column(String(50), nullable=True)
    polyphen_prediction = Column(String(50), nullable=True)
    
    # Additional VEP fields
    consequence = Column(String(100), nullable=True, index=True)  # e.g., missense_variant
    impact = Column(String(20), nullable=True, index=True)  # HIGH, MODERATE, LOW, MODIFIER
    symbol = Column(String(50), nullable=True, index=True)  # Gene symbol
    gene_id = Column(String(50), nullable=True)  # Ensembl gene ID
    transcript_id = Column(String(50), nullable=True)  # Ensembl transcript ID
    biotype = Column(String(50), nullable=True)  # protein_coding, etc.
    
    # Protein change
    amino_acids = Column(String(20), nullable=True)  # e.g., "A/T"
    codons = Column(String(20), nullable=True)  # e.g., "Gcc/Acc"
    protein_position = Column(Integer, nullable=True)
    
    # Additional predictions
    sift_score = Column(Float, nullable=True)
    polyphen_score = Column(Float, nullable=True)
    
    # Clinical (additional to existing clinical_significance)
    clin_sig = Column(String(100), nullable=True)  # pathogenic, benign, etc.
    pubmed_ids = Column(ARRAY(String), nullable=True)
    
    # Metadata
    source_version = Column(String(50), nullable=True)
    annotated_at = Column(DateTime(timezone=True), server_default=func.now())

    variant: Mapped["Variant"] = relationship("Variant", back_populates="annotations", foreign_keys=[variant_id])

    __table_args__ = (
        Index("ix_variant_annotations_variant_source", "variant_id", "annotation_source"),
    )


class GeneBurden(Base):
    """Aggregated per-gene burden metrics for a VariantSet."""

    __tablename__ = "gene_burdens"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    variant_set_id = Column(UUID(as_uuid=True), ForeignKey("variant_sets.id"), nullable=False, index=True)

    gene_symbol = Column(String(100), nullable=True, index=True)
    feature_id = Column(UUID(as_uuid=True), ForeignKey("features.id"), nullable=True, index=True)

    n_variants = Column(Integer, nullable=True)
    n_pathogenic = Column(Integer, nullable=True)
    n_vus = Column(Integer, nullable=True)
    n_benign = Column(Integer, nullable=True)
    burden_score = Column(Float, nullable=True)

    variant_set: Mapped["VariantSet"] = relationship("VariantSet", back_populates="burdens", foreign_keys=[variant_set_id])
    feature: Mapped[Optional["Feature"]] = relationship("Feature", foreign_keys=[feature_id])

    __table_args__ = (
        UniqueConstraint("variant_set_id", "gene_symbol", name="uq_gene_burden_set_gene"),
        Index("ix_gene_burdens_set_gene", "variant_set_id", "gene_symbol"),
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


class DigestSchedule(Base):
    """Scheduled execution of notebook digests (Papermill -> HTML render)."""

    __tablename__ = "digest_schedules"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    program_id = Column(UUID(as_uuid=True), ForeignKey("programs.id", ondelete="CASCADE"), nullable=False, index=True)
    notebook_path = Column(String(500), nullable=False)
    schedule_cron = Column(String(100), nullable=False)
    recipients = Column(JSONB, nullable=False)  # list of emails or user ids
    last_run_at = Column(DateTime(timezone=True), nullable=True)
    last_status = Column(String(20), nullable=True)  # success|failed
    enabled = Column(Boolean, default=True, nullable=False)

    program: Mapped["Program"] = relationship("Program")


class NotebookReview(Base):
    """Notebook review record with signed decision and audit trail."""

    __tablename__ = "notebook_reviews"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    notebook_path = Column(String(500), nullable=False, index=True)
    version_hash = Column(String(64), nullable=False, index=True)
    reviewer_id = Column(UUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"), nullable=True)
    status = Column(String(32), nullable=False, default="pending", index=True)
    comments = Column(Text, nullable=True)
    signature = Column(String(128), nullable=True)
    reviewed_at = Column(DateTime(timezone=True), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    reviewer: Mapped[Optional["User"]] = relationship("User", foreign_keys=[reviewer_id])
    threads: Mapped[list["ReviewThread"]] = relationship("ReviewThread", back_populates="review", cascade="all, delete-orphan")
    snapshot: Mapped[Optional["NotebookSnapshot"]] = relationship("NotebookSnapshot", back_populates="review", uselist=False, cascade="all, delete-orphan")

    __table_args__ = (Index("ix_notebook_reviews_notebook_status", "notebook_path", "status"),)


class ReviewThread(Base):
    """Discussion thread anchored to a notebook review, optionally to a specific cell."""

    __tablename__ = "review_threads"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    review_id = Column(UUID(as_uuid=True), ForeignKey("notebook_reviews.id", ondelete="CASCADE"), nullable=False)
    cell_index = Column(Integer, nullable=True)  # null = general, int = cell-specific
    title = Column(String(255), nullable=False)
    status = Column(String(20), nullable=False, default="open")  # open/resolved/wontfix
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), onupdate=lambda: datetime.now(timezone.utc), nullable=False)

    # Relationships
    review: Mapped["NotebookReview"] = relationship("NotebookReview", back_populates="threads")
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])
    comments: Mapped[list["ReviewComment"]] = relationship("ReviewComment", back_populates="thread", cascade="all, delete-orphan")

    __table_args__ = (Index("ix_review_threads_review_status", "review_id", "status"),)


class ReviewComment(Base):
    """Individual comment in a review thread with support for nested replies."""

    __tablename__ = "review_comments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    thread_id = Column(UUID(as_uuid=True), ForeignKey("review_threads.id", ondelete="CASCADE"), nullable=False)
    parent_id = Column(UUID(as_uuid=True), ForeignKey("review_comments.id", ondelete="SET NULL"), nullable=True)
    content = Column(Text, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), onupdate=lambda: datetime.now(timezone.utc), nullable=False)

    # Relationships
    thread: Mapped["ReviewThread"] = relationship("ReviewThread", back_populates="comments")
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])
    replies: Mapped[list["ReviewComment"]] = relationship("ReviewComment", back_populates="parent", cascade="all, delete-orphan")
    parent: Mapped[Optional["ReviewComment"]] = relationship("ReviewComment", back_populates="replies", remote_side=[id])


class NotebookSnapshot(Base):
    """Stores notebook content at review time for diff computation."""

    __tablename__ = "notebook_snapshots"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    review_id = Column(UUID(as_uuid=True), ForeignKey("notebook_reviews.id", ondelete="CASCADE"), nullable=False, unique=True)
    content_json = Column(JSON, nullable=False)  # Full notebook content
    cell_count = Column(Integer, nullable=False, default=0)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    # Relationships
    review: Mapped["NotebookReview"] = relationship("NotebookReview", back_populates="snapshot")


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


class RepositoryNotification(Base):
    """Notification generated from activity events or repository subscription matches."""

    __tablename__ = "repository_notifications"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    subscription_id = Column(UUID(as_uuid=True), ForeignKey("repository_subscriptions.id"), nullable=True)
    dataset_id = Column(UUID(as_uuid=True), ForeignKey("datasets.id"), nullable=True)
    activity_event_id = Column(UUID(as_uuid=True), ForeignKey("activity_events.id"), nullable=True)
    recipient_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True, index=True)
    notification_type = Column(String(50), default="discovery", nullable=False)
    is_read = Column(Boolean, default=False, nullable=False)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)

    subscription: Mapped[Optional["RepositorySubscription"]] = relationship("RepositorySubscription")
    dataset: Mapped[Optional["Dataset"]] = relationship("Dataset")
    activity_event: Mapped["ActivityEvent"] = relationship("ActivityEvent", back_populates="repository_notifications")
    recipient: Mapped[Optional["User"]] = relationship("User", foreign_keys=[recipient_id])


class ActivityEvent(Base):
    """Activity event that tracks user actions and system events."""

    __tablename__ = "activity_events"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    event_type = Column(String(50), nullable=False, index=True)
    actor_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True, index=True)
    target_type = Column(String(50), nullable=True, index=True)
    target_id = Column(UUID(as_uuid=True), nullable=True, index=True)
    target_name = Column(String(255), nullable=True)  # snapshot for orphan handling
    program_id = Column(UUID(as_uuid=True), ForeignKey("programs.id"), nullable=True, index=True)
    event_metadata = Column(JSON, nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False, index=True)

    # Relationships
    actor: Mapped["User"] = relationship("User")
    program: Mapped["Program"] = relationship("Program")
    repository_notifications: Mapped[List["RepositoryNotification"]] = relationship("RepositoryNotification", back_populates="activity_event")


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


class ExtractionJob(Base):
    """Batch document extraction job for file-based ingestion."""

    __tablename__ = "extraction_jobs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    batch_id = Column(UUID(as_uuid=True), nullable=True, index=True)

    file_count = Column(Integer, nullable=False)
    completed_count = Column(Integer, nullable=False, default=0)
    status = Column(String(50), nullable=False, default="pending")  # pending|running|completed|failed

    started_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    completed_at = Column(DateTime(timezone=True), nullable=True)

    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        onupdate=lambda: datetime.now(timezone.utc),
        nullable=False,
    )

    documents: Mapped[List["ExtractedDocument"]] = relationship(
        "ExtractedDocument", back_populates="job", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_extraction_jobs_status", "status"),
        Index("ix_extraction_jobs_batch_status", "batch_id", "status"),
    )


class ExtractedDocument(Base):
    """Extracted document metadata for a given ExtractionJob."""

    __tablename__ = "extracted_documents"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    job_id = Column(UUID(as_uuid=True), ForeignKey("extraction_jobs.id"), nullable=False, index=True)

    file_path = Column(String(500), nullable=False)
    original_filename = Column(String(500), nullable=False)
    doc_type = Column(String(50), nullable=False, index=True)  # pdf|docx|pptx|csv|xlsx

    extracted_entities = Column(JSONB, nullable=False)
    extraction_config = Column(JSONB, nullable=True)

    status = Column(String(50), nullable=False, default="pending", index=True)  # pending|completed|failed
    error_log = Column(Text, nullable=True)

    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        onupdate=lambda: datetime.now(timezone.utc),
        nullable=False,
    )

    job: Mapped["ExtractionJob"] = relationship("ExtractionJob", back_populates="documents", foreign_keys=[job_id])

    __table_args__ = (
        Index("ix_extracted_documents_job_status", "job_id", "status"),
        Index("ix_extracted_documents_job_doc_type", "job_id", "doc_type"),
    )


class SyncJob(Base):
    """A synchronization job against an external source (ChEMBL/PubChem/UniProt/KEGG)."""

    __tablename__ = "sync_jobs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)

    source = Column(String(50), nullable=False, index=True)  # chembl/pubchem/uniprot/kegg
    sync_type = Column(String(50), nullable=False)  # full/incremental
    status = Column(String(50), nullable=False, default="pending", index=True)  # pending/running/completed/failed

    records_synced = Column(Integer, nullable=False, default=0)
    records_updated = Column(Integer, nullable=False, default=0)
    records_new = Column(Integer, nullable=False, default=0)
    conflicts_detected = Column(Integer, nullable=False, default=0)

    started_at = Column(DateTime(timezone=True), nullable=True)
    completed_at = Column(DateTime(timezone=True), nullable=True)
    error_log = Column(Text, nullable=True)

    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        onupdate=lambda: datetime.now(timezone.utc),
        nullable=False,
    )

    records: Mapped[List["SyncRecord"]] = relationship(
        "SyncRecord", back_populates="job", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_sync_jobs_source_status", "source", "status"),
    )


class SyncRecord(Base):
    """A unique mapping of a source external_id to an internal entity (best-effort)."""

    __tablename__ = "sync_records"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)

    job_id = Column(UUID(as_uuid=True), ForeignKey("sync_jobs.id"), nullable=True, index=True)
    source = Column(String(50), nullable=False, index=True)
    external_id = Column(String(200), nullable=False, index=True)  # e.g. CHEMBL123456

    entity_type = Column(String(50), nullable=False)  # compound/activity/assay/target
    entity_id = Column(UUID(as_uuid=True), nullable=True)

    checksum = Column(String(32), nullable=False)  # MD5 hash
    synced_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    # NOTE: `metadata` is reserved on SQLAlchemy declarative models; use `metadata_` attribute name.
    metadata_ = Column("metadata", JSONB, nullable=False)

    job: Mapped[Optional["SyncJob"]] = relationship("SyncJob", back_populates="records", foreign_keys=[job_id])
    conflicts: Mapped[List["SyncConflict"]] = relationship(
        "SyncConflict", back_populates="record", cascade="all, delete-orphan"
    )

    __table_args__ = (
        UniqueConstraint("source", "external_id", name="uq_sync_record_source_external_id"),
        Index("ix_sync_records_source_external_id", "source", "external_id"),
    )


class SyncConflict(Base):
    """A conflict detected during sync (e.g., local edits or schema mismatch)."""

    __tablename__ = "sync_conflicts"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    record_id = Column(UUID(as_uuid=True), ForeignKey("sync_records.id"), nullable=False, index=True)

    conflict_type = Column(String(50), nullable=False)  # local_edit/schema_mismatch/duplicate
    local_value = Column(JSONB, nullable=False)
    external_value = Column(JSONB, nullable=False)

    resolution_status = Column(
        String(50),
        nullable=False,
        default="pending",  # pending/auto_merged/manual_override/ignored
        index=True,
    )
    resolved_at = Column(DateTime(timezone=True), nullable=True)

    record: Mapped["SyncRecord"] = relationship("SyncRecord", back_populates="conflicts", foreign_keys=[record_id])

    __table_args__ = (
        Index("ix_sync_conflicts_record_status", "record_id", "resolution_status"),
    )


class ChEMBLActivity(Base):
    """ChEMBL activity measurement linked to an optional internal Compound."""

    __tablename__ = "chembl_activities"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=True, index=True)

    chembl_molecule_id = Column(String(50), nullable=False, index=True)
    chembl_assay_id = Column(String(50), nullable=False, index=True)

    activity_type = Column(String(50), nullable=True)  # IC50/Ki/Kd/EC50
    value = Column(Float, nullable=True)
    units = Column(String(50), nullable=True)
    relation = Column(String(10), nullable=True)

    target_chembl_id = Column(String(50), nullable=True, index=True)
    target_name = Column(String(500), nullable=True)
    assay_type = Column(String(50), nullable=True)  # binding/functional

    synced_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    compound: Mapped[Optional["Compound"]] = relationship("Compound", foreign_keys=[compound_id])

    __table_args__ = (
        Index("ix_chembl_activities_molecule_assay", "chembl_molecule_id", "chembl_assay_id"),
    )


class PubChemBioassay(Base):
    """PubChem bioassay record linked to an optional internal Compound."""

    __tablename__ = "pubchem_bioassays"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=True, index=True)

    pubchem_cid = Column(Integer, nullable=False, index=True)
    aid = Column(Integer, nullable=False, index=True)

    activity_outcome = Column(String(50), nullable=True)  # active/inactive/inconclusive
    activity_score = Column(Float, nullable=True)
    assay_name = Column(String(500), nullable=True)

    synced_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    compound: Mapped[Optional["Compound"]] = relationship("Compound", foreign_keys=[compound_id])

    __table_args__ = (
        Index("ix_pubchem_bioassays_cid_aid", "pubchem_cid", "aid"),
    )


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

    NOTE: GraphNode table enhancement tracked in ROADMAP - would provide node-level
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


class IDMapping(Base):
    """Cached ID mappings from external sources (UniProt, KEGG, etc.)."""
    
    __tablename__ = "id_mappings"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    source_type = Column(String(50), nullable=False)  # gene, protein, metabolite
    source_id = Column(String(200), nullable=False)
    target_type = Column(String(50), nullable=False)  # uniprot, kegg_gene, kegg_compound
    target_id = Column(String(200), nullable=False)
    organism = Column(String(50), default="human")
    confidence = Column(Float, nullable=True)
    expires_at = Column(DateTime(timezone=True), nullable=True)  # For TTL-based caching
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc))
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), onupdate=lambda: datetime.now(timezone.utc))
    
    __table_args__ = (
        UniqueConstraint("source_type", "source_id", "target_type", "organism", name="uq_id_mapping_source_target"),
        Index("ix_id_mappings_lookup", "source_type", "source_id", "target_type"),
    )


class MappingRefreshLog(Base):
    """Log of ID mapping refresh operations for tracking sync timestamps."""
    
    __tablename__ = "mapping_refresh_logs"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    source = Column(String(50), nullable=False)  # uniprot, kegg
    status = Column(String(20), nullable=False)  # success, failed
    records_processed = Column(Integer, default=0)
    started_at = Column(DateTime(timezone=True), nullable=False)
    completed_at = Column(DateTime(timezone=True), nullable=True)
    error_message = Column(Text, nullable=True)
    
    __table_args__ = (
        Index("ix_mapping_refresh_logs_source_status", "source", "status"),
    )


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
    Company,
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
)
from amprenta_rag.database.models_flow_cytometry import (  # noqa: F401
    FlowCytometryDataset,
    FlowCytometryParameter,
    FlowCytometryGate,
    FlowCytometryPopulation,
)
from amprenta_rag.database.models_biophysical import (  # noqa: F401
    SPRExperiment,
    SPRSensorgram,
    MSTExperiment,
    MSTDoseResponse,
    DSCExperiment,
    DSCScan,
)
from amprenta_rag.models.eln import (  # noqa: F401
    Protocol,
    ExperimentProtocol,
    LabNotebookEntry,
    LabNotebookEntryAssociation,
)


class ModelMonitoringLog(Base):
    """Model monitoring log for tracking predictions and performance."""

    __tablename__ = "model_monitoring_logs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    model_id = Column(UUID(as_uuid=True), ForeignKey("ml_models.id"), index=True, nullable=False)
    model_version = Column(String(50), index=True, nullable=False)
    prediction_id = Column(UUID(as_uuid=True), unique=True, index=True, nullable=False)
    prediction = Column(Float, nullable=False)
    ground_truth = Column(Float, nullable=True)
    input_hash = Column(String(64), index=True)
    feature_summary = Column(JSON, nullable=False)
    created_at = Column(DateTime(timezone=True), default=datetime.utcnow, index=True)


class BackupRecord(Base):
    """Database backup record for tracking backup operations and metadata."""

    __tablename__ = "backup_records"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    backup_type = Column(String(20), nullable=False, index=True)  # 'full' or 'incremental'
    status = Column(String(20), nullable=False, index=True)  # 'pending', 'running', 'completed', 'failed'
    file_path = Column(String(500), nullable=True)  # S3 key or local file path
    file_size_bytes = Column(BigInteger, nullable=True)
    checksum_sha256 = Column(String(64), nullable=True, index=True)
    started_at = Column(DateTime(timezone=True), nullable=True, index=True)
    completed_at = Column(DateTime(timezone=True), nullable=True, index=True)
    error_message = Column(Text, nullable=True)
    backup_metadata = Column(JSONB, nullable=True)  # JSON: tables, row_counts, pg_version, etc.
    created_at = Column(DateTime(timezone=True), default=datetime.utcnow, nullable=False, index=True)
    updated_at = Column(DateTime(timezone=True), default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Indexes for common queries
    __table_args__ = (
        Index("ix_backup_records_status_created", "status", "created_at"),
        Index("ix_backup_records_type_status", "backup_type", "status"),
    )


class ProjectExport(Base):
    """Temporary project export storage for download."""
    
    __tablename__ = "project_exports"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    file_path = Column(String, nullable=False)  # S3 key or local path
    file_size_bytes = Column(BigInteger, nullable=True)
    checksum_sha256 = Column(String(64), nullable=True)
    entities_summary = Column(String, nullable=True)  # JSON summary
    created_by = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    expires_at = Column(DateTime(timezone=True), nullable=False)
    created_at = Column(DateTime(timezone=True), default=datetime.utcnow, nullable=False, index=True)
    
    # Relationship to user
    creator = relationship("User", backref="project_exports")
    
    # Index for cleanup queries
    __table_args__ = (
        Index("ix_project_exports_expires_at", "expires_at"),
    )


class EntityVersion(Base):
    """Version snapshot for entity provenance tracking."""
    
    __tablename__ = "entity_versions"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    entity_type = Column(String(50), nullable=False, index=True)
    entity_id = Column(UUID(as_uuid=True), nullable=False, index=True)
    version_number = Column(Integer, nullable=False)
    data_snapshot = Column(JSONB, nullable=False)
    checksum_sha256 = Column(String(64), nullable=False)
    created_by = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    change_summary = Column(String(500), nullable=True)
    
    creator = relationship("User", foreign_keys=[created_by])
    
    __table_args__ = (
        UniqueConstraint("entity_type", "entity_id", "version_number", name="uq_entity_version"),
        Index("ix_entity_versions_lookup", "entity_type", "entity_id"),
    )


class ReviewCycle(Base):
    """Defines recurring review schedules for entities."""
    
    __tablename__ = "review_cycles"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(100), nullable=False)
    entity_type = Column(String(50), nullable=False)  # dataset/experiment/compound/signature
    frequency = Column(String(20), nullable=False)  # weekly/monthly/quarterly/yearly
    day_of_week = Column(Integer, nullable=True)  # 0=Mon, 6=Sun (for weekly)
    day_of_month = Column(Integer, nullable=True)  # 1-28 (for monthly/quarterly)
    next_run_at = Column(DateTime(timezone=True), nullable=True)
    reviewer_pool = Column(JSON, nullable=True)  # List of reviewer UUIDs
    is_active = Column(Boolean, default=True, nullable=False)
    program_id = Column(UUID(as_uuid=True), ForeignKey("programs.id", ondelete="SET NULL"), nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    
    program: Mapped[Optional["Program"]] = relationship("Program")
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])
    
    __table_args__ = (
        Index("ix_review_cycles_entity_type", "entity_type"),
        Index("ix_review_cycles_next_run", "next_run_at"),
    )


class ReviewSLA(Base):
    """Defines SLA rules for reviews with escalation chains."""
    
    __tablename__ = "review_slas"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(100), nullable=False)
    entity_type = Column(String(50), nullable=True)  # null = applies to all
    max_review_hours = Column(Integer, nullable=False, default=120)  # 5 days default
    warning_threshold_pct = Column(Integer, nullable=False, default=75)
    escalation_chain = Column(JSON, nullable=True)  # List of user UUIDs
    is_default = Column(Boolean, default=False, nullable=False)
    is_active = Column(Boolean, default=True, nullable=False)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    
    __table_args__ = (
        Index("ix_review_slas_entity_type", "entity_type"),
    )
