"""
SQLAlchemy models for the multi-omics platform.

These models map domain models to Postgres tables. They maintain compatibility
with Notion IDs during migration by storing notion_page_id for dual-write support.
"""

from __future__ import annotations

import uuid
from datetime import datetime

from sqlalchemy import JSON, Boolean, Column, DateTime, Float, ForeignKey, Integer, String, Table, Text, func
from sqlalchemy.dialects.postgresql import ARRAY, UUID, TSVECTOR
from sqlalchemy.orm import relationship

from amprenta_rag.database.base import Base


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
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # External identifiers (stored as JSON for flexibility)
    external_ids = Column(JSON, nullable=True)

    # Audit
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_by = relationship("User", foreign_keys=[created_by_id])

    # Relationships
    experiments = relationship("Experiment", secondary=program_experiment_assoc, back_populates="programs")
    datasets = relationship("Dataset", secondary=program_dataset_assoc, back_populates="programs")
    signatures = relationship("Signature", secondary=program_signature_assoc, back_populates="programs")


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
    created_by = relationship("User", foreign_keys=[created_by_id])

    # Project
    project_id = Column(UUID(as_uuid=True), ForeignKey("projects.id"), nullable=True)

    # Relationships
    programs = relationship("Program", secondary=program_experiment_assoc, back_populates="experiments")
    datasets = relationship("Dataset", secondary=experiment_dataset_assoc, back_populates="experiments")


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
    datasets = relationship("Dataset", secondary=dataset_feature_assoc, back_populates="features")
    signatures = relationship("Signature", secondary=signature_feature_assoc, back_populates="features")


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
    created_by = relationship("User", foreign_keys=[created_by_id])

    # Relationships
    programs = relationship("Program", secondary=program_dataset_assoc, back_populates="datasets")
    experiments = relationship("Experiment", secondary=experiment_dataset_assoc, back_populates="datasets")
    features = relationship("Feature", secondary=dataset_feature_assoc, back_populates="datasets")
    signatures = relationship("Signature", secondary=dataset_signature_assoc, back_populates="datasets")


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
    signature = relationship("Signature", back_populates="components")
    feature = relationship("Feature")


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

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # Audit
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_by = relationship("User", foreign_keys=[created_by_id])

    # Relationships
    components = relationship("SignatureComponent", back_populates="signature", cascade="all, delete-orphan")
    features = relationship("Feature", secondary=signature_feature_assoc, back_populates="signatures")
    datasets = relationship("Dataset", secondary=dataset_signature_assoc, back_populates="signatures")
    programs = relationship("Program", secondary=program_signature_assoc, back_populates="signatures")


class Literature(Base):
    """Literature item (from Zotero or other sources)."""

    __tablename__ = "literature"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    title = Column(String(1000), nullable=False, index=True)
    source_type = Column(String(100), nullable=True)  # Paper, Book, etc.
    abstract = Column(Text, nullable=True)
    doi = Column(String(500), nullable=True, index=True)
    url = Column(String(1000), nullable=True)
    date = Column(String(100), nullable=True)
    journal = Column(String(500), nullable=True)
    year = Column(Integer, nullable=True)
    tags = Column(ARRAY(String), nullable=True)

    # Zotero-specific fields
    zotero_item_key = Column(String(50), nullable=True, unique=True, index=True)
    zotero_item_type = Column(String(100), nullable=True)

    # Embedding status
    embedding_status = Column(String(50), nullable=True, default="Not Embedded")
    last_ingested_at = Column(DateTime, nullable=True)

    # Semantic metadata (diseases, targets, modality, lipid signatures, etc.)
    semantic_metadata = Column(JSON, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # External identifiers
    external_ids = Column(JSON, nullable=True)

    # Relationships
    chunks = relationship("RAGChunk", back_populates="literature", cascade="all, delete-orphan")


class Email(Base):
    """Email or Note item."""

    __tablename__ = "emails"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    title = Column(String(1000), nullable=False, index=True)
    from_sender = Column(String(500), nullable=True)
    item_type = Column(String(100), nullable=True)  # Email, Note, etc.
    tags = Column(ARRAY(String), nullable=True)

    # Email content (full text)
    content = Column(Text, nullable=True)

    # Semantic metadata (stored as JSON for flexibility)
    semantic_metadata = Column(JSON, nullable=True)  # Diseases, targets, modality, etc.

    # Embedding status
    embedding_status = Column(String(50), nullable=True, default="Not Embedded")
    last_ingested_at = Column(DateTime, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID (for backward compatibility only)
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # External identifiers
    external_ids = Column(JSON, nullable=True)

    # Relationships
    chunks = relationship("RAGChunk", back_populates="email", cascade="all, delete-orphan")


class RAGChunk(Base):
    """RAG chunk from literature, datasets, or other sources."""

    __tablename__ = "rag_chunks"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    chunk_id = Column(String(500), nullable=False, unique=True, index=True)  # Original chunk ID format
    chunk_text = Column(Text, nullable=False)
    chunk_index = Column(Integer, nullable=False, default=0)
    snippet = Column(String(500), nullable=True)  # Short preview

    # Source information
    source_type = Column(String(100), nullable=False)  # Literature, Dataset, Email, etc.
    source_id = Column(UUID(as_uuid=True), nullable=True, index=True)  # UUID of source (literature, dataset, etc.)
    source_name = Column(String(500), nullable=True)

    # Zotero-specific fields
    zotero_item_key = Column(String(50), nullable=True, index=True)
    note_key = Column(String(50), nullable=True)
    attachment_key = Column(String(50), nullable=True)
    note_hash = Column(String(32), nullable=True)  # MD5 hash
    attachment_hash = Column(String(32), nullable=True)  # MD5 hash

    # Metadata (stored as JSON for flexibility)
    chunk_metadata = Column(JSON, nullable=True)  # Renamed from 'metadata' (reserved in SQLAlchemy)
    search_vector = Column(TSVECTOR, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # Relationships
    literature_id = Column(UUID(as_uuid=True), ForeignKey("literature.id"), nullable=True, index=True)
    literature = relationship("Literature", back_populates="chunks")

    email_id = Column(UUID(as_uuid=True), ForeignKey("emails.id"), nullable=True, index=True)
    email = relationship("Email", back_populates="chunks")


class Compound(Base):
    """Chemical compound."""

    __tablename__ = "compounds"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(String(200), nullable=False, unique=True, index=True)  # Internal compound ID
    smiles = Column(Text, nullable=False, index=True)
    inchi_key = Column(String(50), nullable=True, unique=True, index=True)
    canonical_smiles = Column(Text, nullable=True)
    molecular_formula = Column(String(100), nullable=True)
    molecular_weight = Column(Float, nullable=True)
    logp = Column(Float, nullable=True)
    hbd_count = Column(Integer, nullable=True)  # Hydrogen bond donor count
    hba_count = Column(Integer, nullable=True)  # Hydrogen bond acceptor count
    rotatable_bonds = Column(Integer, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # External identifiers
    external_ids = Column(JSON, nullable=True)

    # Relationships
    hts_campaigns = relationship("HTSCampaign", back_populates="library")
    hts_results = relationship("HTSResult", back_populates="compound")
    biochemical_results = relationship("BiochemicalResult", back_populates="compound")
    programs = relationship("Program", secondary="compound_program", back_populates="compounds")


class BiochemicalAssay(Base):
    """Biochemical assay definition."""

    __tablename__ = "biochemical_assays"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(200), nullable=False)
    target = Column(String(200), nullable=True)  # e.g., "EGFR", "CDK4"
    assay_type = Column(String(50), nullable=True)  # e.g., "IC50", "EC50", "Ki"
    unit = Column(String(20), nullable=True)  # e.g., "nM", "uM"
    description = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    results = relationship("ActivityResult", back_populates="assay")


class ActivityResult(Base):
    """Activity measurement for a compound in an assay."""

    __tablename__ = "activity_results"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=False)
    assay_id = Column(UUID(as_uuid=True), ForeignKey("biochemical_assays.id"), nullable=False)
    value = Column(Float, nullable=False)  # The measured value
    qualifier = Column(String(5), nullable=True)  # ">", "<", "=", "~"
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)

    # Relationships
    compound = relationship("Compound", backref="activity_results")
    assay = relationship("BiochemicalAssay", back_populates="results")
    created_by = relationship("User")


class ADMEResult(Base):
    """ADME assay result for a compound."""

    __tablename__ = "adme_results"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=False)
    assay_type = Column(String(50), nullable=False)  # permeability, stability, cyp_inhibition
    value = Column(Float, nullable=False)
    unit = Column(String(20), nullable=True)
    conditions = Column(JSON, nullable=True)  # e.g., {"species": "human", "matrix": "microsomes"}
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)

    compound = relationship("Compound", backref="adme_results")
    created_by = relationship("User")


class PKStudy(Base):
    """Pharmacokinetic study result."""

    __tablename__ = "pk_studies"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=False)
    species = Column(String(50), nullable=False)  # mouse, rat, dog, human
    route = Column(String(50), nullable=True)  # IV, PO, IP
    dose = Column(Float, nullable=True)
    dose_unit = Column(String(20), nullable=True)
    auc = Column(Float, nullable=True)  # Area under curve
    c_max = Column(Float, nullable=True)  # Max concentration
    t_max = Column(Float, nullable=True)  # Time to max
    half_life = Column(Float, nullable=True)  # t1/2
    bioavailability = Column(Float, nullable=True)  # F%
    clearance = Column(Float, nullable=True)  # CL
    vd = Column(Float, nullable=True)  # Volume of distribution
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)

    compound = relationship("Compound", backref="pk_studies")
    created_by = relationship("User")


class ToxicologyResult(Base):
    """Toxicology assay result."""

    __tablename__ = "toxicology_results"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=False)
    assay_type = Column(String(50), nullable=False)  # herg, ames, cytotoxicity
    value = Column(Float, nullable=True)
    unit = Column(String(20), nullable=True)
    result = Column(String(50), nullable=True)  # positive, negative, inconclusive
    conditions = Column(JSON, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)

    compound = relationship("Compound", backref="toxicology_results")
    created_by = relationship("User")


class SavedQuestion(Base):
    """Saved user question for tracking Q&A."""

    __tablename__ = "saved_questions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    question_text = Column(Text, nullable=False)
    tags = Column(JSON, nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    is_archived = Column(Boolean, default=False, nullable=False)

    answers = relationship("SavedAnswer", back_populates="question", cascade="all, delete-orphan")
    created_by = relationship("User")


class SavedAnswer(Base):
    """Saved answer linked to a question with evidence and model metadata."""

    __tablename__ = "saved_answers"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    question_id = Column(UUID(as_uuid=True), ForeignKey("saved_questions.id"), nullable=False)
    answer_text = Column(Text, nullable=False)
    evidence = Column(JSON, nullable=True)
    model_used = Column(String(200), nullable=True)
    version = Column(Integer, default=1, nullable=False)
    confidence_score = Column(Float, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    question = relationship("SavedQuestion", back_populates="answers")


class HTSCampaign(Base):
    """High-throughput screening campaign."""

    __tablename__ = "hts_campaigns"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    campaign_id = Column(String(200), nullable=False, unique=True, index=True)
    campaign_name = Column(String(500), nullable=False)
    description = Column(Text, nullable=True)
    assay_type = Column(String(200), nullable=True)
    target = Column(String(500), nullable=True)
    library_id = Column(String(200), nullable=True)
    total_wells = Column(Integer, nullable=True)
    hit_count = Column(Integer, nullable=True, default=0)
    run_date = Column(DateTime, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # Relationships
    library_compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=True)
    library = relationship("Compound", back_populates="hts_campaigns")
    results = relationship("HTSResult", back_populates="campaign", cascade="all, delete-orphan")
    programs = relationship("Program", secondary="hts_campaign_program", back_populates="hts_campaigns")


class HTSResult(Base):
    """High-throughput screening result."""

    __tablename__ = "hts_results"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    result_id = Column(String(200), nullable=False, unique=True, index=True)
    campaign_id = Column(UUID(as_uuid=True), ForeignKey("hts_campaigns.id"), nullable=False, index=True)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=False, index=True)
    well_position = Column(String(50), nullable=True)
    raw_value = Column(Float, nullable=True)
    normalized_value = Column(Float, nullable=True)
    z_score = Column(Float, nullable=True)
    hit_flag = Column(Boolean, nullable=True, default=False)
    hit_category = Column(String(100), nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    campaign = relationship("HTSCampaign", back_populates="results")
    compound = relationship("Compound", back_populates="hts_results")


class BiochemicalResult(Base):
    """Biochemical assay result."""

    __tablename__ = "biochemical_results"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    result_id = Column(String(200), nullable=False, unique=True, index=True)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=False, index=True)
    assay_name = Column(String(500), nullable=False)
    target = Column(String(500), nullable=True)
    ic50 = Column(Float, nullable=True)
    ec50 = Column(Float, nullable=True)
    ki = Column(Float, nullable=True)
    kd = Column(Float, nullable=True)
    activity_type = Column(String(100), nullable=True)
    units = Column(String(50), nullable=True)
    run_date = Column(DateTime, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # Relationships
    compound = relationship("Compound", back_populates="biochemical_results")
    programs = relationship("Program", secondary="biochemical_result_program", back_populates="biochemical_results")


# Association tables for many-to-many relationships
compound_program = Table(
    "compound_program",
    Base.metadata,
    Column("compound_id", UUID(as_uuid=True), ForeignKey("compounds.id"), primary_key=True),
    Column("program_id", UUID(as_uuid=True), ForeignKey("programs.id"), primary_key=True),
    Column("status", String(100), nullable=True),
    Column("notes", Text, nullable=True),
    Column("created_at", DateTime, default=datetime.utcnow, nullable=False),
)

hts_campaign_program = Table(
    "hts_campaign_program",
    Base.metadata,
    Column("campaign_id", UUID(as_uuid=True), ForeignKey("hts_campaigns.id"), primary_key=True),
    Column("program_id", UUID(as_uuid=True), ForeignKey("programs.id"), primary_key=True),
)

biochemical_result_program = Table(
    "biochemical_result_program",
    Base.metadata,
    Column("result_id", UUID(as_uuid=True), ForeignKey("biochemical_results.id"), primary_key=True),
    Column("program_id", UUID(as_uuid=True), ForeignKey("programs.id"), primary_key=True),
)


class LabNotebookEntry(Base):
    """Electronic Lab Notebook entry for tracking experiments, observations, and notes."""

    __tablename__ = "lab_notebook_entries"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    title = Column(String(500), nullable=False)
    content = Column(Text, nullable=False)
    entry_type = Column(String(50), nullable=True)  # e.g., "experiment", "observation", "protocol", "note"
    tags = Column(ARRAY(String), nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False, index=True)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Optional attachments/metadata
    attachments = Column(JSON, nullable=True)  # Store file paths or references
    entry_metadata = Column(
        JSON, nullable=True
    )  # Additional structured data (renamed from 'metadata' - reserved in SQLAlchemy)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # Relationships (many-to-many with various entities via association table)
    linked_entities = relationship(
        "LabNotebookEntryAssociation",
        back_populates="entry",
        cascade="all, delete-orphan",
    )


class LabNotebookEntryAssociation(Base):
    """Association model for linking lab notebook entries to various entities."""

    __tablename__ = "lab_notebook_entry_assoc"

    entry_id = Column(UUID(as_uuid=True), ForeignKey("lab_notebook_entries.id"), primary_key=True)
    entity_type = Column(String(50), primary_key=True)
    entity_id = Column(UUID(as_uuid=True), primary_key=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    entry = relationship("LabNotebookEntry", back_populates="linked_entities")


class User(Base):
    """User account for authentication."""

    __tablename__ = "users"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    username = Column(String(100), unique=True, nullable=False, index=True)
    email = Column(String(255), unique=True, nullable=False)
    password_hash = Column(String(255), nullable=False)
    role = Column(String(50), default="researcher")  # admin, researcher, viewer
    is_active = Column(Boolean, default=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    last_login = Column(DateTime, nullable=True)


class Team(Base):
    """Team/organization for collaboration."""

    __tablename__ = "teams"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False, unique=True)
    description = Column(Text, nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    members = relationship("TeamMember", back_populates="team", cascade="all, delete-orphan")
    projects = relationship("Project", back_populates="team", cascade="all, delete-orphan")


class TeamMember(Base):
    """Team membership with role."""

    __tablename__ = "team_members"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    team_id = Column(UUID(as_uuid=True), ForeignKey("teams.id"), nullable=False)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    role = Column(String(20), default="member")  # owner, admin, member, viewer
    joined_at = Column(DateTime(timezone=True), server_default=func.now())

    team = relationship("Team", back_populates="members")
    user = relationship("User", backref="team_memberships")


class Project(Base):
    """Project within a team."""

    __tablename__ = "projects"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    team_id = Column(UUID(as_uuid=True), ForeignKey("teams.id"), nullable=False)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    is_public = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    team = relationship("Team", back_populates="projects")
    experiments = relationship("Experiment", backref="project")


class Feedback(Base):
    """User feedback and feature requests."""

    __tablename__ = "feedback"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    feedback_type = Column(String(20), nullable=False)  # bug, feature
    title = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    status = Column(String(20), default="new")  # new, reviewing, planned, done, wontfix
    priority = Column(String(20), default="medium")  # low, medium, high, critical
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())

    user = relationship("User", backref="feedback")


class UserFavorite(Base):
    """User favorite pages/bookmarks."""

    __tablename__ = "user_favorites"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    page_name = Column(String(100), nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", backref="favorites")


class Notification(Base):
    """User notifications."""

    __tablename__ = "notifications"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    title = Column(String(255), nullable=False)
    message = Column(Text, nullable=True)
    notification_type = Column(String(50), default="info")  # info, success, warning, discovery
    is_read = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", backref="notifications")


class EmailSubscription(Base):
    """Email subscription preferences for users."""

    __tablename__ = "email_subscriptions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False, index=True)
    subscription_type = Column(String(50), nullable=False)  # digest, alerts, shares
    frequency = Column(String(50), nullable=False)  # daily, weekly, immediate
    is_active = Column(Boolean, default=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", backref="email_subscriptions")


class Bookmark(Base):
    """User bookmarks for entities."""

    __tablename__ = "bookmarks"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    entity_type = Column(String(50), nullable=False)  # experiment, compound, signature
    entity_id = Column(UUID(as_uuid=True), nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", backref="bookmarks")


class ExperimentTemplate(Base):
    """Reusable experiment templates for quick creation."""

    __tablename__ = "experiment_templates"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    design_type = Column(String(50), nullable=True)
    organism = Column(ARRAY(String), nullable=True)
    sample_groups = Column(JSON, nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by = relationship("User", foreign_keys=[created_by_id])


class Note(Base):
    """User notes attached to entities."""

    __tablename__ = "notes"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    entity_type = Column(String(50), nullable=False)  # experiment, compound, signature, dataset, etc.
    entity_id = Column(UUID(as_uuid=True), nullable=False)
    content = Column(Text, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by = relationship("User", foreign_keys=[created_by_id])


class SavedFilter(Base):
    """Saved filter presets for entity lists."""

    __tablename__ = "saved_filters"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    entity_type = Column(String(50), nullable=False)  # experiment, compound, etc.
    filters = Column(JSON, nullable=False)  # Filter criteria as dict
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", foreign_keys=[user_id])


class WorkflowRule(Base):
    """Workflow automation rules."""

    __tablename__ = "workflow_rules"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    trigger_type = Column(String(100), nullable=False)  # e.g., "experiment_created", "compound_registered"
    trigger_config = Column(JSON, nullable=True)  # Filter conditions
    action_type = Column(String(100), nullable=False)  # e.g., "send_notification", "add_tag"
    action_config = Column(JSON, nullable=True)  # Action parameters
    is_active = Column(Boolean, default=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by = relationship("User", foreign_keys=[created_by_id])
    executions = relationship("WorkflowExecution", back_populates="rule", cascade="all, delete-orphan")


class WorkflowExecution(Base):
    """Workflow rule execution history."""

    __tablename__ = "workflow_executions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    rule_id = Column(UUID(as_uuid=True), ForeignKey("workflow_rules.id"), nullable=False)
    trigger_context = Column(JSON, nullable=True)  # What triggered it
    status = Column(String(50), nullable=False, default="pending")  # pending, success, failed
    result = Column(JSON, nullable=True)  # Execution result/error
    triggered_at = Column(DateTime(timezone=True), server_default=func.now())
    completed_at = Column(DateTime(timezone=True), nullable=True)

    rule = relationship("WorkflowRule", back_populates="executions")


class LiteratureCritique(Base):
    """Critical analysis of scientific literature."""

    __tablename__ = "literature_critiques"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    source_type = Column(String(50), nullable=False)  # text, rag_chunk, literature
    source_id = Column(UUID(as_uuid=True), nullable=True)  # Reference to source if applicable
    source_text = Column(Text, nullable=False)  # Original text analyzed
    critique_type = Column(String(50), nullable=False)  # critique, questions, contradiction
    content = Column(JSON, nullable=False)  # The analysis result
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by = relationship("User", foreign_keys=[created_by_id])


class TargetProductProfile(Base):
    """Target Product Profile (TPP) for candidate selection."""

    __tablename__ = "target_product_profiles"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    criteria = Column(JSON, nullable=False)  # Array of {property, min, max, weight, unit}
    is_active = Column(Boolean, default=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by = relationship("User", foreign_keys=[created_by_id])
    nominations = relationship("CandidateNomination", back_populates="tpp", cascade="all, delete-orphan")


class CandidateNomination(Base):
    """Candidate compound nomination against a TPP."""

    __tablename__ = "candidate_nominations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=False)
    tpp_id = Column(UUID(as_uuid=True), ForeignKey("target_product_profiles.id"), nullable=False)
    scores = Column(JSON, nullable=True)  # Per-criterion results {property: {value, score, traffic_light}}
    overall_score = Column(Float, nullable=True)
    status = Column(String(50), default="nominated")  # nominated, approved, rejected
    notes = Column(Text, nullable=True)
    nominated_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    compound = relationship("Compound", foreign_keys=[compound_id])
    tpp = relationship("TargetProductProfile", back_populates="nominations")
    nominated_by = relationship("User", foreign_keys=[nominated_by_id])


class AuditLog(Base):
    """Audit trail for user actions."""

    __tablename__ = "audit_logs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    username = Column(String(100), nullable=True)  # Cached for when user deleted
    action = Column(String(100), nullable=False)  # login, logout, create, update, delete
    entity_type = Column(String(100), nullable=True)  # experiment, dataset, signature, etc.
    entity_id = Column(String(100), nullable=True)  # UUID or identifier
    details = Column(JSON, nullable=True)  # Additional context
    ip_address = Column(String(50), nullable=True)
    timestamp = Column(DateTime, default=datetime.utcnow, index=True)

    user = relationship("User", backref="audit_logs")


class Protocol(Base):
    """Reusable experimental protocol template."""

    __tablename__ = "protocols"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    version = Column(Integer, default=1)
    category = Column(String(100), nullable=True)  # extraction, assay, analysis, etc.
    description = Column(Text, nullable=True)
    steps = Column(JSON, nullable=True)  # [{"order": 1, "title": "...", "instructions": "...", "duration": "..."}]
    materials = Column(JSON, nullable=True)  # [{"name": "...", "quantity": "...", "unit": "..."}]
    parameters = Column(JSON, nullable=True)  # {"temperature": "37C", "duration": "2h", ...}
    is_template = Column(Boolean, default=True)
    is_active = Column(Boolean, default=True)
    parent_id = Column(UUID(as_uuid=True), ForeignKey("protocols.id"), nullable=True)  # For versioning
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    parent = relationship("Protocol", remote_side=[id], backref="versions")
    created_by = relationship("User", foreign_keys=[created_by_id])


class ExperimentProtocol(Base):
    """Link between experiments and protocols used."""

    __tablename__ = "experiment_protocols"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=False)
    protocol_id = Column(UUID(as_uuid=True), ForeignKey("protocols.id"), nullable=False)
    executed_at = Column(DateTime, nullable=True)
    notes = Column(Text, nullable=True)
    deviations = Column(JSON, nullable=True)  # Any changes from standard protocol

    experiment = relationship("Experiment", backref="protocol_links")
    protocol = relationship("Protocol", backref="experiment_links")


class DiscoveryJob(Base):
    """Track automated repository discovery jobs."""

    __tablename__ = "discovery_jobs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    repository = Column(String(50), nullable=False)  # GEO, MW, PRIDE, ArrayExpress
    query = Column(String(500), nullable=False)  # Search query used
    status = Column(String(50), default="pending")  # pending, running, completed, failed
    studies_found = Column(Integer, default=0)
    studies_imported = Column(Integer, default=0)
    error_message = Column(Text, nullable=True)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)

    created_by = relationship("User", foreign_keys=[created_by_id])


class DiscoveredStudy(Base):
    """Studies found during discovery jobs."""

    __tablename__ = "discovered_studies"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    job_id = Column(UUID(as_uuid=True), ForeignKey("discovery_jobs.id"), nullable=False)
    study_id = Column(String(100), nullable=False)  # GSE123456, MTBLS123, etc.
    repository = Column(String(50), nullable=False)
    title = Column(String(500), nullable=True)
    description = Column(Text, nullable=True)
    organism = Column(String(200), nullable=True)
    omics_type = Column(String(100), nullable=True)
    status = Column(String(50), default="new")  # new, reviewed, imported, skipped
    imported_experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True)
    discovered_at = Column(DateTime, default=datetime.utcnow)

    job = relationship("DiscoveryJob", backref="discovered_studies")
    imported_experiment = relationship("Experiment")


class HarvestSchedule(Base):
    """Scheduled automated harvesting from repositories."""

    __tablename__ = "harvest_schedules"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    repository = Column(String(50), nullable=False)  # geo, metabolights
    query = Column(String(500), nullable=False)
    interval_hours = Column(Integer, default=24)
    is_active = Column(Boolean, default=True)
    last_run = Column(DateTime(timezone=True), nullable=True)
    next_run = Column(DateTime(timezone=True), nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by = relationship("User", foreign_keys=[created_by_id])


class StorageLocation(Base):
    """Physical storage hierarchy for samples (freezer/shelf/box/rack)."""

    __tablename__ = "storage_locations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    location_type = Column(String(50), nullable=True)  # freezer, shelf, box, rack
    parent_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    temperature = Column(String(50), nullable=True)
    capacity = Column(Integer, nullable=True)
    description = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)

    parent = relationship("StorageLocation", remote_side=[id], backref="children")


class Sample(Base):
    """Biological sample with storage and lineage tracking."""

    __tablename__ = "samples"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    sample_type = Column(String(100), nullable=True)
    barcode = Column(String(255), unique=True, nullable=True, index=True)
    storage_location_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    position = Column(String(100), nullable=True)  # e.g., box position
    parent_sample_id = Column(UUID(as_uuid=True), ForeignKey("samples.id"), nullable=True)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True)
    quantity = Column(Float, nullable=True)
    unit = Column(String(50), nullable=True)
    status = Column(String(50), default="available")  # available, depleted, reserved
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    notes = Column(Text, nullable=True)

    storage_location = relationship("StorageLocation", backref="samples")
    parent_sample = relationship("Sample", remote_side=[id], backref="child_samples")
    experiment = relationship("Experiment", backref="samples")
    created_by = relationship("User", foreign_keys=[created_by_id])


class SampleTransfer(Base):
    """Track movement of samples between storage locations."""

    __tablename__ = "sample_transfers"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    sample_id = Column(UUID(as_uuid=True), ForeignKey("samples.id"), nullable=False)
    from_location_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    to_location_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    transferred_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    transferred_at = Column(DateTime, default=datetime.utcnow)
    notes = Column(Text, nullable=True)

    sample = relationship("Sample", backref="transfers")
    from_location = relationship("StorageLocation", foreign_keys=[from_location_id], backref="outgoing_transfers")
    to_location = relationship("StorageLocation", foreign_keys=[to_location_id], backref="incoming_transfers")
    transferred_by = relationship("User", foreign_keys=[transferred_by_id])


# Add relationship to LabNotebookEntry after LabNotebookEntryAssociation is defined
LabNotebookEntry.linked_entities = relationship(
    "LabNotebookEntryAssociation",
    back_populates="entry",
    cascade="all, delete-orphan",
)


# Add relationships to Program model after all models are defined
# This avoids forward reference issues
Program.compounds = relationship("Compound", secondary=compound_program, back_populates="programs", viewonly=False)
Program.hts_campaigns = relationship(
    "HTSCampaign", secondary=hts_campaign_program, back_populates="programs", viewonly=False
)
Program.biochemical_results = relationship(
    "BiochemicalResult", secondary=biochemical_result_program, back_populates="programs", viewonly=False
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
