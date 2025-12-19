"""Chemistry domain models split from the main models module."""

from __future__ import annotations

from datetime import datetime

from sqlalchemy import Boolean, Column, DateTime, Float, ForeignKey, Integer, JSON, String, Table, Text, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import relationship

from amprenta_rag.database.base import Base
import uuid


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


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
    aromatic_rings = Column(Integer, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Migration support: Notion page ID
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    # External identifiers
    external_ids = Column(JSON, nullable=True)

    # Version for concurrent editing safety
    version = Column(Integer, default=1, nullable=False)

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


class GenericAssayResult(Base):
    """Generic assay result that can be linked to experiments or compounds."""

    __tablename__ = "generic_assay_results"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    assay_name = Column(String(200), nullable=False, index=True)
    assay_type = Column(String(100), nullable=False, index=True)  # e.g., "biochemical", "cellular", "in_vivo", "custom"
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True, index=True)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=True, index=True)
    result_data = Column(JSON, nullable=False)  # Flexible JSON structure for any assay data
    assay_metadata = Column(JSON, nullable=True)  # Additional metadata (conditions, units, etc.)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    experiment = relationship("Experiment", backref="generic_assay_results")
    compound = relationship("Compound", backref="generic_assay_results")
    created_by = relationship("User", foreign_keys=[created_by_id])


__all__ = [
    "Compound",
    "BiochemicalAssay",
    "ActivityResult",
    "ADMEResult",
    "PKStudy",
    "ToxicologyResult",
    "HTSCampaign",
    "HTSResult",
    "BiochemicalResult",
    "TargetProductProfile",
    "CandidateNomination",
    "GenericAssayResult",
    "compound_program",
    "hts_campaign_program",
    "biochemical_result_program",
]

