"""
SQLite schema definitions for chemistry and HTS data.

This module defines the database schema for storing:
- Compounds (SMILES, InChIKeys, internal IDs)
- Libraries (compound collections)
- HTS Campaigns (screening campaigns)
- HTS Results (assay results from screening)
- Biochemical Results (detailed biochemical assays)
- Compound-Program relationships
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

# SQLite schema definitions as CREATE TABLE statements

COMPOUNDS_TABLE = """
CREATE TABLE IF NOT EXISTS compounds (
    compound_id TEXT PRIMARY KEY,
    corporate_id TEXT UNIQUE,
    smiles TEXT NOT NULL,
    inchi_key TEXT UNIQUE,
    canonical_smiles TEXT,
    molecular_formula TEXT,
    molecular_weight REAL,
    logp REAL,
    hbd_count INTEGER,
    hba_count INTEGER,
    rotatable_bonds INTEGER,
    salt_form TEXT,
    batch_number TEXT,
    parent_compound_id TEXT REFERENCES compounds(compound_id),
    registration_date TIMESTAMP,
    registered_by TEXT,
    duplicate_of_id TEXT REFERENCES compounds(compound_id),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
"""

LIBRARIES_TABLE = """
CREATE TABLE IF NOT EXISTS libraries (
    library_id TEXT PRIMARY KEY,
    library_name TEXT NOT NULL,
    description TEXT,
    vendor TEXT,
    compound_count INTEGER DEFAULT 0,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
"""

LIBRARY_COMPOUNDS_TABLE = """
CREATE TABLE IF NOT EXISTS library_compounds (
    library_id TEXT NOT NULL,
    compound_id TEXT NOT NULL,
    vendor_id TEXT,
    vendor_name TEXT,
    PRIMARY KEY (library_id, compound_id),
    FOREIGN KEY (library_id) REFERENCES libraries(library_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id)
);
"""

HTS_CAMPAIGNS_TABLE = """
CREATE TABLE IF NOT EXISTS hts_campaigns (
    campaign_id TEXT PRIMARY KEY,
    campaign_name TEXT NOT NULL,
    description TEXT,
    assay_type TEXT,
    target TEXT,
    library_id TEXT,
    total_wells INTEGER,
    hit_count INTEGER DEFAULT 0,
    run_date DATE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (library_id) REFERENCES libraries(library_id)
);
"""

HTS_RESULTS_TABLE = """
CREATE TABLE IF NOT EXISTS hts_results (
    result_id TEXT PRIMARY KEY,
    campaign_id TEXT NOT NULL,
    compound_id TEXT NOT NULL,
    well_position TEXT,
    raw_value REAL,
    normalized_value REAL,
    z_score REAL,
    hit_flag INTEGER DEFAULT 0,
    hit_category TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (campaign_id) REFERENCES hts_campaigns(campaign_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id)
);
"""

BIOCHEMICAL_RESULTS_TABLE = """
CREATE TABLE IF NOT EXISTS biochemical_results (
    result_id TEXT PRIMARY KEY,
    compound_id TEXT NOT NULL,
    assay_name TEXT NOT NULL,
    target TEXT,
    ic50 REAL,
    ec50 REAL,
    ki REAL,
    kd REAL,
    activity_type TEXT,
    units TEXT,
    run_date DATE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id)
);
"""

COMPOUND_PROGRAM_TABLE = """
CREATE TABLE IF NOT EXISTS compound_program (
    compound_id TEXT NOT NULL,
    program_id TEXT NOT NULL,
    role TEXT,
    status TEXT,
    notes TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (compound_id, program_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id)
);
"""

COMPOUND_SIGNATURE_TABLE = """
CREATE TABLE IF NOT EXISTS compound_signature (
    compound_id TEXT NOT NULL,
    signature_id TEXT NOT NULL,
    effect_type TEXT,
    correlation REAL,
    p_value REAL,
    evidence_source TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (compound_id, signature_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id)
);
"""

COMPOUND_SEQUENCE_TABLE = """
CREATE TABLE IF NOT EXISTS compound_sequence (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    prefix TEXT DEFAULT 'AMP',
    next_number INTEGER DEFAULT 1
);
"""

# Indexes for performance
INDEXES = [
    "CREATE INDEX IF NOT EXISTS idx_compounds_inchi_key ON compounds(inchi_key);",
    "CREATE INDEX IF NOT EXISTS idx_compounds_smiles ON compounds(smiles);",
    "CREATE UNIQUE INDEX IF NOT EXISTS idx_compounds_corporate_id ON compounds(corporate_id);",
    "CREATE INDEX IF NOT EXISTS idx_hts_results_campaign ON hts_results(campaign_id);",
    "CREATE INDEX IF NOT EXISTS idx_hts_results_compound ON hts_results(compound_id);",
    "CREATE INDEX IF NOT EXISTS idx_hts_results_hit ON hts_results(hit_flag);",
    "CREATE INDEX IF NOT EXISTS idx_biochemical_compound ON biochemical_results(compound_id);",
    "CREATE INDEX IF NOT EXISTS idx_biochemical_target ON biochemical_results(target);",
    "CREATE INDEX IF NOT EXISTS idx_compound_program_program ON compound_program(program_id);",
    "CREATE INDEX IF NOT EXISTS idx_compound_signature_signature ON compound_signature(signature_id);",
]


@dataclass
class Compound:
    """Represents a chemical compound."""
    compound_id: str
    corporate_id: Optional[str] = None
    smiles: str
    inchi_key: Optional[str] = None
    canonical_smiles: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    hbd_count: Optional[int] = None
    hba_count: Optional[int] = None
    rotatable_bonds: Optional[int] = None
    salt_form: Optional[str] = None
    batch_number: Optional[str] = None
    parent_compound_id: Optional[str] = None
    registration_date: Optional[str] = None
    registered_by: Optional[str] = None
    duplicate_of_id: Optional[str] = None


@dataclass
class HTSCampaign:
    """Represents an HTS screening campaign."""
    campaign_id: str
    campaign_name: str
    description: Optional[str] = None
    assay_type: Optional[str] = None
    target: Optional[str] = None
    library_id: Optional[str] = None
    total_wells: Optional[int] = None
    hit_count: int = 0
    run_date: Optional[str] = None


@dataclass
class HTSResult:
    """Represents a single HTS assay result."""
    result_id: str
    campaign_id: str
    compound_id: str
    well_position: Optional[str] = None
    raw_value: Optional[float] = None
    normalized_value: Optional[float] = None
    z_score: Optional[float] = None
    hit_flag: int = 0
    hit_category: Optional[str] = None


@dataclass
class BiochemicalResult:
    """Represents a biochemical assay result."""
    result_id: str
    compound_id: str
    assay_name: str
    target: Optional[str] = None
    ic50: Optional[float] = None
    ec50: Optional[float] = None
    ki: Optional[float] = None
    kd: Optional[float] = None
    activity_type: Optional[str] = None
    units: Optional[str] = None
    run_date: Optional[str] = None


@dataclass
class CompoundSignatureLink:
    """
    Represents a link between a compound and a Postgres signature.

    This table is used to record compound–signature relationships discovered
    via screening, HTS, or analysis workflows.
    """

    compound_id: str
    signature_id: str
    effect_type: Optional[str] = None  # 'reverses', 'mimics', 'partial', 'unknown'
    correlation: Optional[float] = None  # -1.0 to 1.0
    p_value: Optional[float] = None
    evidence_source: Optional[str] = None
    created_at: Optional[str] = None
    updated_at: Optional[str] = None


