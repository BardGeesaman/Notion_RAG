# Multi-Modality Therapeutics Platform Expansion (v2)

**Updated based on Reviewer feedback**

Strategic 2-3 year roadmap to extend Amprenta RAG from small molecules to all major therapeutic modalities with equivalent data management, analysis, and discovery capabilities.

---

## Current State (Small Molecules)

The platform currently supports:
- **Data Models**: Compound, SMILES, molecular descriptors, assay results
- **Analysis**: ADMET prediction, SAR analysis, similarity search, fingerprinting
- **Workflows**: HTS screening, hit-to-lead, lead optimization
- **Integrations**: ChEMBL, PubChem, RDKit

---

## Revised Phase Structure (Per Reviewer)

```
Phase 0 (Foundation) - Q1
    │
    ├── Phase 1A: Peptides        ← Bridge between chemistry and biology
    │       │
    │       ↓
    ├── Phase 1B: Antibodies      ← Core biologics (defer ADCs)
    │       │
    │       ↓
    ├── Phase 1C: Formulation     ← Cross-cutting for all biologics
    │       │
    │       ↓
    ├── Phase 2A: Nucleic Acids   ← Reuses sequence infrastructure
    │       │
    │       ↓
    ├── Phase 2B: Vaccines        ← Uses epitope mapper from biologics
    │       │
    │       ↓
    ├── Phase 2C: Cell Therapy    ← Manufacturing + patient tracking
    │       │
    │       ↓
    ├── Phase 3A: ADCs            ← Multi-modal (antibody + compound)
    │       │
    │       ↓
    └── Phase 3B: Gene Therapy    ← Most complex, last
```

---

## Phase 0: Extensible Foundation (Q1) - REVISED

### 0.1 Database Architecture: Joined-Table Inheritance (P1 FIX)

**Use SQLAlchemy joined-table inheritance** (NOT single-table):

```python
# amprenta_rag/models/therapeutics.py (NEW)

from sqlalchemy import Column, String, DateTime, ForeignKey
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import declared_attr
from amprenta_rag.database.base import Base

class TherapeuticEntity(Base):
    """Base class for all therapeutic modalities using joined-table inheritance."""
    __tablename__ = "therapeutic_entities"
    
    id = Column(UUID(as_uuid=True), primary_key=True)
    modality = Column(String(50), nullable=False, index=True)  # Discriminator
    name = Column(String(500), nullable=False)
    description = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    program_id = Column(UUID(as_uuid=True), ForeignKey("programs.id"), nullable=True)
    lifecycle_status = Column(String(20), default="active", index=True)
    
    __mapper_args__ = {
        "polymorphic_on": modality,
        "polymorphic_identity": "base",
    }


class SmallMolecule(TherapeuticEntity):
    """Small molecule compounds (migrated from Compound)."""
    __tablename__ = "small_molecules"
    
    id = Column(UUID(as_uuid=True), ForeignKey("therapeutic_entities.id"), primary_key=True)
    smiles = Column(String(2000), nullable=False)
    inchi_key = Column(String(27), nullable=True, index=True)
    mol_weight = Column(Float, nullable=True)
    # ... existing Compound fields
    
    __mapper_args__ = {"polymorphic_identity": "small_molecule"}
```

### 0.2 Migration Strategy for Existing Compound Data (P1 FIX)

**Alembic Migration Plan**:

```python
# alembic/versions/xxx_migrate_compound_to_therapeutic_entity.py

def upgrade():
    # Step 1: Create therapeutic_entities table
    op.create_table("therapeutic_entities", ...)
    
    # Step 2: Create small_molecules table (inherits from therapeutic_entities)
    op.create_table("small_molecules", ...)
    
    # Step 3: Migrate existing Compound data
    op.execute("""
        INSERT INTO therapeutic_entities (id, modality, name, created_at, created_by_id, program_id)
        SELECT id, 'small_molecule', name, created_at, created_by_id, program_id
        FROM compounds
    """)
    
    op.execute("""
        INSERT INTO small_molecules (id, smiles, inchi_key, mol_weight, ...)
        SELECT id, smiles, inchi_key, mol_weight, ...
        FROM compounds
    """)
    
    # Step 4: Update foreign keys in dependent tables
    # (AssayResult, ActivityEvent, etc.)
    
    # Step 5: Rename compounds -> compounds_legacy (keep for rollback)
    op.rename_table("compounds", "compounds_legacy")

def downgrade():
    # Reverse migration
    op.rename_table("compounds_legacy", "compounds")
    op.drop_table("small_molecules")
    op.drop_table("therapeutic_entities")
```

### 0.3 Plugin Architecture: ModalityAnalyzer Protocol (P1 FIX)

**Define interface contract before implementation**:

```python
# amprenta_rag/analysis/registry.py (NEW)

from typing import Protocol, Dict, Any, List
from amprenta_rag.models.therapeutics import TherapeuticEntity

class PropertyPrediction(BaseModel):
    """Standard prediction result."""
    property_name: str
    value: float
    unit: Optional[str]
    confidence: Optional[float]
    model_version: str

class ValidationError(BaseModel):
    """Validation error details."""
    field: str
    message: str
    severity: str  # error, warning, info

class SimilarEntity(BaseModel):
    """Similarity search result."""
    entity_id: UUID
    similarity_score: float
    modality: str


class ModalityAnalyzer(Protocol):
    """Protocol for modality-specific analyzers."""
    
    modality: str
    
    def get_descriptors(self, entity: TherapeuticEntity) -> Dict[str, Any]:
        """Calculate modality-specific descriptors."""
        ...
    
    def predict_properties(self, entity: TherapeuticEntity) -> List[PropertyPrediction]:
        """Predict properties relevant to this modality."""
        ...
    
    def validate(self, entity: TherapeuticEntity) -> List[ValidationError]:
        """Validate entity data for this modality."""
        ...
    
    def similarity_search(
        self, 
        entity: TherapeuticEntity, 
        db: Session, 
        limit: int = 10
    ) -> List[SimilarEntity]:
        """Find similar entities within this modality."""
        ...
    
    def get_external_references(self, entity: TherapeuticEntity) -> List[ExternalReference]:
        """Fetch external database references."""
        ...


class AnalyzerRegistry:
    """Registry for modality-specific analyzers."""
    
    _analyzers: Dict[str, ModalityAnalyzer] = {}
    
    @classmethod
    def register(cls, analyzer: ModalityAnalyzer) -> None:
        cls._analyzers[analyzer.modality] = analyzer
    
    @classmethod
    def get(cls, modality: str) -> Optional[ModalityAnalyzer]:
        return cls._analyzers.get(modality)
    
    @classmethod
    def analyze(cls, entity: TherapeuticEntity) -> Dict[str, Any]:
        analyzer = cls.get(entity.modality)
        if not analyzer:
            raise ValueError(f"No analyzer registered for {entity.modality}")
        return analyzer.get_descriptors(entity)
```

### 0.4 Shared Sequence Support

```python
# amprenta_rag/models/sequences.py (NEW)

class Sequence(BaseModel):
    """Unified sequence representation."""
    sequence: str
    sequence_type: str  # protein, dna, rna
    modifications: Optional[List[SequenceModification]] = None
    
    def validate_alphabet(self) -> bool:
        """Validate sequence uses correct alphabet."""
        ...
    
    def to_fasta(self) -> str:
        """Export as FASTA format."""
        ...


class SequenceModification(BaseModel):
    """Chemical modification at specific position."""
    position: int
    modification_type: str  # e.g., "2'-OMe", "phosphorothioate", "N-methyl"
    description: Optional[str]
```

---

## Phase 1A: Peptides Module (Q2) - ELEVATED

**Rationale**: Peptides bridge small molecules and proteins. They share chemistry concepts (modifications, cyclization) while introducing sequence-based analysis.

### Data Models

```python
class Peptide(TherapeuticEntity):
    __tablename__ = "peptides"
    
    id = Column(UUID, ForeignKey("therapeutic_entities.id"), primary_key=True)
    sequence = Column(Text, nullable=False)  # Standard + non-natural AAs
    sequence_type = Column(String(20), default="linear")  # linear, cyclic, stapled
    modifications = Column(JSON)  # Non-natural AAs, D-forms, N-methylation
    cyclization_type = Column(String(50))  # head-to-tail, disulfide, staple
    molecular_weight = Column(Float)
    
    __mapper_args__ = {"polymorphic_identity": "peptide"}
```

### Analysis Capabilities
- Stability prediction (proteolytic cleavage sites)
- Cell permeability prediction (Rule of 5 for peptides)
- Helicity calculation (for stapled peptides)
- Non-natural amino acid handling (P2 gap fix)
- PK modeling integration (P2 gap fix)

### External Integrations
- PeptideAtlas
- APD (Antimicrobial Peptide Database)
- CycPeptMPDB (cyclic peptide membrane permeability)

---

## Phase 1B: Antibodies Module (Q2-Q3)

**Split from "Biologics" per Reviewer - defer ADCs to Phase 3A**

### Data Models

```python
class Antibody(TherapeuticEntity):
    __tablename__ = "antibodies"
    
    id = Column(UUID, ForeignKey("therapeutic_entities.id"), primary_key=True)
    heavy_chain = Column(Text, nullable=False)
    light_chain = Column(Text, nullable=False)
    isotype = Column(String(20))  # IgG1, IgG2, IgG4, etc.
    format = Column(String(50))  # full-length, Fab, scFv, bispecific
    cdr_regions = Column(JSON)  # CDR1-3 for H and L chains
    framework_regions = Column(JSON)
    species_origin = Column(String(50))  # human, humanized, chimeric
    target_antigen = Column(String(200))
    
    __mapper_args__ = {"polymorphic_identity": "antibody"}


class FusionProtein(TherapeuticEntity):
    __tablename__ = "fusion_proteins"
    
    id = Column(UUID, ForeignKey("therapeutic_entities.id"), primary_key=True)
    domains = Column(JSON)  # List of domain definitions
    linker_sequences = Column(JSON)
    fc_region = Column(Text)
    
    __mapper_args__ = {"polymorphic_identity": "fusion_protein"}
```

### Analysis Capabilities
- CDR identification (IMGT, Kabat, Chothia numbering)
- Developability assessment
  - Aggregation propensity
  - Viscosity prediction
  - Chemical liability hotspots
- Immunogenicity prediction (T-cell epitopes)
- PTM prediction (glycosylation, deamidation)
- Manufacturability scoring (P2 gap fix)
- Biosimilar comparison tools (P2 gap fix)

### External Integrations
- UniProt
- IMGT (immunogenetics)
- SAbDab (antibody structures)
- PDB
- **TherapeuticTarget Database** (P2 gap fix)
- **BindingDB** (P2 gap fix)

---

## Phase 1C: Formulation Module (Q3) - ELEVATED

**Rationale**: Critical cross-cutting capability for all biologics per Reviewer.

### Data Models

```python
class Formulation(Base):
    __tablename__ = "formulations"
    
    id = Column(UUID, primary_key=True)
    therapeutic_entity_id = Column(UUID, ForeignKey("therapeutic_entities.id"))
    name = Column(String(200))
    form = Column(String(50))  # liquid, lyophilized, spray-dried
    excipients = Column(JSON)  # Buffer, surfactant, cryoprotectant
    ph = Column(Float)
    concentration = Column(Float)
    storage_conditions = Column(JSON)  # Temperature, humidity
    stability_data = Column(JSON)  # Accelerated, real-time
    
    # Cold chain tracking (P2 gap fix)
    cold_chain_required = Column(Boolean, default=True)
    min_temp = Column(Float)
    max_temp = Column(Float)
    excursion_tolerance = Column(JSON)
```

### Capabilities
- Stability profiling (thermal, oxidative, photolytic)
- Aggregation prediction
- Viscosity modeling
- Cold chain requirement analysis (P2 gap fix)
- Container-closure compatibility

---

## Phase 2A: Nucleic Acid Therapeutics (Q4 - Year 2 Q1)

### Data Models

```python
class NucleicAcidTherapeutic(TherapeuticEntity):
    __tablename__ = "nucleic_acids"
    
    id = Column(UUID, ForeignKey("therapeutic_entities.id"), primary_key=True)
    sequence = Column(Text, nullable=False)
    oligo_type = Column(String(20))  # ASO, siRNA, mRNA, aptamer
    target_gene = Column(String(100))
    target_sequence = Column(Text)
    modifications = Column(JSON)  # Backbone, sugar, base modifications
    delivery_vehicle = Column(String(50))  # LNP, GalNAc, exosome (P2 gap fix)
    
    __mapper_args__ = {"polymorphic_identity": "nucleic_acid"}


class DeliveryVehicle(Base):
    """LNP, GalNAc, exosomes, etc. (P2 gap fix)"""
    __tablename__ = "delivery_vehicles"
    
    id = Column(UUID, primary_key=True)
    name = Column(String(200))
    vehicle_type = Column(String(50))  # LNP, conjugate, nanoparticle
    components = Column(JSON)  # Lipid composition, targeting ligand
    target_tissue = Column(String(100))
    encapsulation_efficiency = Column(Float)
```

### Analysis Capabilities
- Secondary structure prediction (RNAfold)
- Off-target prediction (seed matching, BLAST)
- Modification impact analysis
- Tissue distribution modeling
- Synthesis cost estimation (P2 gap fix)

### External Integrations
- Ensembl
- miRBase
- RNAcentral
- **NCBI Gene** (P2 gap fix)

---

## Phase 2B: Vaccines Module (Year 2 Q2)

### Data Models

```python
class VaccineAntigen(TherapeuticEntity):
    __tablename__ = "vaccine_antigens"
    
    id = Column(UUID, ForeignKey("therapeutic_entities.id"), primary_key=True)
    antigen_type = Column(String(50))  # protein, peptide, mRNA, viral_vector
    sequence = Column(Text)
    epitopes = Column(JSON)  # B-cell, T-cell epitopes
    expression_system = Column(String(100))
    pathogen_target = Column(String(200))
    
    __mapper_args__ = {"polymorphic_identity": "vaccine_antigen"}


class VaccineFormulation(Base):
    __tablename__ = "vaccine_formulations"
    
    id = Column(UUID, primary_key=True)
    antigen_id = Column(UUID, ForeignKey("vaccine_antigens.id"))
    adjuvant_id = Column(UUID, ForeignKey("adjuvants.id"))
    formulation_type = Column(String(50))
    dose = Column(Float)
    route = Column(String(50))  # IM, SC, intranasal


class Adjuvant(Base):
    """Reference table of adjuvants (P2 gap fix)"""
    __tablename__ = "adjuvants"
    
    id = Column(UUID, primary_key=True)
    name = Column(String(200))
    adjuvant_type = Column(String(100))  # alum, oil-in-water, TLR agonist
    mechanism = Column(Text)
    approved_indications = Column(JSON)
```

### Analysis Capabilities
- Epitope prediction (B-cell, MHC-I, MHC-II)
- Antigen design optimization
- Adjuvant selection guidance
- Cold chain modeling (P2 gap fix)

### External Integrations
- IEDB (Immune Epitope Database)
- **Alternative to NetMHC** - use open-source MHCflurry (P1 licensing fix)

---

## Phase 2C: Cell Therapy Module (Year 2 Q3-Q4)

### Data Models

```python
class CellProduct(TherapeuticEntity):
    __tablename__ = "cell_products"
    
    id = Column(UUID, ForeignKey("therapeutic_entities.id"), primary_key=True)
    cell_type = Column(String(100))  # CAR-T, CAR-NK, TIL, MSC
    source = Column(String(50))  # autologous, allogeneic
    phenotype = Column(JSON)  # Surface markers, cytokine profile
    
    __mapper_args__ = {"polymorphic_identity": "cell_product"}


class CARConstruct(Base):
    __tablename__ = "car_constructs"
    
    id = Column(UUID, primary_key=True)
    cell_product_id = Column(UUID, ForeignKey("cell_products.id"))
    scfv_sequence = Column(Text)
    hinge_domain = Column(String(50))
    transmembrane_domain = Column(String(50))
    costimulatory_domains = Column(JSON)  # 4-1BB, CD28
    signaling_domain = Column(String(50))  # CD3zeta


class ManufacturingBatch(Base):
    """Track manufacturing (P2 gap fix)"""
    __tablename__ = "manufacturing_batches"
    
    id = Column(UUID, primary_key=True)
    cell_product_id = Column(UUID, ForeignKey("cell_products.id"))
    batch_number = Column(String(100))
    patient_id = Column(UUID, ForeignKey("patients.id"), nullable=True)
    apheresis_date = Column(DateTime)
    manufacturing_start = Column(DateTime)
    manufacturing_end = Column(DateTime)
    release_date = Column(DateTime)
    qc_metrics = Column(JSON)  # Viability, transduction efficiency, potency
    status = Column(String(50))


class PotencyAssay(Base):
    """Standardized QC (P2 gap fix)"""
    __tablename__ = "potency_assays"
    
    id = Column(UUID, primary_key=True)
    batch_id = Column(UUID, ForeignKey("manufacturing_batches.id"))
    assay_type = Column(String(100))
    result = Column(Float)
    specification = Column(Float)
    pass_fail = Column(Boolean)
```

### Analysis Capabilities
- CAR construct optimization
- Expansion/persistence modeling
- CRS risk prediction
- Patient logistics tracking (P2 gap fix)

### External Integrations
- FlowRepository
- ClinicalTrials.gov
- **ImmPort** (P2 gap fix)

---

## Phase 3A: ADCs (Multi-Modal) (Year 3 Q1)

**Deferred and split per Reviewer - requires both antibody and compound infrastructure**

### Data Models

```python
class ADC(TherapeuticEntity):
    """Antibody-Drug Conjugate - multi-modal entity"""
    __tablename__ = "adcs"
    
    id = Column(UUID, ForeignKey("therapeutic_entities.id"), primary_key=True)
    
    # Cross-modality references (P2 multi-modal fix)
    antibody_id = Column(UUID, ForeignKey("antibodies.id"), nullable=False)
    payload_id = Column(UUID, ForeignKey("small_molecules.id"), nullable=False)
    
    linker_chemistry = Column(String(200))
    linker_type = Column(String(50))  # cleavable, non-cleavable
    conjugation_site = Column(String(100))  # lysine, cysteine, site-specific
    dar = Column(Float)  # Drug-to-antibody ratio
    dar_distribution = Column(JSON)  # Distribution of DAR species
    
    __mapper_args__ = {"polymorphic_identity": "adc"}
    
    # Relationships
    antibody = relationship("Antibody")
    payload = relationship("SmallMolecule")
```

### Analysis Capabilities
- DAR optimization
- Linker stability prediction
- Bystander effect modeling
- Integrated PK/PD (antibody + payload)

---

## Phase 3B: Gene Therapy Module (Year 3 Q2-Q4)

**Last phase per Reviewer - most complex, least public data**

### Data Models

```python
class GeneTherapyVector(TherapeuticEntity):
    __tablename__ = "gene_therapy_vectors"
    
    id = Column(UUID, ForeignKey("therapeutic_entities.id"), primary_key=True)
    vector_type = Column(String(50))  # AAV, lentiviral, adenoviral
    serotype = Column(String(20))  # AAV2, AAV9, etc.
    tropism = Column(JSON)  # Target tissues
    transgene_id = Column(UUID, ForeignKey("transgenes.id"))
    
    __mapper_args__ = {"polymorphic_identity": "gene_therapy_vector"}


class Transgene(Base):
    __tablename__ = "transgenes"
    
    id = Column(UUID, primary_key=True)
    gene_name = Column(String(100))
    sequence = Column(Text)
    promoter = Column(String(100))
    regulatory_elements = Column(JSON)
    expression_profile = Column(JSON)


class VectorProduction(Base):
    """Manufacturing QC (P2 gap fix)"""
    __tablename__ = "vector_productions"
    
    id = Column(UUID, primary_key=True)
    vector_id = Column(UUID, ForeignKey("gene_therapy_vectors.id"))
    batch_number = Column(String(100))
    titer = Column(Float)  # vg/mL
    purity = Column(Float)
    potency = Column(Float)
    empty_full_ratio = Column(Float)
```

### Analysis Capabilities
- Capsid engineering
- Promoter strength prediction
- Immunogenicity (anti-AAV antibodies)
- Expression durability modeling (P2 gap fix)
- Biodistribution prediction

### External Integrations
- AAV capsid databases
- GTEx (tissue expression)

---

## Cross-Cutting Features (All Phases)

### Translational Biomarkers (Phase 2)

Per Reviewer - connects entities to patient outcomes.

```python
class TranslationalBiomarker(Base):
    __tablename__ = "translational_biomarkers"
    
    id = Column(UUID, primary_key=True)
    therapeutic_entity_id = Column(UUID, ForeignKey("therapeutic_entities.id"))
    biomarker_name = Column(String(200))
    biomarker_type = Column(String(50))  # PD, efficacy, safety, resistance
    measurement_method = Column(String(200))
    clinical_correlation = Column(Text)
    patient_stratification = Column(JSON)
```

### Cost Modeling (Phase 2)

COGS estimation per modality.

```python
class CostModel(Base):
    __tablename__ = "cost_models"
    
    id = Column(UUID, primary_key=True)
    therapeutic_entity_id = Column(UUID, ForeignKey("therapeutic_entities.id"))
    manufacturing_cost = Column(Float)
    raw_materials_cost = Column(Float)
    qc_testing_cost = Column(Float)
    estimated_cogs = Column(Float)
    scale = Column(String(50))  # pilot, commercial
```

### Regulatory Pathway Mapping (Phase 2)

IND/BLA/NDA requirements.

```python
class RegulatoryPathway(Base):
    __tablename__ = "regulatory_pathways"
    
    id = Column(UUID, primary_key=True)
    therapeutic_entity_id = Column(UUID, ForeignKey("therapeutic_entities.id"))
    pathway_type = Column(String(50))  # IND, BLA, NDA, 510k
    designation = Column(String(100))  # breakthrough, fast_track, orphan
    required_studies = Column(JSON)
    timeline_estimate = Column(Integer)  # months
```

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| **Migration breaks existing workflows** | Phased migration with `compounds_legacy` table for rollback |
| **NetMHC licensing** | Use open-source MHCflurry instead |
| **Polymorphic query performance** | Add modality-specific indexes, consider materialized views |
| **Data scarcity for CGT** | Partner with academic institutions, use synthetic data for testing |
| **Domain expertise** | Hire/consult modality specialists for each phase |

---

## Success Metrics

| Metric | Target |
|--------|--------|
| Modalities supported | 7 (small molecule, peptide, antibody, nucleic acid, vaccine, cell, gene) |
| Multi-modal entities | ADCs properly reference both antibody and compound |
| Shared code reuse | >60% across modalities |
| API consistency | Unified CRUD via TherapeuticEntity base |
| Migration success | 100% compound data migrated without data loss |
| Test coverage | >80% per module |
| Plugin compliance | All analyzers implement ModalityAnalyzer protocol |

---

## Timeline Summary

| Phase | Content | Duration |
|-------|---------|----------|
| **Phase 0** | Foundation (inheritance, migration, plugins) | Q1 |
| **Phase 1A** | Peptides | Q2 |
| **Phase 1B** | Antibodies | Q2-Q3 |
| **Phase 1C** | Formulation | Q3 |
| **Phase 2A** | Nucleic Acids | Q4-Y2Q1 |
| **Phase 2B** | Vaccines | Y2Q2 |
| **Phase 2C** | Cell Therapy | Y2Q3-Q4 |
| **Phase 3A** | ADCs (multi-modal) | Y3Q1 |
| **Phase 3B** | Gene Therapy | Y3Q2-Q4 |

**Total**: ~3 years for full implementation

---

## Immediate Next Steps

1. Save this plan to Strategic Backlog in ROADMAP.md
2. Begin Phase 0 design when ready to implement
3. Identify domain experts for modality-specific guidance

