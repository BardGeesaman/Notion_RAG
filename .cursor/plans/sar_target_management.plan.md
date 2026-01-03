# SAR Audit + Target Management

## Overview
Two-phase plan: (1) Audit and polish existing SAR functionality with test coverage and documentation, (2) Build Target Management as a new first-class entity with druggability, assay linkage, and competitive landscape.

**Reviewer Status:** ✅ APPROVED WITH P1 FIXES (applied below)
**Timeline:** Phase 1: 3-5 days, Phase 2: 2-3 weeks

---

## Phase 1: SAR Audit & Polish (3-5 days)

### Current SAR State
Existing functionality in:
- `amprenta_rag/chemistry/sar_analysis.py` - Activity cliffs, Lipinski, SAR matrix
- `amprenta_rag/chemistry/rgroup.py` - MCS, R-group decomposition
- `amprenta_rag/api/routers/sar.py` - 9 API endpoints
- `scripts/dashboard/pages/chemistry/sar_analysis.py` - Dashboard

### Phase 1 Tasks

1. **Test Coverage Audit**
   - Run existing SAR tests, identify gaps
   - Add tests for untested edge cases
   - Target: 80%+ coverage for SAR modules

2. **API Documentation**
   - Add OpenAPI descriptions to all SAR endpoints
   - Document request/response schemas

3. **Dashboard Polish**
   - Verify all SAR dashboard tabs work end-to-end
   - Fix any broken visualizations
   - Add loading states and error handling

4. **Update ROADMAP**
   - Mark SAR Analysis Module as complete
   - Document any remaining gaps as P2/P3

---

## Phase 2: Target Management (2-3 weeks)

### Problem
Targets are currently stored as strings (`Experiment.targets`, `BiochemicalResult.target_name`). No first-class entity for:
- Druggability assessment
- Assay linkage
- Competitive landscape
- Validation tracking

### Data Model

```python
# Association table for compound-target relationships with activity data
compound_target = Table(
    'compound_target',
    Base.metadata,
    Column('compound_id', UUID, ForeignKey('compounds.id'), primary_key=True),
    Column('target_id', UUID, ForeignKey('targets.id'), primary_key=True),
    Column('activity_type', String(50)),  # IC50, Ki, Kd, EC50
    Column('activity_value', Float),
    Column('activity_units', String(20)),  # nM, uM
    Column('assay_id', UUID, ForeignKey('experiments.id'), nullable=True)
)

class Target(Base):
    __tablename__ = "targets"
    
    id = Column(UUID, primary_key=True)
    name = Column(String(200), nullable=False, unique=True)
    description = Column(Text)  # P1 FIX: Added
    gene_symbol = Column(String(50), index=True)
    uniprot_id = Column(String(20), index=True)
    
    # Classification
    target_class = Column(String(100))  # kinase, gpcr, ion_channel, etc.
    target_family = Column(String(100))  # e.g., "Tyrosine kinase"
    
    # Druggability
    druggability_score = Column(Float)  # 0-1
    druggability_source = Column(String(50))  # manual, cansar, opentargets
    pocket_count = Column(Integer)
    
    # Validation
    validation_status = Column(String(50))  # hypothesis, validated, clinical
    validation_evidence = Column(JSON)  # Links to papers, experiments
    
    # Lifecycle (P1 FIX: Added)
    lifecycle_status = Column(String(50), default="active")  # active, archived, deprecated
    is_synthetic = Column(Boolean, default=False)  # Mutant, fusion protein, etc.
    
    # External IDs
    chembl_id = Column(String(50))
    ensembl_id = Column(String(50))
    
    # Relationships
    assays = relationship("Experiment", back_populates="target")
    compounds = relationship("Compound", secondary=compound_target)
    
    created_at = Column(DateTime)
    updated_at = Column(DateTime)
```

### API Endpoints (8)

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | /targets | List targets with filters |
| GET | /targets/{id} | Get target details |
| POST | /targets | Create target |
| PUT | /targets/{id} | Update target |
| GET | /targets/{id}/assays | Get linked assays |
| GET | /targets/{id}/compounds | Get compounds with activity |
| GET | /targets/{id}/landscape | Competitive landscape |
| POST | /targets/{id}/druggability | Calculate druggability |

### Dashboard (4 tabs)

1. **Target List** - Browse, search, filter targets
2. **Target Detail** - Profile, druggability, validation status
3. **Assay Linkage** - Link experiments/assays to target
4. **Competitive Landscape** - Known drugs, clinical compounds, literature

### External Integrations

**Phase 2A (MVP):**
- **UniProt** - Gene info, protein sequence, function
- **ChEMBL** - Approved drugs, clinical candidates, bioactivity

**Phase 2B (Defer):**
- **Open Targets** - Complex GraphQL API, disease associations
- **CanSAR** - Druggability (subscription unclear)

---

## Implementation Order

```
Phase 1: SAR Audit (3-5 days)
├── Day 1: Test coverage audit + gap analysis
├── Day 2-3: Add missing tests
├── Day 4: Dashboard polish + docs
└── Day 5: ROADMAP update

Phase 2: Target Management (2-3 weeks)
├── Week 1:
│   ├── Day 1-2: Target model + migration
│   ├── Day 3-4: API endpoints (CRUD)
│   └── Day 5: Service layer + tests
├── Week 2:
│   ├── Day 1-2: Dashboard (4 tabs)
│   ├── Day 3-4: External integrations (UniProt, ChEMBL)
│   └── Day 5: Druggability calculation
└── Week 3:
    ├── Day 1-2: E2E tests (6)
    ├── Day 3: Documentation
    └── Day 4-5: Buffer for issues
```

## Success Criteria

**Phase 1:**
- SAR test coverage >= 80%
- All SAR dashboard tabs functional
- ROADMAP updated

**Phase 2:**
- Target model with migration
- 8 API endpoints with tests
- 4-tab dashboard
- UniProt/ChEMBL integration
- 30+ tests (12 service + 12 API + 6 E2E)

---

## TODOs

- [ ] Phase 1: SAR test coverage audit and gap analysis
- [ ] Phase 1: Add missing SAR tests (target 80%+ coverage)
- [ ] Phase 1: Dashboard polish and documentation
- [ ] Phase 2: Target model + Alembic migration
- [ ] Phase 2: Target API endpoints (8) + service layer
- [ ] Phase 2: Target Management dashboard (4 tabs)
- [ ] Phase 2: UniProt/ChEMBL integrations
- [ ] Phase 2: Tests (30+) and documentation

