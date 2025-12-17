# Test Data Seeding Guide

**Last Updated**: December 2025

---

## 1. Overview

The test data seeding system provides **idempotent, deterministic scripts** for populating the Amprenta RAG database with realistic sample data across all supported domains:

- **Transcriptomics** - Gene expression datasets, RNA-seq experiments
- **Proteomics** - Protein abundance, mass spec data
- **Metabolomics** - Small molecule profiling, metabolite features
- **Lipidomics** - Lipid species identification and quantification
- **HTS/SAR** - High-throughput screening campaigns, structure-activity relationships
- **Signatures** - Multi-omics signatures, feature lists, biological programs

### Purpose

- **Development**: Quickly populate a local database for feature development
- **Testing**: Generate consistent test datasets for integration and E2E tests
- **Demos**: Create realistic demo environments for stakeholders
- **CI/CD**: Seed ephemeral test environments in automated pipelines

### Key Features

- **Idempotent**: Safe to run multiple times; uses upsert patterns
- **Deterministic**: Seeded random number generators ensure reproducible data
- **Size-aware**: Small/medium/large presets for different use cases
- **Domain-specific**: Each seeder understands its domain's realistic distributions
- **Fast**: Bulk inserts optimized for performance

---

## 2. Quick Start

### Prerequisites

- PostgreSQL database running and configured in `.env`
- Python environment with `requirements.txt` installed
- Database schema migrations applied (`alembic upgrade head`)

### Seed Everything (Default)

```bash
python scripts/seed_all.py
```

This runs all seeders with **small** size preset.

### Seed with Custom Size

```bash
# Medium dataset
python scripts/seed_all.py --size medium

# Large dataset (warning: may take 5-10 minutes)
python scripts/seed_all.py --size large
```

### Reset and Reseed

```bash
# Clear all data and reseed
python scripts/seed_all.py --reset --size medium
```

### Dry Run (Preview)

```bash
# See what would be created without making changes
python scripts/seed_all.py --dry-run --size large
```

---

## 3. CLI Reference

### seed_all.py (Orchestrator)

```bash
python scripts/seed_all.py [OPTIONS]
```

#### Flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--size` | choice | `small` | Size preset: `small`, `medium`, or `large` |
| `--reset` | flag | `False` | Truncate all tables before seeding (⚠️ destructive) |
| `--seed` | int | `1234` | Random seed for deterministic generation |
| `--dry-run` | flag | `False` | Preview what would be created without executing |
| `--only` | text | `None` | Run only specific seeder (e.g., `transcriptomics` or `hts`) |

#### Examples

```bash
# Seed only transcriptomics
python scripts/seed_all.py --only transcriptomics

# Large dataset with custom random seed
python scripts/seed_all.py --size large --seed 12345

# Reset and seed medium dataset
python scripts/seed_all.py --reset --size medium

# Dry run to preview large dataset
python scripts/seed_all.py --dry-run --size large
```

### Individual Seeders

Each domain has its own seeder script with the same CLI interface:

```bash
python scripts/seed_transcriptomics_data.py --size medium --seed 1234
python scripts/seed_proteomics_data.py --size medium --seed 1234
python scripts/seed_metabolomics_data.py --size medium --seed 1234
python scripts/seed_lipidomics_data.py --size medium --seed 1234
python scripts/seed_hts_sar_data.py --size medium --seed 1234
python scripts/seed_signatures_data.py --size medium --seed 1234
```

---

## 4. Size Presets

### Small (Default)

**Use case**: Rapid development, unit tests, quick demos

| Domain | Programs | Experiments | Datasets | Features | Signatures | Compounds |
|--------|----------|-------------|----------|----------|------------|-----------|
| Transcriptomics | 2 | 4 | 8 | ~500 | 2 | - |
| Proteomics | 2 | 4 | 8 | ~300 | 2 | - |
| Metabolomics | 2 | 4 | 8 | ~200 | 2 | 100 |
| Lipidomics | 1 | 2 | 4 | ~150 | 1 | 50 |
| HTS/SAR | 1 | 2 | 1 campaign | - | - | 200 |

**Total time**: ~10-20 seconds

### Medium

**Use case**: Integration testing, realistic demos, local development

| Domain | Programs | Experiments | Datasets | Features | Signatures | Compounds |
|--------|----------|-------------|----------|----------|------------|-----------|
| Transcriptomics | 5 | 15 | 30 | ~5,000 | 10 | - |
| Proteomics | 5 | 15 | 30 | ~3,000 | 10 | - |
| Metabolomics | 5 | 15 | 30 | ~2,000 | 10 | 500 |
| Lipidomics | 3 | 10 | 20 | ~1,500 | 5 | 300 |
| HTS/SAR | 3 | 8 | 5 campaigns | - | - | 1,000 |

**Total time**: ~1-2 minutes

### Large

**Use case**: Performance testing, stress testing, production-like environments

| Domain | Programs | Experiments | Datasets | Features | Signatures | Compounds |
|--------|----------|-------------|----------|----------|------------|-----------|
| Transcriptomics | 10 | 50 | 100 | ~20,000 | 25 | - |
| Proteomics | 10 | 50 | 100 | ~15,000 | 25 | - |
| Metabolomics | 10 | 50 | 100 | ~10,000 | 25 | 2,000 |
| Lipidomics | 5 | 25 | 50 | ~5,000 | 15 | 1,000 |
| HTS/SAR | 5 | 20 | 15 campaigns | - | - | 5,000 |

**Total time**: ~5-10 minutes

---

## 5. Individual Seeders

### Transcriptomics Seeder

**Script**: `scripts/seed_transcriptomics_data.py`

**Generates**:
- Gene expression datasets (RNA-seq, microarray)
- Experiments with tissue/cell type metadata
- Features with HGNC gene symbols, Ensembl IDs
- Differential expression signatures (up/down regulation)

**Realistic distributions**:
- Log-normal expression values
- Batch effects simulation
- Fold-change and p-value correlations

---

### Proteomics Seeder

**Script**: `scripts/seed_proteomics_data.py`

**Generates**:
- Protein abundance datasets (LC-MS/MS, TMT)
- Experiments with sample preparation metadata
- Features with UniProt IDs, protein names
- Phosphorylation and PTM signatures

**Realistic distributions**:
- Log-normal abundance values
- Missing value patterns (MNAR/MCAR)
- Peptide-to-protein aggregation

---

### Metabolomics Seeder

**Script**: `scripts/seed_metabolomics_data.py`

**Generates**:
- Metabolite profiling datasets (GC-MS, LC-MS)
- Experiments with ionization mode metadata
- Features with HMDB IDs, KEGG compound IDs, m/z values
- Pathway-enriched signatures

**Realistic distributions**:
- Log-normal metabolite concentrations
- Retention time clustering
- Adduct and isotope patterns

---

### Lipidomics Seeder

**Script**: `scripts/seed_lipidomics_data.py`

**Generates**:
- Lipid species datasets (LIPID MAPS categories)
- Experiments with extraction method metadata
- Features with systematic lipid names, sum compositions
- Lipid class signatures

**Realistic distributions**:
- Bimodal distributions per lipid class
- Chain length and saturation patterns
- Class-specific abundance ranges

---

### HTS/SAR Seeder

**Script**: `scripts/seed_hts_sar_data.py`

**Generates**:
- HTS campaigns with plate layouts
- Compound libraries with SMILES, properties (MW, LogP, TPSA)
- Dose-response curves (IC50, EC50)
- SAR series with matched molecular pairs

**Realistic distributions**:
- Activity distributions (hit rates 1-5%)
- Druglike property distributions (Lipinski's Rule of Five)
- Hill slopes and curve quality metrics

---

### Signatures Seeder

**Script**: `scripts/seed_signatures_data.py`

**Generates**:
- Cross-omics signatures linking multiple datasets
- Biological programs (pathways, cellular processes)
- Feature membership with weights and directions

**Realistic distributions**:
- Signature sizes (10-500 features)
- Sparse feature overlaps
- Weight distributions (log-normal)

---

## 6. Troubleshooting

### Common Issues

#### Error: "psycopg2.errors.UndefinedTable"

**Cause**: Database schema not initialized

**Fix**:
```bash
alembic upgrade head
```

#### Error: "ModuleNotFoundError: No module named 'faker'"

**Cause**: Missing dependencies

**Fix**:
```bash
pip install -r requirements.txt
```

#### Error: "UNIQUE constraint violation"

**Cause**: Data already exists and seeder is not fully idempotent

**Fix**:
```bash
python scripts/seed_all.py --reset --size small
```

#### Warning: "Seeding took longer than expected"

**Cause**: Large dataset or slow database connection

**Fix**:
- Use `--size small` for faster seeding
- Check database connection (localhost vs remote)
- Ensure indexes are created (run migrations)

#### Error: "Schema drift detected"

**Cause**: Seeder expects columns/tables that don't exist

**Fix**:
1. Check if migrations are up to date: `alembic current`
2. Apply pending migrations: `alembic upgrade head`
3. If seeder is newer than schema, update schema first

---

## 7. For Developers

### Adding a New Seeder

Follow this pattern to create a new domain seeder:

#### 1. Create seeder script

```python
# scripts/seed_my_domain.py
import argparse
from sqlalchemy.orm import Session
from amprenta_rag.database.connection import get_db_session
from amprenta_rag.models.database import MyModel
import random

SIZE_PRESETS = {
    "small": {"entities": 10},
    "medium": {"entities": 50},
    "large": {"entities": 200},
}

def seed_my_domain(size: str, reset: bool, seed: int, dry_run: bool):
    """Seed My Domain data."""
    random.seed(seed)
    config = SIZE_PRESETS[size]
    
    if dry_run:
        print(f"[DRY RUN] Would create {config['entities']} entities")
        return
    
    with get_db_session() as session:
        if reset:
            session.query(MyModel).delete()
            session.commit()
        
        # Generate and insert data
        for i in range(config["entities"]):
            entity = MyModel(name=f"Entity {i}", value=random.random())
            session.add(entity)
        
        session.commit()
        print(f"✓ Created {config['entities']} entities")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Seed My Domain data")
    parser.add_argument("--size", choices=["small", "medium", "large"], default="small",
                        help="Size preset")
    parser.add_argument("--reset", action="store_true",
                        help="Truncate tables before seeding")
    parser.add_argument("--seed", type=int, default=1234,
                        help="Random seed")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview without executing")
    
    args = parser.parse_args()
    seed_my_domain(args.size, args.reset, args.seed, args.dry_run)
```

#### 2. Register in seed_all.py

```python
# scripts/seed_all.py
from seed_my_domain import seed_my_domain

SEEDERS = {
    "my_domain": seed_my_domain,
    # ... other seeders
}
```

#### 3. Test the seeder

```bash
# Test dry run
python scripts/seed_my_domain.py --dry-run --size small

# Test actual seeding
python scripts/seed_my_domain.py --size small --seed 42

# Test idempotency (run twice)
python scripts/seed_my_domain.py --size small --seed 42
python scripts/seed_my_domain.py --size small --seed 42
```

#### 4. Update documentation

Add entry to this file's "Individual Seeders" section.

---

### Best Practices

1. **Use upserts**: Check for existing records before inserting
2. **Bulk inserts**: Use `session.bulk_insert_mappings()` for large datasets
3. **Foreign keys**: Seed dependencies first (Programs → Experiments → Datasets)
4. **Deterministic**: Use seeded RNG for all randomness
5. **Realistic**: Study real data distributions before implementing
6. **Fast**: Target < 10 seconds for small preset
7. **Documented**: Add clear comments for complex generation logic

---

## Reference

- **Plan**: `.cursor/plans/comprehensive_test_data_seeding_suite_87ab8efa.plan.md`
- **Related**: `docs/TESTING_GUIDE.md`, `docs/LOCAL_SETUP.md`
- **Database Schema**: `alembic/versions/`

