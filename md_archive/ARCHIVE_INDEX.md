# Archive Index

**Last Updated**: December 9, 2025  
**Purpose**: Historical documentation preserved for reference

This archive contains completed work, historical status reports, and superseded documentation.

---

## What's Archived

### Migration & Transition (22 files)
Historical documentation from the Notion → Postgres migration:
- Notion migration plans, status, and completion reports
- Postgres migration strategies and transition guides
- Migration timeout resolution documents

**Status**: Migration complete as of December 2025 [[memory:11923692]]

### Completion Reports (21 files)
Implementation completion documents:
- Feature linking, genomics pipeline, signature integration
- Tier 3/4 priority completions
- Gap filling, protocol implementation
- Optional enhancements completion

### Test Results (13 files)
Test documentation and results:
- Repository ingestion tests (GEO, ENA, PRIDE, MetaboLights)
- Tier 3 testing suite results
- Post-incident verification reports

### Historical Implementation (17 files)
Implementation-specific historical docs:
- ALS study searches and inspections
- GEO repository rewrites and fixes
- MetaboLights extraction development
- Feature extraction optimization learnings

### Redundant Status (15 files)
Superseded by current documentation:
- Old status files (ACTUAL_STATUS, CURRENT_STATUS)
- Comprehensive roadmap status (superseded by UNIFIED_STRATEGIC_ROADMAP)
- API implementation status
- Repository/testing status snapshots

### Notion-Specific (2 files)
Notion implementation docs (now optional):
- Database setup (now Postgres-first)
- Notion disabled by default configuration

---

## Total Archived: 91 files

---

## Active Documentation Location

**See**: `/docs/` directory (47 active files)

### Categories of Active Docs:
- **User Guides**: API Reference, Configuration, Usage Examples, User Guide
- **Feature Guides**: Auto-Linking, Feature Caching, Metadata Editing, Quality Checks, Statistical Analysis, Visualizations
- **Setup**: Local Setup, Install Postgres, Workstation Setup, Deployment Guide
- **Architecture**: Architecture, FastAPI API, Cross-Omics Reasoning, Dashboard Enhancements
- **Repository Guides**: Repository Import Guide, Repository Ingestion Guide, Quick Start Repository Import
- **Testing**: Testing Guide, Testing Strategy
- **Protocols**: Master Bioinformatics Protocol, Genomics Pipeline Protocol, Salmon Installation

---

## Why These Were Archived

**Criteria for Archiving**:
1. ✅ Work is complete and superseded by current implementation
2. ✅ Status reports that are no longer current
3. ✅ Test results from completed validation efforts
4. ✅ Migration documentation after successful transition
5. ✅ Redundant status files replaced by unified documentation

**Not Archived**:
- Current user-facing guides
- Active implementation protocols
- Strategic planning documents (in `/context/`)
- Current testing/deployment guides

---

## Accessing Archive

All archived files remain available in `/md_archive/` for historical reference.

**To view an archived file**:
```bash
# Example
cat /Users/bard/Documents/RAG/md_archive/NOTION_MIGRATION_FINAL_STATUS.md
```

---

## Archive Organization

Files are organized by original purpose:

```
md_archive/
├── Migration files (COMPLETE_*, NOTION_MIGRATION_*, POSTGRES_MIGRATION_*)
├── Completion reports (*_COMPLETE.md, FINAL_*, PRIORITY_*)
├── Test results (*_TEST_*, TIER3_TESTING_*)
├── Historical implementation (GEO_*, METABOLIGHTS_*, ALS_*)
├── Status snapshots (*_STATUS.md, *_SUMMARY.md)
└── This index (ARCHIVE_INDEX.md)
```

---

**Questions?** See active documentation in `/docs/` or current roadmap in `/context/UNIFIED_STRATEGIC_ROADMAP.md`

