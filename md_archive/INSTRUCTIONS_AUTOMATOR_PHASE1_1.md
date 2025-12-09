# Instructions for Automator Agent - Phase 1.1: Verify Dependencies

## Task Overview

**Agent:** Automator  
**Phase:** 1.1  
**Objective:** Verify dependencies before deleting Notion modules  
**Priority:** Critical - Must complete before Phase 1.2 (deletion)

## Background

Before deleting Notion-specific modules, we need to verify:
1. Which files import each Notion module
2. What specific functions/classes are imported
3. That all imports will be handled in Phase 2 (import removal)

## Modules to Verify

The following Notion modules are candidates for deletion:

1. `amprenta_rag.clients.notion_client`
2. `amprenta_rag.ingestion.notion_pages`
3. `amprenta_rag.ingestion.dataset_notion_utils`
4. `amprenta_rag.ingestion.signature_notion_crud`
5. `amprenta_rag.analysis.pathway_notion_integration`
6. `amprenta_rag.chemistry.notion_integration`
7. `amprenta_rag.migration.dual_write`

## Instructions

### Step 1: Search for Imports

For each module listed above, search the codebase for:
- `from <module> import ...`
- `import <module>`
- `from <module>.<submodule> import ...`

**Search locations:**
- `amprenta_rag/` directory (all Python files)
- `scripts/` directory (all Python files)
- Any test files in `amprenta_rag/tests/`

**Tools to use:**
- `grep` or `ripgrep` for pattern matching
- Search for exact module paths
- Search for partial matches (e.g., `notion_client`, `notion_pages`)

### Step 2: Document Findings

For each Notion module, create a section in your report with:

**Module:** `amprenta_rag.clients.notion_client` (example)

**Files that import it:**
- `amprenta_rag/query/cross_omics/helpers.py` - imports `notion_headers`, `get_page_text`
- `amprenta_rag/ingestion/notion_pages.py` - imports `notion_headers`
- ... (list all files)

**Specific imports used:**
- `notion_headers()` - used in 15 files
- `get_page_text()` - used in 8 files
- `_fetch_chunk_text_property()` - used internally only

**Notes:**
- Any special considerations
- Circular dependencies
- Test files that mock these modules

### Step 3: Verify Phase 2 Coverage

Check that all identified imports are covered in Phase 2 tasks:

**Example verification:**
- ✅ `amprenta_rag/query/cross_omics/helpers.py` → Covered in Phase 2.7
- ✅ `amprenta_rag/ingestion/notion_pages.py` → This file will be deleted, so no action needed
- ⚠️ `amprenta_rag/some_other_file.py` → NOT covered in Phase 2 (needs to be added)

### Step 4: Generate Report

Create a markdown report file: `docs/NOTION_DEPENDENCY_REPORT.md`

**Report structure:**

```markdown
# Notion Module Dependency Report

**Generated:** [Date]
**Purpose:** Verify dependencies before deleting Notion modules

## Summary

- Total modules to delete: 7
- Total files importing Notion modules: [count]
- Files requiring import removal in Phase 2: [count]
- Files safe to ignore (will be deleted): [count]

## Module-by-Module Analysis

### 1. amprenta_rag.clients.notion_client

**Files importing:**
- [list files]

**Functions/classes imported:**
- [list imports]

**Phase 2 coverage:**
- [verify each file is covered]

### 2. amprenta_rag.ingestion.notion_pages
[... repeat for each module]

## Recommendations

- [ ] All imports are covered in Phase 2
- [ ] No unexpected dependencies found
- [ ] Safe to proceed with Phase 1.2 (deletion)

## Files Requiring Attention

[List any files that import Notion modules but are NOT covered in Phase 2 tasks]
```

## Search Commands

Use these commands to search for imports:

```bash
# Search for notion_client imports
grep -r "from amprenta_rag.clients.notion_client" amprenta_rag/ scripts/
grep -r "import.*notion_client" amprenta_rag/ scripts/

# Search for notion_pages imports
grep -r "from amprenta_rag.ingestion.notion_pages" amprenta_rag/ scripts/
grep -r "import.*notion_pages" amprenta_rag/ scripts/

# Search for dataset_notion_utils imports
grep -r "from amprenta_rag.ingestion.dataset_notion_utils" amprenta_rag/ scripts/
grep -r "import.*dataset_notion_utils" amprenta_rag/ scripts/

# Search for signature_notion_crud imports
grep -r "from amprenta_rag.ingestion.signature_notion_crud" amprenta_rag/ scripts/
grep -r "import.*signature_notion_crud" amprenta_rag/ scripts/

# Search for pathway_notion_integration imports
grep -r "from amprenta_rag.analysis.pathway_notion_integration" amprenta_rag/ scripts/
grep -r "import.*pathway_notion_integration" amprenta_rag/ scripts/

# Search for chemistry.notion_integration imports
grep -r "from amprenta_rag.chemistry.notion_integration" amprenta_rag/ scripts/
grep -r "from amprenta_rag.chemistry import.*notion" amprenta_rag/ scripts/

# Search for dual_write imports
grep -r "from amprenta_rag.migration.dual_write" amprenta_rag/ scripts/
grep -r "import.*dual_write" amprenta_rag/ scripts/
```

## Additional Patterns to Search

Also search for:
- Partial imports: `notion_client`, `notion_pages`, `notion_utils`
- Function calls: `notion_headers()`, `get_page_text()`, `fetch_notion_page()`
- Class references: `DualWriteManager`

## Output Requirements

1. **Dependency Report:** `docs/NOTION_DEPENDENCY_REPORT.md`
   - Complete analysis of all imports
   - Verification of Phase 2 coverage
   - Recommendations for proceeding

2. **Summary:** Brief summary in your response:
   - Total files found importing Notion modules
   - Any files NOT covered in Phase 2
   - Recommendation: Safe to proceed or needs attention

## Acceptance Criteria

- [x] All 7 Notion modules searched
- [x] All importing files identified
- [x] All imports documented
- [x] Phase 2 coverage verified
- [x] Report generated and saved
- [x] Any gaps identified and documented

## Next Steps

After completing this task:
- If all dependencies are covered → Proceed to Phase 1.2 (deletion)
- If gaps found → Update Phase 2 plan to include missing files
- Share report with Architect for review before deletion

## Questions?

If you encounter:
- Circular dependencies → Document them
- Unclear imports → List them for Architect review
- Test mocks → Note them (will be handled in Phase 12)

