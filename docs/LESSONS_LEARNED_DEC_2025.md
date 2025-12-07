# Lessons Learned - December 7, 2025

## What Went Wrong

### 1. Two Diverging Repositories
- **Issue:** `Documents/Notion RAG` (cloud) and `Documents/RAG` (local) diverged
- **Impact:** Code changes in one repo weren't reflected in the other
- **Root cause:** Manual copy instead of proper git workflow

### 2. Nested Folders Not Synced
- **Issue:** When syncing repos, nested module folders were missing:
  - `amprenta_rag/agent/`
  - `amprenta_rag/domain/`
  - `amprenta_rag/analysis/pathway/`
  - `amprenta_rag/analysis/program_maps/`
  - And 10+ more
- **Impact:** Dashboard pages failed to import, broken functionality
- **Root cause:** Git worktree/cherry-pick didn't carry forward all directories

### 3. Context Files in .gitignore
- **Issue:** Important context files were excluded from git:
  ```
  ChatGPT context/
  md_archive/
  ```
- **Impact:** AI agents started sessions without critical context/roadmap info
- **Root cause:** Files added to .gitignore for cloud sync reasons

### 4. No Commits During Bulk Changes
- **Issue:** 50+ files modified without intermediate commits
- **Impact:** When things broke, no clean checkpoint to return to
- **Root cause:** Rushing through batch changes without testing

### 5. Dashboard/Backend Version Mismatch
- **Issue:** Restored old dashboard expecting Notion functions that were just stubbed
- **Impact:** Cascading import errors across all pages
- **Root cause:** Mixing code from different points in time

---

## Recommendations for Agents

### For ALL Agents:

1. **Commit frequently**
   - After each successful batch of changes
   - Use descriptive commit messages
   - Create checkpoint before major changes

2. **Test after each change**
   - Run imports: `python -c "from module import *"`
   - Run dashboard: Check pages load
   - Run API: Check endpoints respond

3. **Check git status before starting**
   - `git status` - Are there uncommitted changes?
   - `git log -1` - What was the last commit?
   - Are you on the right branch?

4. **Verify file existence before modifying**
   - Don't assume modules exist
   - Check imports work before changing them

### For Implementor Agent:

1. **Before making changes:**
   ```bash
   git status
   git stash  # if needed
   ```

2. **After each batch:**
   ```bash
   python -c "from amprenta_rag.module import *"  # test imports
   git add -A && git commit -m "Batch X: description"
   ```

3. **If something breaks:**
   ```bash
   git diff  # see what changed
   git checkout -- file.py  # revert single file
   git reset --hard HEAD  # nuclear option
   ```

### For Reviewer Agent:

1. **Verify imports work** before approving
2. **Check for missing dependencies** 
3. **Confirm git checkpoint exists**

### For Architect Agent:

1. **Maintain single source of truth** for codebase location
2. **Document what's NOT in git** (credentials, large files, context)
3. **Plan for cross-workstation workflows**

---

## Cross-Workstation Workflow

### Setup (One Time)

1. **Consolidate to ONE repository**
   - Primary: `Documents/RAG`
   - Archive: `Documents/Notion RAG` (read-only reference)

2. **Remove context from .gitignore** (or create sync script)
   ```bash
   # Option A: Remove from .gitignore
   # Edit .gitignore, remove:
   # ChatGPT context/
   # md_archive/
   
   # Option B: Create sync script (see below)
   ```

3. **Document non-git files**
   - Credentials (copy manually)
   - Large data files (copy or re-download)
   - Reference genomes (copy or re-download)

### Daily Workflow

#### Starting Work (ANY Workstation)
```bash
cd ~/Documents/RAG
git fetch origin
git status  # check for conflicts
git pull origin main
```

#### During Work
```bash
# Commit frequently
git add -A
git commit -m "Description of changes"
```

#### Ending Work Session
```bash
git add -A
git commit -m "End of session: summary"
git push origin main
```

### Sync Script for Non-Git Files

Create `scripts/sync_workspace.sh`:
```bash
#!/bin/bash
# Sync non-git files between workstations

SOURCE="${1:-/Volumes/Backup/RAG-sync}"
DEST="${2:-.}"

# Context files (if not in git)
rsync -av "$SOURCE/ChatGPT context/" "$DEST/ChatGPT context/"
rsync -av "$SOURCE/md_archive/" "$DEST/md_archive/"

# Credentials (sensitive)
rsync -av "$SOURCE/credentials/" "$DEST/credentials/"

# Large data files (optional)
# rsync -av "$SOURCE/reference_data/" "$DEST/reference_data/"
# rsync -av "$SOURCE/salmon_index/" "$DEST/salmon_index/"

echo "Sync complete!"
```

---

## Checklist: Before Starting Any Session

- [ ] `git pull` - Get latest changes
- [ ] `git status` - Clean working directory?
- [ ] Dashboard loads? `streamlit run scripts/run_dashboard.py`
- [ ] API runs? `uvicorn amprenta_rag.api.main:app`
- [ ] Context files present? `ls context/MASTER_CONTEXT*.md`

---

## Files That Should Be in Git (But Weren't)

| File/Folder | Purpose | Action Taken |
|-------------|---------|--------------|
| `amprenta_rag/agent/` | Chat agent | Now synced |
| `amprenta_rag/domain/` | Domain models | Now synced |
| `amprenta_rag/analysis/pathway/` | Pathway analysis | Now synced |
| `amprenta_rag/analysis/program_maps/` | Program maps | Now synced |
| `scripts/dashboard/` | 27-page dashboard | Now synced |

## Files That Are NOT in Git (By Design)

| File/Folder | Reason | Sync Method |
|-------------|--------|-------------|
| `credentials/` | Sensitive | Manual copy |
| `reference_data/` | Large (1GB+) | Manual copy or re-download |
| `salmon_index/` | Large, generated | Re-generate or copy |
| `geo_cache/` | Cache, temporary | Auto-downloads |
| `.env` | Sensitive | Manual copy |
| `ChatGPT context/` | Historical | Removed from .gitignore |
| `md_archive/` | Historical | Removed from .gitignore |

---

## Summary

**Human oversight caught what automation missed.**

The key insight: Git is great for code, but workspace setup requires explicit documentation and verification of ALL components - not just what's tracked in git.

---

*Document created: December 7, 2025*
*Last updated: December 7, 2025*

