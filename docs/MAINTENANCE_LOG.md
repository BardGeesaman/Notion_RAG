# Maintenance Log

## 2025-12-24 — Dependency + Dead Code Audit

### Scope

- **Security audit**: `pip-audit` (reports under `reports/security_audit.*`)
- **Outdated dependency audit**: `pip list --outdated` (reports under `reports/outdated.*`)
- **Dead code scan**: `vulture` + `ruff F401` (reports under `reports/dead_code_*`)

### Security findings (4 packages)

From `reports/security_audit.md`:

| Package | Version | Finding |
|---|---:|---|
| biopython | 1.86 | CVE-2025-68463 |
| nbconvert | 7.16.6 | CVE-2025-53000 |
| pdfminer-six | 20251107 | GHSA-f83h-ghpp-7wcc |
| py | 1.11.0 | PYSEC-2022-42969 |

**Notes**
- `pip-audit` did not report fix versions for these findings in this environment. We need to confirm available patched versions upstream before upgrading.

### Outdated packages (33 total)

From `reports/outdated_summary.md`:

- **Total outdated**: 33
- **Notable major-version jumps** (higher risk / more testing needed):
  - `numpy` 1.26.4 → 2.2.6
  - (others may exist; see `reports/outdated_summary.md` table)

### Dead code / hygiene findings

From `reports/dead_code_summary.md`:

- **Vulture candidates (>= 60%)**: 680
- **High-confidence candidates (>= 80%)**: 99
- **Unused imports (ruff F401)**: 115

**Notes**
- Many high-confidence hits are "unused variable" warnings (often intentional throwaways).
- Framework-driven code (FastAPI, Pydantic, pytest fixtures) is frequently flagged without a whitelist; we used `scripts/maintenance/vulture_whitelist.py` to reduce false positives, but triage is still required.

## Recommended actions (prioritized)

### 1) CRITICAL — Address security vulnerabilities (4 packages)

- **biopython**: investigate CVE-2025-68463 applicability (XXE via `Bio.Entrez`); upgrade to patched version if available, and validate that any XML parsing is hardened.
- **nbconvert**: CVE-2025-53000 affects Windows search path behavior; if nbconvert is used in Windows environments, prioritize upgrade and/or harden conversion paths.
- **pdfminer-six**: GHSA-f83h-ghpp-7wcc is severe; assess whether `pdfminer-six` is used in production paths. If yes, upgrade and consider sandboxing PDF parsing jobs.
- **py**: PYSEC-2022-42969; this package often comes as an indirect dependency (historically via pytest tooling). Prefer upgrading/removing if not needed.

### 2) HIGH — Major-version updates and dependency drift

- Separate upgrades into **small, testable batches**:
  - Batch A: security-fix upgrades
  - Batch B: major-version upgrades (e.g., numpy 2.x)
  - Batch C: minor/patch updates
- Run full test suite + smoke tests after each batch (API + dashboard + E2E).

### 3) MEDIUM — Dead code triage (99 high-confidence)

- Triage the `reports/dead_code_summary.md` list:
  - Confirm whether flagged endpoints are truly unused or only "dynamically referenced".
  - For "unused variables", decide whether to:
    - keep (intentional, improves readability), or
    - replace with `_` to clarify intent.

### 4) LOW — Unused imports (115)

- Address `ruff F401` findings:
  - Safe mechanical cleanup in most cases.
  - Prefer autofix (`ruff check --fix --select F401`) after confirming no runtime import side-effects.

---

## Actions Completed (2025-12-24)

### Unused Imports Cleanup
- Removed 115 unused imports across 64 files
- Method: ruff check --select F401 --fix
- Verification: 0 F401 errors remaining
- Tests: PASS

### Security Vulnerabilities
- 4 packages identified with CVEs
- Status: No patches available yet (monitoring)
- Packages pinned to latest versions in requirements.txt

### Dead Code Removal
- Triaged 99 high-confidence vulture candidates
- Removed 23 true dead code items (unused args, locals, stubs)
- Skipped 69 false positives (FastAPI routes, Pydantic, fixtures)
- 7 uncertain candidates deferred for manual review
- Files modified: 17
- Tests: PASS

## Actions Completed (2025-01-01)

### Dead Code Final Triage
- Re-scanned codebase with vulture (25 candidates found)
- Triaged 25 remaining candidates
- Removed 8 dead code items:
  - 3 unused upload request schemas (DSCUploadRequest, MSTUploadRequest, SPRUploadRequest)
  - 1 unused error response schema (ErrorResponse)
  - 3 unused plotly imports (make_subplots)
  - 1 unreachable code block
  - Fixed 1 syntax error (indentation)
  - Fixed 1 unreachable else condition
- Added 4 items to whitelist with justification:
  - Signal handler parameters (signum)
  - Stub function parameters (db_id)
  - Conditional dashboard variables (touch_friendly, max_concurrent)
- Remaining candidates: 13 (mostly API schemas that may be used for type hints)
- Codebase significantly cleaner (68% reduction in vulture findings)
- Tests: PASS

