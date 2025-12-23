# Maintenance Summary (2025-12-24)

## Executive summary

We ran a baseline maintenance audit across dependencies and dead-code hygiene.

- **Security**: 4 vulnerable packages flagged by `pip-audit`.
- **Outdated dependencies**: 33 packages outdated (includes at least one major-version jump: `numpy` 1.x → 2.x).
- **Dead code signals**: 99 high-confidence vulture findings; 115 unused imports (ruff F401).

Artifacts:
- `reports/security_audit.md`, `reports/security_audit.json`
- `reports/outdated_summary.md`, `reports/outdated.json`
- `reports/dead_code_vulture.txt`, `reports/dead_code_summary.md`

## Prioritized action items

### 1) CRITICAL — Security vulnerabilities (4 packages)

Packages flagged:
- **biopython 1.86**: CVE-2025-68463
- **nbconvert 7.16.6**: CVE-2025-53000
- **pdfminer-six 20251107**: GHSA-f83h-ghpp-7wcc
- **py 1.11.0**: PYSEC-2022-42969

Actions:
- Confirm whether each package is used in production/runtime paths vs dev-only.
- Identify patched versions upstream and upgrade in a dedicated PR.
- Add/strengthen tests around the impacted functionality (PDF parsing, nbconvert, Entrez workflows).

### 2) HIGH — Major version outdated packages

Actions:
- Split upgrades into batches and gate with tests.
- Treat major upgrades (notably `numpy` 2.x) as a mini-migration:
  - run unit + integration + E2E
  - validate any scientific stack compatibility (scipy, pandas, statsmodels, etc.)

### 3) MEDIUM — Dead code removal candidates

Signals:
- 99 high-confidence vulture hits (>= 80% confidence)

Actions:
- Triage top offenders first:
  - unused FastAPI endpoints/services (often false positives due to dynamic wiring)
  - unused variables in ingestion/util scripts (likely real)
- For confirmed dead code:
  - remove or consolidate
  - add regression tests where needed

### 4) LOW — Minor version updates + unused imports

Signals:
- 115 unused imports (ruff F401)

Actions:
- Apply `ruff --fix --select F401` in a mechanical cleanup PR (after spot-checking any intentional import side-effects).
- Roll minor/patch upgrades after security/major upgrades stabilize.

## Recommended sequencing (lowest risk)

1. Security upgrades (4 packages) + tests
2. Unused import cleanup (F401) (mechanical)
3. Dead-code triage/removals (small PRs)
4. Major-version upgrades (numpy 2.x and any dependency chain) with full regression testing


