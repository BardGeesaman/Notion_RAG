# OWASP Top 10 Security Audit

Systematic review of all 10 OWASP categories with remediation. Estimated effort: 3 weeks.

## Approach: Scripted Audit

Use automated grep/script-based detection for bulk endpoint analysis, manual review only for flagged items. This reduces timeline from 4.5 weeks (full manual) to 3 weeks.

## Current Security Posture

Already implemented (can mark as PASS with verification):
- **A03 Injection**: SQL parameterization, HTML sanitization (`sanitize.py`)
- **A07 Auth Failures**: Rate limiting, account lockout (`lockout.py`)  
- **A10 SSRF**: URL validation in extraction.py (private IP blocking)
- Security headers, secrets management

## Audit Structure

Each category will produce:
1. **Audit checklist** with pass/fail for each control
2. **Remediation code** for gaps identified
3. **Test coverage** for security controls
4. **Documentation update** to SECURITY.md

---

## Batch 1: Access Control (A01) - Days 1-5

**A01: Broken Access Control** - Most critical, highest impact

### Approach: Scripted Detection
1. **Automated scan**: Grep all `@router` endpoints, flag those missing `get_current_user`
2. **Manual review**: Only flagged endpoints (estimated 10-20% of 434)
3. **IDOR tests**: Focus on entity-access endpoints (datasets, compounds, experiments)

### Script: `scripts/audit_auth_deps.py`
```python
# Scan all routers for endpoints missing auth dependency
# Output: CSV of endpoint, file, has_auth, needs_review
```

### Deliverables
1. **Endpoint Authorization Matrix** - Auto-generated CSV with manual annotations
2. **Missing Auth Remediation** - Add `Depends(get_current_user)` where missing
3. **IDOR Test Suite** - Tests for top 20 high-risk endpoints
4. **Horizontal Privilege Tests** - User A accessing User B's data

### Key Files
- `amprenta_rag/api/routers/` (84 files, 434 endpoints)
- `amprenta_rag/auth/permissions.py` (entity-level checks: can_view_entity, can_edit_entity)
- `amprenta_rag/api/dependencies.py` (get_current_user dependency)

---

## Batch 2: Cryptography & Authentication (A02, A07) - Days 6-7

**A02: Cryptographic Failures**
- JWT token signing algorithm audit (HS256 vs RS256)
- Password hashing strength (bcrypt rounds)
- Secrets in logs/error messages

**A07: Identification & Authentication Failures**
- Session management audit
- Password policy enforcement
- Token expiration verification

### Key Files
- `amprenta_rag/auth/password.py`
- `amprenta_rag/auth/session.py`
- `amprenta_rag/auth/signatures.py`

---

## Batch 3: Injection Hardening (A03) - Days 8-9

**A03: Injection** - Expand beyond current coverage

### Scope
- Audit 73 outbound HTTP calls for SSRF (beyond extraction.py)
- Command injection in subprocess calls (MAGeCK, fpocket, Vina)
- Path traversal in file operations
- NoSQL/ORM injection patterns

### Key Files to Audit
- `amprenta_rag/ingestion/` - external API calls (33 files with requests.get/post)
- `amprenta_rag/structural/` - subprocess calls (fpocket, vina)
- `amprenta_rag/crispr/` - MAGeCK subprocess
- `amprenta_rag/api/routers/imaging.py` - file uploads

### Specific Concerns
- `connectivity/string_client.py` - STRING API calls
- `ingestion/genomics/vep_client.py` - VEP API calls
- `ingestion/papers/*.py` - PubMed, OpenAlex, Semantic Scholar
- `sync/adapters/*.py` - ChEMBL, PubChem, UniProt

---

## Batch 4: Design & Configuration (A04, A05) - Days 10-11

**A04: Insecure Design**
- Business logic flaws (mass assignment, race conditions)
- Trust boundary review
- Fail-secure defaults

**A05: Security Misconfiguration**
- CORS policy audit
- Debug mode in production check
- Error response information leakage
- Default credentials audit

### Key Files
- `amprenta_rag/api/main.py` - CORS, debug settings
- `amprenta_rag/config/` - configuration audit
- `.env.example` - default values review

### Specific Checks
- CORS `allow_origins` not set to `*` in production
- `DISABLE_AUTH` env var not enabled in production
- Error responses don't leak stack traces
- No hardcoded credentials in code

---

## Batch 5: Components & Integrity (A06, A08) - Days 12-13

**A06: Vulnerable Components**
- Address 4 known CVEs:
  - biopython
  - nbconvert
  - pdfminer-six
  - py (pytest dependency)
- Dependency update plan
- pip-audit integration in CI

**A08: Software & Data Integrity**
- Pickle/joblib deserialization audit (ML models)
- CI/CD pipeline security
- Unsigned artifact detection

### Key Files
- `requirements.txt`
- `amprenta_rag/ml/registry.py` - model loading with joblib
- `.github/workflows/ci.yml` - CI security

### Remediation
- Add `pip-audit` to CI pipeline
- Document acceptable risk for CVEs without patches
- Add model signature verification (stretch goal)

---

## Batch 6: Logging & Monitoring (A09) - Days 14-15

**A09: Security Logging & Monitoring Failures**

### Deliverables
1. **Security Logger Module** - `amprenta_rag/auth/security_logger.py`
2. **Event Taxonomy** - Define loggable security events
3. **Alerting Integration** - Critical event notifications
4. **Log Review Dashboard** - Security event viewer (Streamlit)

### Events to Log
- Authentication failures/successes
- Authorization failures (403s)
- Rate limit violations
- Input validation failures
- Account lockouts
- Suspicious patterns (enumeration, scanning)

### Implementation
```python
# amprenta_rag/auth/security_logger.py
import logging
from enum import Enum

class SecurityEvent(Enum):
    AUTH_SUCCESS = "auth_success"
    AUTH_FAILURE = "auth_failure"
    AUTHZ_FAILURE = "authz_failure"
    RATE_LIMIT = "rate_limit"
    ACCOUNT_LOCKED = "account_locked"
    VALIDATION_FAILURE = "validation_failure"
    SUSPICIOUS_ACTIVITY = "suspicious_activity"

security_logger = logging.getLogger("security")

def log_security_event(
    event: SecurityEvent,
    user_id: Optional[UUID],
    ip_address: str,
    details: dict
) -> None:
    """Log security event with structured data."""
    ...
```

---

## Batch 7: Documentation & Final Review - Day 16

1. Update `docs/SECURITY.md` with audit findings
2. Create `docs/OWASP_AUDIT_REPORT.md` with compliance matrix
3. Update `docs/ROADMAP.md` to mark audit complete
4. Run full security test suite
5. Git commit and push

---

## Timeline Summary (3 weeks)

| Batch | Scope | Days |
|-------|-------|------|
| 1 | A01 Access Control (scripted + manual) | 1-5 |
| 2 | A02/A07 Crypto & Auth | 6-7 |
| 3 | A03 Injection (73 HTTP calls) | 8-9 |
| 4 | A04/A05 Design & Config | 10-11 |
| 5 | A06/A08 Components & CVEs | 12-13 |
| 6 | A09 Security Logging | 14-15 |
| 7 | Documentation | 16 |

---

## Test Coverage Target

| Category | New Tests | Focus Areas |
|----------|-----------|-------------|
| A01 Access Control | 10-12 | IDOR, horizontal privilege (top 20 endpoints) |
| A02/A07 Crypto/Auth | 6-8 | JWT, password, session |
| A03 Injection | 8-10 | SSRF, command injection |
| A04/A05 Design/Config | 4-6 | CORS, error leakage |
| A06/A08 Components | 3-4 | CVE verification |
| A09 Logging | 4-6 | Event logging |
| **Total** | **35-46** | |

---

## Success Criteria

- [ ] All 10 OWASP categories audited with documented findings
- [ ] Critical/High severity issues remediated
- [ ] Medium severity issues documented in backlog
- [ ] Security test suite expanded by 35+ tests
- [ ] SECURITY.md updated with audit results
- [ ] OWASP_AUDIT_REPORT.md created with compliance matrix

---

## Risk Assessment

### Timeline Risk (MITIGATED)
- 434 endpoints audited via scripted detection
- Manual review limited to flagged endpoints (~50-80)
- 5 days allocated for Batch 1 (vs original 3)

### Scope Risk  
- 73 HTTP clients may have inconsistent SSRF protection
- Mitigation: Create shared `safe_request()` utility
- Fallback: Document gaps in backlog if time-constrained

### Dependency Risk
- CVEs may not have patches available
- Mitigation: Document accepted risk, pin versions
- Fallback: Add to monitoring pipeline

