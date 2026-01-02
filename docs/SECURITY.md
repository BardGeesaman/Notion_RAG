# Security Implementation Guide

## Overview

This document describes the security patterns and implementations used in Amprenta RAG.

## Input Validation & Injection Prevention

### SQL Injection Protection

**Pattern**: Use parameterized queries for all database operations.

**Implementation**: `amprenta_rag/auth/company_context.py`

```python
# SECURE: Parameterized query
cursor.execute(
    "SELECT id, name FROM companies WHERE user_id = %s",
    (user_id,)
)

# INSECURE: String concatenation (never do this)
# cursor.execute(f"SELECT id FROM companies WHERE user_id = {user_id}")
```

**Key Points**:
- Always use `%s` placeholders with parameter tuples
- Never use f-strings or `.format()` for SQL queries
- Apply to all raw SQL operations

### HTML Sanitization

**Pattern**: Sanitize all user-generated content before storage or display.

**Implementation**: `amprenta_rag/utils/sanitize.py`

```python
import nh3

def sanitize_html(content: str) -> str:
    """Sanitize HTML content using nh3 library."""
    return nh3.clean(content)
```

**Usage**:
```python
# Before storing user content
clean_content = sanitize_html(user_input)
```

**Key Points**:
- Use `nh3` library (Rust-based, secure)
- Sanitize before storage, not just display
- Apply to all user-generated HTML content

### Schema Validation

**Pattern**: Use Pydantic strict mode for high-risk endpoints.

**Implementation**: `amprenta_rag/api/schemas.py`

```python
from pydantic import BaseModel, ConfigDict

class StrictBaseSchema(BaseModel):
    """Base schema with strict validation enabled."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

class CompoundCreate(StrictBaseSchema):
    smiles: str
    name: Optional[str] = None
```

**Key Points**:
- `str_strip_whitespace=True`: Remove leading/trailing whitespace
- `extra='forbid'`: Reject unknown fields
- Gradual rollout to high-risk endpoints first

### Request Size Limits

**Pattern**: Limit request payload sizes to prevent DoS attacks.

**Implementation**: `amprenta_rag/api/middleware.py`

```python
class RequestSizeLimitMiddleware:
    def __init__(self, app, max_upload_size: int = 10 * 1024 * 1024):  # 10MB
        self.app = app
        self.max_upload_size = max_upload_size

    async def __call__(self, scope, receive, send):
        if scope["type"] == "http":
            content_length = int(scope.get("headers", {}).get(b"content-length", b"0"))
            if content_length > self.max_upload_size:
                # Return 413 Payload Too Large
                ...
```

**Configuration**:
- Standard endpoints: 10MB limit
- File upload endpoints: 100MB limit
- Configurable per endpoint as needed

## Authentication & Authorization

### Rate Limiting

**Pattern**: Per-user and per-IP rate limiting on sensitive endpoints.

**Implementation**: slowapi middleware with database-backed storage

**Limits**:
- Signatures endpoint: 5 requests/minute
- Share-links endpoint: 30 requests/minute
- Authentication endpoints: Account lockout after 5 failed attempts (15-minute lockout)

### Account Lockout

**Pattern**: Temporary account lockout after repeated failed authentication attempts.

**Implementation**:
- Track failed attempts in database
- 5 failures → 15-minute lockout
- Clear counter on successful login

## Secrets Management

### Environment Variables

**Pattern**: Use AWS Secrets Manager for production, `.env` files for development.

**Implementation**: `amprenta_rag/config/secrets.py`

**Key Points**:
- Never hardcode secrets in source code
- Use `.env.example` template for development setup
- Rotate secrets quarterly or after team changes
- See `docs/SECRETS_ROTATION.md` for procedures

### GitHub Secrets

**Pattern**: Store CI/CD secrets in GitHub repository secrets.

**Implementation**: See `docs/GITHUB_SECRETS.md`

**Required Secrets**:
- `OPENAI_API_KEY`: LLM integration
- `POSTGRES_PASSWORD`: Database access
- `JWT_SECRET_KEY`: Authentication tokens
- `LIGHTSAIL_SSH_KEY`: Deployment access

## Network Security

### HTTPS/TLS

**Requirements**:
- All production traffic over HTTPS
- TLS 1.2+ minimum
- Proper certificate validation

### CORS Configuration

**Pattern**: Restrict cross-origin requests to trusted domains.

**Implementation**: FastAPI CORS middleware

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["https://yourdomain.com"],  # Specific domains only
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE"],
    allow_headers=["*"],
)
```

### SSRF Protection

**Pattern**: Use safe HTTP client for all external requests.

**Implementation**: `amprenta_rag/utils/safe_requests.py`

```python
from amprenta_rag.utils.safe_requests import safe_get, safe_post

# Safe external API calls with allowlist validation
response = safe_get(
    url="https://api.trusted-service.com/data",
    allowed_domains=["api.trusted-service.com"],
    timeout=30
)

# Automatic protection against:
# - Private IP ranges (RFC 1918)
# - Localhost/loopback addresses
# - Cloud metadata endpoints
# - Non-HTTP/HTTPS schemes
```

**Key Features**:
- URL allowlist validation
- Private IP range blocking
- Timeout enforcement
- Request size limits
- Redirect following controls

## Security Headers

**Implemented Headers**: All security headers are configured via `SecurityHeadersMiddleware` in `amprenta_rag/api/middleware.py`:

- `Content-Security-Policy`: `"default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'"` - Prevent XSS attacks
- `Strict-Transport-Security`: `"max-age=31536000; includeSubDomains"` - Force HTTPS (production only)
- `X-Frame-Options`: `"DENY"` - Prevent clickjacking
- `X-Content-Type-Options`: `"nosniff"` - Prevent MIME sniffing
- `Referrer-Policy`: `"strict-origin-when-cross-origin"` - Control referrer information
- `X-XSS-Protection`: `"1; mode=block"` - Legacy XSS protection

**Key Features**:
- Environment-aware HSTS (production-only to avoid localhost issues)
- Streamlit XSRF protection with meta tags
- FastAPI and Streamlit dual configuration

## Monitoring & Logging

### Security Logging

**Implementation**: `amprenta_rag/utils/security_logger.py`

**Structured Security Logger**:
```python
from amprenta_rag.utils.security_logger import log_security_event

# Log authentication events
log_security_event(
    event_type="auth_failure",
    user_id=user_id,
    ip_address=client_ip,
    details={"reason": "invalid_password", "endpoint": "/login"}
)

# Log access control violations
log_security_event(
    event_type="access_denied",
    user_id=user_id,
    resource="sensitive_endpoint",
    details={"required_role": "admin", "user_role": "user"}
)
```

**Monitored Events**:
- Failed authentication attempts
- Rate limit violations
- Account lockouts
- Unusual access patterns
- Input validation failures
- Authorization failures
- Privilege escalation attempts

**Log Format**: Structured JSON with timestamp, event type, user context, and security-relevant details.

### Audit Trail

**Pattern**: Log all sensitive operations with user context.

**Implementation**: Integrated with Provenance Ledger

**Tracked Operations**:
- Data access/modification
- Configuration changes
- User management
- Secret rotations

## Vulnerability Management

### Dependency Scanning

**Tools**:
- `pip-audit` for Python package vulnerabilities
- GitHub Dependabot for automated updates
- `detect-secrets` for secret scanning

**Process**:
1. Weekly vulnerability scans
2. Immediate patching for critical CVEs
3. Monthly dependency updates

### Code Review

**Security Checklist**:
- [ ] Input validation on all user inputs
- [ ] Parameterized queries for database operations
- [ ] Proper authentication/authorization checks
- [ ] No hardcoded secrets
- [ ] Appropriate error handling (no information leakage)

## Incident Response

### Security Breach Response

1. **Immediate Actions**:
   - Identify scope of compromise
   - Rotate all potentially affected secrets
   - Block malicious IP addresses
   - Notify team immediately

2. **Investigation**:
   - Review access logs
   - Identify attack vectors
   - Assess data exposure

3. **Recovery**:
   - Patch vulnerabilities
   - Update security controls
   - Document lessons learned

### Contact Information

**Security Team**: [Contact details]
**Emergency Escalation**: [Emergency contacts]

## Compliance

### Data Protection

**GDPR Compliance**:
- Data export functionality (`docs/SECRETS_ROTATION.md`)
- Data deletion workflows
- User consent tracking
- Data processing audit logs

### Industry Standards

**Alignment with**:
- OWASP Top 10 security risks
- NIST Cybersecurity Framework
- ISO 27001 security controls

## Security Testing

### Regular Testing

**Automated**:
- Static code analysis (security rules)
- Dependency vulnerability scanning
- Secret detection in commits

**Manual**:
- Quarterly security reviews
- Annual penetration testing
- Code review security checklist

## Security Audits

### OWASP Top 10 Audit (2025-01-02)

A comprehensive security audit was completed covering all OWASP Top 10 2021 categories. See `docs/OWASP_AUDIT_REPORT.md` for:

- Detailed findings for each category
- Test coverage summary (56 security tests)
- Remediation roadmap and priorities
- Compliance status and metrics

**Key Outcomes**:
- 7/10 categories fully compliant
- 3/10 categories partial compliance with remediation plan
- 56 security tests added (100% pass rate)
- Security logging and SSRF protection implemented

## References

- [OWASP Top 10](https://owasp.org/www-project-top-ten/)
- [NIST Cybersecurity Framework](https://www.nist.gov/cyberframework)
- [Pydantic Security](https://docs.pydantic.dev/latest/concepts/validators/)
- [FastAPI Security](https://fastapi.tiangolo.com/tutorial/security/)
- [OWASP Audit Report](docs/OWASP_AUDIT_REPORT.md)
