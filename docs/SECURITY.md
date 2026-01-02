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

## Security Headers

**Recommended Headers** (Future implementation):
- `Content-Security-Policy`: Prevent XSS attacks
- `Strict-Transport-Security`: Force HTTPS
- `X-Frame-Options`: Prevent clickjacking
- `X-Content-Type-Options`: Prevent MIME sniffing
- `Referrer-Policy`: Control referrer information

## Monitoring & Logging

### Security Events

**Log These Events**:
- Failed authentication attempts
- Rate limit violations
- Account lockouts
- Unusual access patterns
- Input validation failures

**Implementation**:
```python
import logging
security_logger = logging.getLogger("security")

# Log security events
security_logger.warning(f"Failed login attempt for user {user_id} from IP {client_ip}")
```

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

## References

- [OWASP Top 10](https://owasp.org/www-project-top-ten/)
- [NIST Cybersecurity Framework](https://www.nist.gov/cyberframework)
- [Pydantic Security](https://docs.pydantic.dev/latest/concepts/validators/)
- [FastAPI Security](https://fastapi.tiangolo.com/tutorial/security/)
