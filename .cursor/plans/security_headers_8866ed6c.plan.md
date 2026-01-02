---
name: Security Headers
overview: Add HTTP security headers (CSP, HSTS, X-Frame-Options, X-Content-Type-Options, Referrer-Policy) to FastAPI and Streamlit for defense-in-depth protection against XSS, clickjacking, and MIME sniffing attacks.
todos:
  - id: fastapi-middleware
    content: Add SecurityHeadersMiddleware to FastAPI with all 6 headers
    status: pending
  - id: streamlit-config
    content: Configure Streamlit XSRF protection and CSP meta tag
    status: pending
  - id: tests
    content: Add security header tests (FastAPI + Streamlit)
    status: pending
  - id: documentation
    content: Update SECURITY.md and ROADMAP
    status: pending
---

# Security Headers Implementation

Complete the security hardening sprint by adding HTTP security headers to both FastAPI API and Streamlit dashboard.

## Scope

| Header | Purpose | Value |
|--------|---------|-------|
| Content-Security-Policy | Prevent XSS | `default-src 'self'; script-src 'self' 'unsafe-inline' 'unsafe-eval'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self' *` |
| Strict-Transport-Security | Force HTTPS | `max-age=31536000; includeSubDomains` |
| X-Frame-Options | Prevent clickjacking | `SAMEORIGIN` |
| X-Content-Type-Options | Prevent MIME sniffing | `nosniff` |
| Referrer-Policy | Control referrer info | `strict-origin-when-cross-origin` |
| Permissions-Policy | Restrict browser features | `geolocation=(), microphone=(), camera=()` |

## Architecture

```mermaid
flowchart LR
    subgraph FastAPI[FastAPI API]
        SecurityHeadersMiddleware --> AllRoutes
    end
    subgraph Streamlit[Streamlit Dashboard]
        StreamlitConfig[".streamlit/config.toml"] --> Headers
    end
```

## Implementation

### Batch 1: FastAPI Security Headers Middleware
- Add `SecurityHeadersMiddleware` to [`amprenta_rag/api/middleware.py`](amprenta_rag/api/middleware.py)
- Register in [`amprenta_rag/api/main.py`](amprenta_rag/api/main.py) (after CORS)
- Environment-aware: Skip HSTS in dev mode
- 4 tests for header presence

### Batch 2: Streamlit Configuration
- Update [`.streamlit/config.toml`](.streamlit/config.toml) with `server.enableXsrfProtection`
- Add CSP meta tag via custom HTML component
- 2 tests for Streamlit headers

### Batch 3: Documentation and Finalization
- Update [`docs/SECURITY.md`](docs/SECURITY.md) to mark headers as implemented
- Update ROADMAP to mark Security Headers complete
- Final test suite run

## Key Files

| File | Changes |
|------|---------|
| `amprenta_rag/api/middleware.py` | Add SecurityHeadersMiddleware |
| `amprenta_rag/api/main.py` | Register middleware |
| `.streamlit/config.toml` | Enable XSRF protection |
| `amprenta_rag/tests/security/test_security_headers.py` | New test file |
| `docs/SECURITY.md` | Update documentation |

## CSP Considerations

Streamlit requires relaxed CSP due to:
- Inline scripts for component rendering
- Inline styles for theming
- WebSocket connections for live updates
- Blob URLs for file downloads

We'll use a permissive but still protective policy.

## Effort

- **Estimated**: 3-5 days
- **Batches**: 3
- **Tests**: ~6 new tests

