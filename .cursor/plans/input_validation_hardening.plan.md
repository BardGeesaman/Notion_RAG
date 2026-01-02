# Input Validation Hardening

## Overview

Harden input validation across the API and dashboard to prevent injection attacks and ensure data integrity. This complements the Rate Limiting feature just completed.

## Current State

### Findings

1. **No Pydantic strict mode** - 510 BaseModel classes across 48 files, none using `strict=True`
2. **SQL injection risk** - `company_context.py` line 24 uses f-string in SQL:
   ```python
   db.execute(text(f"SET LOCAL app.current_company_id = '{current_user.company_id}'"))
   ```
3. **HTML injection risks** - 41 uses of `unsafe_allow_html=True` in Streamlit dashboards
4. **Comment widget XSS** - `_highlight_mentions()` embeds user content in HTML without sanitization

### Positive Findings

- `_filter_to_sql` in vector_store.py uses parameterized queries correctly
- Most SQLAlchemy queries use ORM patterns (safe)

## Implementation Plan

### Batch 1: SQL Parameterization Audit (Implementor)

1. **Fix company_context.py SQL injection**
   ```python
   # Before (UNSAFE)
   db.execute(text(f"SET LOCAL app.current_company_id = '{current_user.company_id}'"))
   
   # After (SAFE)
   db.execute(text("SET LOCAL app.current_company_id = :company_id"), {"company_id": str(current_user.company_id)})
   ```

2. **Audit all raw SQL usage**
   - Search for `execute(text(f"` and `execute(text(f'` patterns
   - Search for `execute(f"` patterns
   - Verify all use parameterized queries

3. **Add SQL injection test**
   - Test that company_id with SQL injection payload is safely escaped

### Batch 2: HTML Sanitization (Implementor)

**P1 FIX: Use nh3 instead of deprecated bleach (Rust-backed, maintained by Mozilla)**

1. **Install nh3 for HTML sanitization**
   ```
   nh3>=0.2.0
   ```

2. **Create sanitization utility** (`amprenta_rag/utils/sanitize.py`)
   ```python
   """HTML sanitization utilities for user-generated content.
   
   Uses nh3 (Rust-backed, maintained by Mozilla) instead of deprecated bleach.
   """
   import nh3
   from typing import Set
   
   # Allowed tags for rich text (comments, descriptions)
   ALLOWED_TAGS: Set[str] = {'b', 'i', 'u', 'em', 'strong', 'a', 'br', 'p', 'span'}
   ALLOWED_ATTRIBUTES: dict[str, Set[str]] = {
       'a': {'href', 'title'},
       'span': {'style'},
   }
   
   def sanitize_html(content: str) -> str:
       """Sanitize HTML content, removing dangerous tags/attributes.
       
       Args:
           content: User-provided HTML content
       
       Returns:
           Sanitized HTML safe for rendering
       """
       return nh3.clean(
           content,
           tags=ALLOWED_TAGS,
           attributes=ALLOWED_ATTRIBUTES,
       )
   
   def escape_html(content: str) -> str:
       """Escape all HTML in content (for plain text display)."""
       return nh3.clean(content, tags=set())
   
   def escape_for_html_display(content: str) -> str:
       """Escape content for safe display in HTML context.
       
       Renamed from sanitize_for_markdown for clarity (P2 fix).
       """
       return escape_html(content)
   ```

3. **Update comment_widget.py**
   ```python
   from amprenta_rag.utils.sanitize import escape_html
   
   def _highlight_mentions(content: str) -> str:
       # FIRST escape HTML to prevent XSS
       safe_content = escape_html(content)
       # THEN apply highlighting
       pattern = r'@([a-zA-Z0-9_-]+)'
       highlighted = re.sub(
           pattern,
           r'<span style="color: #1f77b4; font-weight: bold;">@\1</span>',
           safe_content
       )
       return highlighted
   ```

4. **Audit all unsafe_allow_html usage**
   - Review 41 instances across 18 files
   - Apply sanitization where user content is rendered
   - Document safe uses (static HTML)

### Batch 3: Pydantic Strict Mode (Implementor)

**P1 FIX: Gradual rollout - don't apply strict to all 510 schemas at once**

1. **Keep BaseSchema non-strict (safe default)**
   ```python
   class BaseSchema(BaseModel):
       """Base schema with common configuration."""
       
       model_config = ConfigDict(
           from_attributes=True,
           strict=False,  # Keep non-strict for backwards compatibility
           json_encoders={
               UUID: str,
               datetime: lambda v: v.isoformat(),
           }
       )
   ```

2. **Create StrictBaseSchema for HIGH-RISK endpoints**
   ```python
   class StrictBaseSchema(BaseModel):
       """Strict schema for security-sensitive endpoints."""
       
       model_config = ConfigDict(
           from_attributes=True,
           strict=True,  # Enforce strict type coercion
           json_encoders={
               UUID: str,
               datetime: lambda v: v.isoformat(),
           }
       )
   ```

3. **Apply StrictBaseSchema to:**
   - `SignRequest` (password input)
   - `ScrapeRequest` (URL input)
   - `CompoundCreate` (SMILES input)
   - Any new security-sensitive schemas

4. **Add field-level validators for common patterns**
   ```python
   from pydantic import field_validator, Field
   import re
   
   # Safe string pattern (no control characters, reasonable length)
   def validate_safe_string(v: str, max_length: int = 10000) -> str:
       if not v:
           return v
       # Remove null bytes and control characters
       v = re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', '', v)
       # Enforce max length
       if len(v) > max_length:
           raise ValueError(f"String exceeds maximum length of {max_length}")
       return v
   ```

3. **Add SMILES validation** (chemistry inputs)
   ```python
   @field_validator('smiles')
   @classmethod
   def validate_smiles_safe(cls, v: str) -> str:
       if not v or not v.strip():
           raise ValueError("SMILES cannot be empty")
       v = v.strip()
       # Basic SMILES character validation
       if not re.match(r'^[A-Za-z0-9@+\-\[\]()=#%.$\\/:]+$', v):
           raise ValueError("Invalid SMILES characters")
       if len(v) > 5000:
           raise ValueError("SMILES too long (max 5000 chars)")
       return v
   ```

4. **Run test suite to identify strict mode breaks**
   - Some tests may need updates for stricter validation

### Batch 4: Request Body Size Limits (Implementor)

**P1 FIX: Add per-route exceptions for upload endpoints (100MB)**

1. **Add request size middleware with route exceptions**
   ```python
   # In amprenta_rag/api/middleware.py
   
   MAX_REQUEST_SIZE = int(os.getenv("MAX_REQUEST_SIZE_MB", "10")) * 1024 * 1024
   MAX_UPLOAD_SIZE = int(os.getenv("MAX_UPLOAD_SIZE_MB", "100")) * 1024 * 1024
   
   # Routes that allow larger uploads
   UPLOAD_ROUTES = {
       "/api/v1/datasets/upload",
       "/api/v1/genomics/alignments/upload",
       "/api/v1/imaging/upload",
       "/api/v1/flow-cytometry/upload",
       "/api/extraction/ocr",
   }
   
   class RequestSizeLimitMiddleware(BaseHTTPMiddleware):
       """Limit request body size to prevent resource exhaustion.
       
       Uses higher limit for known upload endpoints.
       """
       
       async def dispatch(self, request: Request, call_next):
           content_length = request.headers.get("content-length")
           if not content_length:
               return await call_next(request)
           
           size = int(content_length)
           path = request.url.path
           
           # Check if this is an upload route
           is_upload = any(path.startswith(route) for route in UPLOAD_ROUTES)
           max_size = MAX_UPLOAD_SIZE if is_upload else MAX_REQUEST_SIZE
           
           if size > max_size:
               limit_mb = max_size // 1024 // 1024
               return JSONResponse(
                   status_code=413,
                   content={
                       "error": "request_too_large",
                       "detail": f"Request exceeds {limit_mb}MB limit"
                   }
               )
           return await call_next(request)
   ```

2. **Configure in main.py**
   ```python
   app.add_middleware(RequestSizeLimitMiddleware)
   ```

### Batch 5: Tests + Documentation (Implementor)

1. **Create `amprenta_rag/tests/security/test_input_validation.py`**
   - SQL injection prevention tests
   - HTML sanitization tests
   - XSS prevention tests
   - Request size limit tests

2. **Update ROADMAP**
   - Mark "Input Validation Hardening" as COMPLETE

3. **Update session-memory**

## Key Files

| File | Action |
|------|--------|
| `requirements.txt` | ADD bleach |
| `amprenta_rag/utils/sanitize.py` | CREATE |
| `amprenta_rag/auth/company_context.py` | FIX SQL injection |
| `amprenta_rag/api/schemas.py` | ADD strict mode + validators |
| `amprenta_rag/api/middleware.py` | ADD request size limit |
| `scripts/dashboard/components/comment_widget.py` | FIX XSS |
| `amprenta_rag/tests/security/test_input_validation.py` | CREATE |

## Configuration

| Variable | Default | Description |
|----------|---------|-------------|
| `MAX_REQUEST_SIZE_MB` | 10 | Maximum request body size in MB |
| `MAX_UPLOAD_SIZE_MB` | 100 | Maximum upload size for file endpoints |

## Effort Estimate

- **Batch 1**: 2 hours (SQL audit + fix)
- **Batch 2**: 2 hours (HTML sanitization)
- **Batch 3**: 2 hours (Pydantic strict mode)
- **Batch 4**: 1 hour (request size limits)
- **Batch 5**: 2 hours (tests + docs)
- **Total**: ~9 hours (1-2 days)

## P2 Observations (Future)

- Content Security Policy (CSP) headers
- Automated SAST scanning in CI
- Input fuzzing tests
- File upload validation (magic bytes, not just extension)
- **Command injection audit** - 38 subprocess calls found (Reviewer P2)
- Gradual strict mode expansion to more schemas

