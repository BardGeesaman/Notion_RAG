"""Custom middleware for API security."""
import asyncio
import os
import logging
from typing import Set
from fastapi import Request
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import JSONResponse

logger = logging.getLogger(__name__)

REQUEST_TIMEOUT_SECONDS = int(os.getenv("REQUEST_TIMEOUT_SECONDS", "60"))

# Request size limits
MAX_REQUEST_SIZE = int(os.getenv("MAX_REQUEST_SIZE_MB", "10")) * 1024 * 1024
MAX_UPLOAD_SIZE = int(os.getenv("MAX_UPLOAD_SIZE_MB", "100")) * 1024 * 1024

# Routes that allow larger uploads
UPLOAD_ROUTES: Set[str] = {
    "/api/v1/datasets/upload",
    "/api/v1/genomics/alignments/upload",
    "/api/v1/imaging/upload",
    "/api/v1/flow-cytometry/upload",
    "/api/extraction/ocr",
    "/api/v1/flow-cytometry/fcs",
}


class TimeoutMiddleware(BaseHTTPMiddleware):
    """Timeout long-running requests to prevent Slowloris attacks.
    
    Slowloris attacks work by opening many connections and sending
    partial requests slowly, tying up server resources. This middleware
    enforces a maximum request duration.
    
    Args:
        app: The FastAPI application
        timeout: Maximum request duration in seconds (default from env)
    """
    
    def __init__(self, app, timeout: int = REQUEST_TIMEOUT_SECONDS):
        super().__init__(app)
        self.timeout = timeout
        logger.info(f"TimeoutMiddleware initialized with {timeout}s timeout")
    
    async def dispatch(self, request: Request, call_next):
        """Process request with timeout enforcement."""
        try:
            return await asyncio.wait_for(
                call_next(request),
                timeout=self.timeout
            )
        except asyncio.TimeoutError:
            client = request.client.host if request.client else "unknown"
            logger.warning(
                f"Request timeout: {request.method} {request.url.path} "
                f"from {client} exceeded {self.timeout}s"
            )
            return JSONResponse(
                status_code=504,
                content={
                    "error": "request_timeout",
                    "detail": f"Request exceeded maximum duration of {self.timeout} seconds",
                }
            )


class RequestSizeLimitMiddleware(BaseHTTPMiddleware):
    """Limit request body size to prevent resource exhaustion.
    
    Uses higher limit for known upload endpoints to allow file uploads
    while protecting other endpoints from oversized requests.
    
    Attributes:
        default_limit: Default max size for regular requests (10MB)
        upload_limit: Max size for upload routes (100MB)
    """
    
    def __init__(self, app, default_limit: int = MAX_REQUEST_SIZE, 
                 upload_limit: int = MAX_UPLOAD_SIZE):
        super().__init__(app)
        self.default_limit = default_limit
        self.upload_limit = upload_limit
    
    async def dispatch(self, request: Request, call_next):
        content_length = request.headers.get("content-length")
        if not content_length:
            return await call_next(request)
        
        try:
            size = int(content_length)
        except ValueError:
            return await call_next(request)
        
        path = request.url.path
        
        # Check if this is an upload route
        is_upload = any(path.startswith(route) for route in UPLOAD_ROUTES)
        max_size = self.upload_limit if is_upload else self.default_limit
        
        if size > max_size:
            limit_mb = max_size // 1024 // 1024
            logger.warning(
                f"Request too large: {size} bytes from {request.client.host} "
                f"to {path} (limit: {limit_mb}MB)"
            )
            return JSONResponse(
                status_code=413,
                content={
                    "error": "request_too_large",
                    "detail": f"Request body exceeds {limit_mb}MB limit",
                    "limit_bytes": max_size,
                    "request_bytes": size,
                }
            )
        return await call_next(request)


class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """Add security headers to all responses."""
    
    async def dispatch(self, request: Request, call_next):
        response = await call_next(request)
        
        # Always add these headers
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "SAMEORIGIN"
        response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
        response.headers["Permissions-Policy"] = "geolocation=(), microphone=(), camera=()"
        
        # CSP - permissive for API but protective
        response.headers["Content-Security-Policy"] = (
            "default-src 'self'; "
            "script-src 'self' 'unsafe-inline' 'unsafe-eval'; "
            "style-src 'self' 'unsafe-inline'; "
            "img-src 'self' data: blob:; "
            "connect-src 'self' *"
        )
        
        # HSTS only in production (check ENVIRONMENT env var)
        import os
        if os.environ.get("ENVIRONMENT", "dev") != "dev":
            response.headers["Strict-Transport-Security"] = "max-age=31536000; includeSubDomains"
        
        return response
