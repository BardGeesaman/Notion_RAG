"""Custom middleware for API security."""
import asyncio
import os
import logging
from fastapi import Request
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import JSONResponse

logger = logging.getLogger(__name__)

REQUEST_TIMEOUT_SECONDS = int(os.getenv("REQUEST_TIMEOUT_SECONDS", "60"))


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
