import os
import jwt
from datetime import datetime, timedelta, timezone


def generate_jupyter_token(username: str, expires_hours: int = 1) -> str:
    """Generate JWT token for JupyterHub authentication."""
    secret = os.environ.get("JWT_SECRET_KEY", "dev-secret-change-me")
    payload = {
        "username": username,
        "exp": datetime.now(timezone.utc) + timedelta(hours=expires_hours),
        "iat": datetime.now(timezone.utc),
    }
    return jwt.encode(payload, secret, algorithm="HS256")


def get_jupyterhub_url(username: str, base_url: str = "http://localhost:8888") -> str:
    """Generate JupyterHub URL with auth token."""
    token = generate_jupyter_token(username)
    return f"{base_url}/hub/login?token={token}"


