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


def get_jupyterhub_url(
    username: str,
    base_url: str = "http://localhost:8888",
    experiment_id: str = None,
    dataset_id: str = None,
    compound_id: str = None
) -> str:
    """Generate JupyterHub URL with auth token and optional context."""
    token = generate_jupyter_token(username)
    url = f"{base_url}/hub/login?token={token}"

    # Add context params for notebooks to use
    if experiment_id:
        url += f"&experiment_id={experiment_id}"
    if dataset_id:
        url += f"&dataset_id={dataset_id}"
    if compound_id:
        url += f"&compound_id={compound_id}"

    return url


