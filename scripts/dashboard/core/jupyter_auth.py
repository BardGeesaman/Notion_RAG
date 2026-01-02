import os
import jwt
from datetime import datetime, timedelta, timezone


def generate_jupyter_token(username: str, expires_hours: int = 1) -> str:
    """Generate JWT token for JupyterHub authentication."""
    secret = os.environ.get("JWT_SECRET_KEY")
    if not secret:
        raise RuntimeError(
            "JWT_SECRET_KEY not configured. "
            "Generate with: python -c \"import secrets; print(secrets.token_hex(32))\""
        )
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


def get_voila_url(
    username: str,
    dashboard: str = "experiment_dashboard.ipynb",
    context: dict = None,
    base_url: str = "http://localhost:8888"
) -> str:
    """Generate Voila URL with optional context parameters."""
    url = f"{base_url}/user/{username}/voila/render/templates/{dashboard}"
    if context:
        import base64
        import json
        ctx_json = json.dumps(context)
        ctx_b64 = base64.urlsafe_b64encode(ctx_json.encode()).decode()
        url += f"?ctx={ctx_b64}"
    return url

