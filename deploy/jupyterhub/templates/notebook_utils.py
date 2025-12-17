"""
Shared utilities for Voila/Jupyter notebooks shipped in `deploy/jupyterhub/templates/`.

These helpers standardize:
- RDKit warning suppression
- API client auto-detection / demo-mode fallback
- Simple UI helpers (demo banner, error display)
- Best-effort RDKit SMILES parsing
"""

from __future__ import annotations


def suppress_rdkit_warnings():
    """Suppress RDKit kekulization and other warnings."""
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        pass

    import logging

    logging.getLogger("amprenta_rag.chemistry").setLevel(logging.WARNING)


def get_api_client(api_url: str = None):
    """Get RAGClient with automatic URL detection and fallback.

    Tries URLs in order:
    1. Provided api_url
    2. AMPRENTA_API_URL or API_URL env vars
    3. http://host.docker.internal:8000
    4. http://localhost:8000

    Returns (client, is_demo_mode) tuple.
    """
    from amprenta_rag.client import RAGClient
    import os

    urls = []
    if api_url:
        urls.append(api_url)
    env_url = os.environ.get("AMPRENTA_API_URL") or os.environ.get("API_URL")
    if env_url:
        urls.append(env_url)
    urls.extend(
        [
            "http://host.docker.internal:8000",
            "http://localhost:8000",
        ]
    )

    for url in urls:
        try:
            client = RAGClient(api_url=url)
            # Test connection
            client.http.get("/health")
            return client, False
        except Exception:
            continue

    return None, True


def demo_mode_banner():
    """Display a banner indicating demo mode is active."""
    from IPython.display import display, HTML

    display(
        HTML(
            """
        <div style="background:#fff3cd; border:1px solid #ffc107; border-radius:4px; padding:10px; margin:10px 0;">
            <strong>⚠️ DEMO MODE</strong> — API unavailable. Showing sample data.
        </div>
    """
        )
    )


def display_error(message: str, details: str = None):
    """Display a styled error message."""
    from IPython.display import display, HTML

    html = '<div style="background:#f8d7da; border:1px solid #f5c6cb; border-radius:4px; padding:10px; margin:10px 0;">'
    html += f"<strong>❌ Error:</strong> {message}"
    if details:
        html += f'<pre style="margin-top:8px; white-space:pre-wrap;">{details}</pre>'
    html += "</div>"
    display(HTML(html))


def mol_from_smiles_safe(smiles: str):
    """Parse SMILES with fallback for kekulization issues."""
    if not smiles:
        return None
    try:
        from rdkit import Chem

        m = Chem.MolFromSmiles(smiles)
        if m is not None:
            return m
        # Fallback for aromatic issues
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return None
        Chem.SanitizeMol(
            m,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        return m
    except Exception:
        return None



