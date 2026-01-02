"""
JupyterHub configuration for a simple single-host deployment.

Notes:
- Uses DummyAuthenticator (dev only).
- Persists hub DB and user home directories via docker-compose volumes.
"""

import os
import sys

sys.path.insert(0, "/etc/jupyterhub")
from token_authenticator import TokenAuthenticator


# Bind JupyterHub to the external port requested by the roadmap/deploy plan.
c.JupyterHub.bind_url = "http://0.0.0.0:8888"
c.JupyterHub.hub_connect_ip = "amprenta-jupyterhub"

# Persist hub state/database
c.JupyterHub.db_url = "sqlite:////srv/jupyterhub/jupyterhub.sqlite"

# Authentication: TokenAuthenticator (JWT from Streamlit)
c.JupyterHub.authenticator_class = TokenAuthenticator

# JWT config
_jwt_secret = os.environ.get("JWT_SECRET_KEY")
if not _jwt_secret:
    raise RuntimeError("JWT_SECRET_KEY environment variable is required for JupyterHub")
c.TokenAuthenticator.jwt_secret = _jwt_secret

# Spawner: DockerSpawner (production-style user isolation)
c.JupyterHub.spawner_class = "dockerspawner.DockerSpawner"
c.DockerSpawner.image = "amprenta-singleuser:latest"
c.DockerSpawner.network_name = "jupyterhub-network"
c.DockerSpawner.remove = False  # Remove containers when stopped (disabled temporarily for debugging)
c.DockerSpawner.volumes = {
    "jupyterhub-user-{username}": "/home/jovyan/work",
}
c.DockerSpawner.environment = {
    "API_URL": os.environ.get("API_URL", "http://host.docker.internal:8000"),
    "JUPYTERHUB_SINGLEUSER_APP": "jupyter_server.serverapp.ServerApp",
    # Y.js WebSocket server URL for real-time collaborative editing
    "JUPYTER_COLLABORATION_YJS_URL": os.environ.get("YJS_SERVER_URL", "ws://yjs-server:1234"),
}
c.DockerSpawner.cmd = None

# Default landing page for user servers
c.Spawner.default_url = "/lab"


