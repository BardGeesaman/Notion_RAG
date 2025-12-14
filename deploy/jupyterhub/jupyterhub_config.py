"""
JupyterHub configuration for a simple single-host deployment.

Notes:
- Uses DummyAuthenticator (dev only).
- Persists hub DB and user home directories via docker-compose volumes.
"""

import os
from pathlib import Path


# Bind JupyterHub to the external port requested by the roadmap/deploy plan.
c.JupyterHub.bind_url = "http://0.0.0.0:8888"

# Persist hub state/database
c.JupyterHub.db_url = "sqlite:////srv/jupyterhub/jupyterhub.sqlite"

# Authentication: DummyAuthenticator (dev)
c.JupyterHub.authenticator_class = "dummy"
c.DummyAuthenticator.password = os.environ.get("JUPYTERHUB_PASSWORD", "dev")

# Spawner: local single-server per user
c.JupyterHub.spawner_class = "jupyterhub.spawner.SimpleLocalProcessSpawner"
c.Spawner.default_url = "/lab"

# Store user notebooks in a persistent volume-backed path (no system users required).
c.Spawner.notebook_dir = "/srv/jupyterhub/users/{username}"


def _ensure_user_dir(spawner):
    username = getattr(spawner.user, "name", "user")
    user_dir = Path("/srv/jupyterhub/users") / username
    user_dir.mkdir(parents=True, exist_ok=True)


c.Spawner.pre_spawn_hook = _ensure_user_dir

# Ensure user servers can see API_URL inside notebooks if needed
c.Spawner.environment = {
    "API_URL": os.environ.get("API_URL", "http://localhost:8000"),
}


