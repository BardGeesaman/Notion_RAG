"""
JupyterHub configuration for a simple single-host deployment.

Notes:
- Uses PAM (Local) authentication.
- For dev convenience, allows any username and will create system users on first login.
- Persists hub DB and user home directories via docker-compose volumes.
"""

import os


# Bind JupyterHub to the external port requested by the roadmap/deploy plan.
c.JupyterHub.bind_url = "http://0.0.0.0:8888"

# Persist hub state/database
c.JupyterHub.db_url = "sqlite:////srv/jupyterhub/jupyterhub.sqlite"

# Authentication: Local (PAM) auth
c.JupyterHub.authenticator_class = "jupyterhub.auth.PAMAuthenticator"

# Dev mode: allow any user and auto-create system users
c.Authenticator.allow_all = True
c.PAMAuthenticator.create_system_users = True

# Spawner: local single-server per user
c.JupyterHub.spawner_class = "jupyterhub.spawner.SimpleLocalProcessSpawner"
c.Spawner.default_url = "/lab"
c.Spawner.notebook_dir = "/home/{username}"

# Ensure user servers can see API_URL inside notebooks if needed
c.Spawner.environment = {
    "API_URL": os.environ.get("API_URL", "http://localhost:8000"),
}


