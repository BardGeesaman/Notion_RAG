"""
JupyterHub configuration for a simple single-host deployment.

Notes:
- Uses DummyAuthenticator (dev only).
- Persists hub DB and user home directories via docker-compose volumes.
"""

import os


# Bind JupyterHub to the external port requested by the roadmap/deploy plan.
c.JupyterHub.bind_url = "http://0.0.0.0:8888"
c.JupyterHub.hub_connect_ip = "amprenta-jupyterhub"

# Persist hub state/database
c.JupyterHub.db_url = "sqlite:////srv/jupyterhub/jupyterhub.sqlite"

# Authentication: DummyAuthenticator (dev)
c.JupyterHub.authenticator_class = "dummy"
c.DummyAuthenticator.password = os.environ.get("JUPYTERHUB_PASSWORD", "dev")

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
}
c.DockerSpawner.cmd = None

# Default landing page for user servers
c.Spawner.default_url = "/lab"


