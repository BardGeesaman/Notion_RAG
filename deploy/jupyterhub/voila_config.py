"""Voila configuration notes for read-only dashboard rendering.

This repository expects Voila to serve notebook dashboards at:
  /voila/render/{notebook_path}

Typical local command (example):
  voila --no-browser --port=8866 --VoilaConfiguration.enable_nbextensions=True

In JupyterHub deployments, Voila is commonly proxied behind JupyterHub or a reverse proxy.
"""

# This file is intentionally minimal; real deployments should configure:
# - base_url / prefix handling
# - authentication / access control
# - allowed notebook paths (templates directory only)


