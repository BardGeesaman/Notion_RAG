FROM jupyter/minimal-notebook:latest

USER root

COPY amprenta_rag/__init__.py /opt/amprenta/amprenta_rag/__init__.py
COPY amprenta_rag/logging_utils.py /opt/amprenta/amprenta_rag/logging_utils.py

# Copy full chemistry module tree for RDKit/Voila notebooks (includes __init__.py)
COPY amprenta_rag/chemistry/ /opt/amprenta/amprenta_rag/chemistry/

COPY amprenta_rag/models/domain.py /opt/amprenta/amprenta_rag/models/domain.py
RUN mkdir -p /opt/amprenta/amprenta_rag/models && \
    echo "" > /opt/amprenta/amprenta_rag/models/__init__.py

COPY amprenta_rag/api/schemas.py /opt/amprenta/amprenta_rag/api/schemas.py
RUN mkdir -p /opt/amprenta/amprenta_rag/api && \
    echo "" > /opt/amprenta/amprenta_rag/api/__init__.py

# Copy full client module tree (notebooks use multiple resource clients)
COPY amprenta_rag/client/ /opt/amprenta/amprenta_rag/client/

# Copy pathway analysis utilities for pathway_enrichment notebook (includes __init__.py)
RUN mkdir -p /opt/amprenta/amprenta_rag/analysis
COPY amprenta_rag/analysis/__init__.py /opt/amprenta/amprenta_rag/analysis/__init__.py
COPY amprenta_rag/analysis/pathway/ /opt/amprenta/amprenta_rag/analysis/pathway/

COPY deploy/jupyterhub/templates/ /home/jovyan/templates/
RUN chown -R ${NB_UID}:${NB_GID} /home/jovyan/templates

ENV PYTHONPATH="/opt/amprenta:${PYTHONPATH}"

# Install JupyterLab/Notebook via conda/mamba to keep Python packages and served static assets aligned.
# Mixing pip-installed JupyterLab/Notebook on top of the base image can leave stale CORE_OUTPUT assets
# in place and trigger federated-extension shared-module version conflicts in the browser.
RUN mamba install -y -c conda-forge \
    "jupyterlab==4.5.*" \
    "notebook==7.5.*" \
    && mamba clean -afy

# Rebuild JupyterLab app assets to ensure the served CORE_OUTPUT bundle matches the installed
# JupyterLab/Notebook versions and federated extensions.
RUN jupyter lab clean && jupyter lab build --dev-build=False --minimize=False

USER ${NB_UID}

RUN pip install --no-cache-dir \
    "numpy<2.0" \
    pydantic \
    httpx \
    pandas \
    seaborn \
    matplotlib \
    rdkit-pypi \
    # Voila 0.5.x uses a JupyterLab/federated-extension frontend path that has repeatedly
    # hit shared-module version mismatches in this repo's image stack. Pinning to 0.4.6
    # yields reliable ipywidgets rendering for our dashboards (SAR Delta Explorer, HTS, etc.).
    "voila==0.4.6" \
    "nbclient==0.7.4" \
    ipywidgets

# Enable voila as a Jupyter Server extension so /voila/render works under JupyterHub.
RUN mkdir -p /etc/jupyter/jupyter_server_config.d && \
    printf '%s\n' \
      '{' \
      '  "ServerApp": {' \
      '    "jpserver_extensions": {' \
      '      "voila": true,' \
      '      "voila.server_extension": true' \
      '    }' \
      '  },' \
      '  "Voila": {' \
      '    "root_dir": "/home/jovyan"' \
      '  },' \
      '  "VoilaConfiguration": {' \
      '    "show_tracebacks": true,' \
      '    "template": "classic",' \
      '    "file_whitelist": [".*"]' \
      '  }' \
      '}' \
    > /etc/jupyter/jupyter_server_config.d/voila.json

WORKDIR /home/jovyan


