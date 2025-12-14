FROM jupyterhub/jupyterhub:latest

# Keep image behavior simple and explicit
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /srv/jupyterhub

# Install Python deps (kernel libs + http client + rdkit)
COPY deploy/jupyterhub/requirements.txt /tmp/requirements.txt
RUN python -m pip install --no-cache-dir --upgrade pip \
    && python -m pip install --no-cache-dir -r /tmp/requirements.txt

# "Install" amprenta-client from local package.
# The client imports other amprenta_rag modules (api.schemas, models.domain),
# so we copy the required Python package and expose it via PYTHONPATH.
COPY amprenta_rag /opt/amprenta_rag/amprenta_rag
ENV PYTHONPATH=/opt/amprenta_rag:${PYTHONPATH}

# JupyterHub configuration
COPY deploy/jupyterhub/jupyterhub_config.py /srv/jupyterhub/jupyterhub_config.py

EXPOSE 8888

CMD ["jupyterhub", "-f", "/srv/jupyterhub/jupyterhub_config.py"]


