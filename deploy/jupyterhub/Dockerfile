FROM jupyterhub/jupyterhub:latest

# Keep image behavior simple and explicit
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /srv/jupyterhub

# Install Python deps (kernel libs + http client)
COPY deploy/jupyterhub/requirements.txt /tmp/requirements.txt
RUN python -m pip install --no-cache-dir --upgrade pip \
    && python -m pip install --no-cache-dir -r /tmp/requirements.txt

# Install Miniconda (architecture-aware) + RDKit
RUN apt-get update && apt-get install -y wget && \
    ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"; \
    elif [ "$ARCH" = "aarch64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"; \
    else \
        echo "Unsupported architecture: $ARCH" && exit 1; \
    fi && \
    wget -q "$MINICONDA_URL" -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda install -y -c conda-forge --override-channels rdkit && \
    /opt/conda/bin/conda clean -afy

ENV PATH="/opt/conda/bin:$PATH"

# "Install" amprenta-client from local package.
# The client imports other amprenta_rag modules (api.schemas, models.domain),
# so we copy the required Python package and expose it via PYTHONPATH.
COPY amprenta_rag /opt/amprenta_rag/amprenta_rag
ENV PYTHONPATH=/opt/amprenta_rag:${PYTHONPATH}

# JupyterHub configuration
COPY deploy/jupyterhub/jupyterhub_config.py /srv/jupyterhub/jupyterhub_config.py

EXPOSE 8888

CMD ["jupyterhub", "-f", "/srv/jupyterhub/jupyterhub_config.py"]


