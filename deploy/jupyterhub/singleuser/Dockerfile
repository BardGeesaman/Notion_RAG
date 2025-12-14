FROM jupyter/minimal-notebook:latest

USER root

COPY amprenta_rag/__init__.py /opt/amprenta/amprenta_rag/__init__.py

COPY amprenta_rag/models/domain.py /opt/amprenta/amprenta_rag/models/domain.py
RUN mkdir -p /opt/amprenta/amprenta_rag/models && \
    echo "" > /opt/amprenta/amprenta_rag/models/__init__.py

COPY amprenta_rag/api/schemas.py /opt/amprenta/amprenta_rag/api/schemas.py
RUN mkdir -p /opt/amprenta/amprenta_rag/api && \
    echo "" > /opt/amprenta/amprenta_rag/api/__init__.py

COPY amprenta_rag/client/__init__.py /opt/amprenta/amprenta_rag/client/__init__.py
COPY amprenta_rag/client/base.py /opt/amprenta/amprenta_rag/client/base.py
COPY amprenta_rag/client/compounds.py /opt/amprenta/amprenta_rag/client/compounds.py
COPY amprenta_rag/client/datasets.py /opt/amprenta/amprenta_rag/client/datasets.py
COPY amprenta_rag/client/experiments.py /opt/amprenta/amprenta_rag/client/experiments.py
COPY amprenta_rag/client/features.py /opt/amprenta/amprenta_rag/client/features.py
COPY amprenta_rag/client/programs.py /opt/amprenta/amprenta_rag/client/programs.py
COPY amprenta_rag/client/signatures.py /opt/amprenta/amprenta_rag/client/signatures.py
COPY amprenta_rag/client/helpers.py /opt/amprenta/amprenta_rag/client/helpers.py

COPY deploy/jupyterhub/templates/ /home/jovyan/templates/
RUN chown -R ${NB_UID}:${NB_GID} /home/jovyan/templates

ENV PYTHONPATH="/opt/amprenta:${PYTHONPATH}"

USER ${NB_UID}

RUN pip install --no-cache-dir \
    "numpy<2.0" \
    pydantic \
    httpx \
    pandas \
    seaborn \
    matplotlib \
    scipy \
    scikit-learn \
    papermill \
    nbconvert \
    rdkit-pypi

WORKDIR /home/jovyan


