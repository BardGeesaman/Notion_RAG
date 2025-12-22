FROM python:3.10-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install
COPY requirements.txt .
RUN python -m pip install --upgrade pip && \
    python -m pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY amprenta_rag/ ./amprenta_rag/
COPY alembic/ ./alembic/
COPY alembic.ini .

# Expose FastAPI port
EXPOSE 8000

# Run FastAPI
CMD ["uvicorn", "amprenta_rag.api.main:app", "--host", "0.0.0.0", "--port", "8000"]


