## Voila + SAR Delta Explorer — Quickstart (portable across workstations)

This playbook is the “don’t think, just run” path to:
- get the SAR API endpoints online
- seed deterministic SAR data (including a guaranteed cliff)
- render `sar_delta_explorer.ipynb` via Voila

### Prereqs
- **Docker Desktop** running
- Project checked out locally
- Your FastAPI environment available (commonly `conda env` named `myenv`)

### 1) Start the API (FastAPI)

From repo root:

```bash
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate myenv
python -m uvicorn amprenta_rag.api.main:app --reload --port 8000
```

Sanity check:

```bash
curl -s http://localhost:8000/health
curl -s http://localhost:8000/api/v1/sar/targets
```

### 2) Seed SAR demo data (idempotent)

```bash
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate myenv
python scripts/seed_sar_data.py
```

If you want a clean slate:

```bash
python scripts/seed_sar_data.py --reset
```

### 3) Build the singleuser image (bakes notebooks + pinned Voila)

```bash
cd deploy/jupyterhub
docker compose build singleuser
```

Notes:
- The singleuser image pins `voila==0.4.6` + `nbclient==0.7.4` for reliable widget rendering.

### 4) Run a standalone Voila smoke test (recommended)

From repo root:

```bash
docker rm -f voila-sar-test 2>/dev/null || true
docker run -d --name voila-sar-test -p 8866:8866 amprenta-singleuser:latest \
  bash -lc "voila --strip_sources=True --port=8866 --no-browser --Voila.ip=0.0.0.0 /home/jovyan/templates/sar_delta_explorer.ipynb"
```

Open:
- `http://localhost:8866/`

Important:
- From inside the Voila container, the host API should be reachable at:
  - `http://host.docker.internal:8000` (macOS/Windows Docker Desktop)
- If you’re on Linux, you may need host-gateway wiring (example):

```bash
docker run --add-host=host.docker.internal:host-gateway ...
```

### 5) API smoke checks (what “good” looks like)

```bash
curl -s http://localhost:8000/api/v1/sar/targets
curl -s http://localhost:8000/api/v1/sar/targets/CDK2/compounds | head
curl -s "http://localhost:8000/api/v1/sar/targets/CDK2/cliffs?similarity_threshold=0.6&fold_change=10&limit=3"
```

### Common failure modes (fast diagnosis)
- **Voila shows blank / widgets don’t render**:
  - Check browser console for JupyterLab shared-module errors.
  - In this repo, pinning `voila==0.4.6` avoids the federated-extension mismatch.
- **Cliffs endpoint returns `[]` even with low thresholds**:
  - It’s usually **RDKit parsing**: if `MolFromSmiles()` fails for all SMILES, you get zero fingerprints.
  - The API uses best-effort parsing (sanitize=False + skip kekulize) for seeded aromatics.

### Debugging discipline (prevents combinatorial rabbit holes)
- Follow `docs/DEBUG_PROTOCOL.md`:
  - change **one axis at a time**
  - run the smallest deterministic smoke test after each change


