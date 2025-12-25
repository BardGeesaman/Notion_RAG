"""AutoML notebook templates API (list + launch + run tracking)."""

from __future__ import annotations

import os
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field


router = APIRouter(prefix="/automl", tags=["AutoML"])


TEMPLATES_DIR = (Path(__file__).resolve().parents[3] / "deploy" / "jupyterhub" / "templates" / "automl").resolve()
JUPYTERHUB_URL = os.environ.get("JUPYTERHUB_URL", os.environ.get("API_URL", "http://localhost:8000")).rstrip("/")

# In-process run store (MVP). In production, persist to DB.
_RUNS: Dict[str, Dict[str, Any]] = {}


def _list_templates() -> List[Dict[str, Any]]:
    if not TEMPLATES_DIR.exists():
        return []
    out: List[Dict[str, Any]] = []
    for p in sorted(TEMPLATES_DIR.glob("*.ipynb")):
        rel = p.relative_to(TEMPLATES_DIR.parents[1])  # .../templates/<rel>
        nb_rel = str(rel).replace("\\", "/")
        tpl_id = p.stem
        out.append(
            {
                "id": tpl_id,
                "title": tpl_id.replace("_", " ").title(),
                "path": nb_rel,
                "jupyter_url": f"{JUPYTERHUB_URL}/notebooks/{nb_rel}",
            }
        )
    return out


class AutoMLLaunchRequest(BaseModel):
    template_id: str
    params: Dict[str, Any] = Field(default_factory=dict)
    run_mode: str = Field("jupyter", description="jupyter|papermill")


@router.get("/templates", response_model=List[Dict[str, Any]])
def list_templates() -> List[Dict[str, Any]]:
    return _list_templates()


@router.post("/launch", response_model=Dict[str, Any])
def launch_template(payload: AutoMLLaunchRequest) -> Dict[str, Any]:
    templates = {t["id"]: t for t in _list_templates()}
    tpl = templates.get(payload.template_id)
    if not tpl:
        raise HTTPException(status_code=404, detail="Template not found")

    if payload.run_mode == "jupyter":
        return {"mode": "jupyter", "template": tpl, "run_id": None}

    if payload.run_mode != "papermill":
        raise HTTPException(status_code=400, detail="Invalid run_mode")

    try:
        import papermill as pm  # type: ignore
    except Exception:
        raise HTTPException(status_code=501, detail="papermill not installed")

    run_id = str(uuid.uuid4())
    root = Path(__file__).resolve().parents[3]
    out_dir = root / "data" / "automl_runs" / run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_ipynb = out_dir / "output.ipynb"

    nb_path = (TEMPLATES_DIR / f"{payload.template_id}.ipynb").resolve()
    if not nb_path.exists():
        raise HTTPException(status_code=404, detail="Notebook file missing on disk")

    _RUNS[run_id] = {
        "id": run_id,
        "template_id": payload.template_id,
        "status": "running",
        "started_at": time.time(),
        "ended_at": None,
        "params": payload.params,
        "output_ipynb": str(out_ipynb),
        "error": None,
    }

    try:
        pm.execute_notebook(str(nb_path), str(out_ipynb), parameters=payload.params)
        _RUNS[run_id]["status"] = "success"
    except Exception as e:  # noqa: BLE001
        _RUNS[run_id]["status"] = "failed"
        _RUNS[run_id]["error"] = str(e)
    finally:
        _RUNS[run_id]["ended_at"] = time.time()

    return {"mode": "papermill", "template": tpl, "run_id": run_id}


@router.get("/runs", response_model=List[Dict[str, Any]])
def list_runs() -> List[Dict[str, Any]]:
    return list(_RUNS.values())[-200:]


@router.get("/runs/{run_id}", response_model=Dict[str, Any])
def get_run(run_id: str) -> Dict[str, Any]:
    r = _RUNS.get(run_id)
    if not r:
        raise HTTPException(status_code=404, detail="Run not found")
    return r


__all__ = ["router"]


