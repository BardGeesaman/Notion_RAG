"""
Experiments resource client.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.api.schemas import Experiment, ExperimentCreate, ExperimentUpdate
from amprenta_rag.client.base import BaseHTTPClient


class ExperimentsClient:
    def __init__(self, http: BaseHTTPClient):
        self._http = http

    def list(
        self,
        *,
        skip: int = 0,
        limit: int = 100,
        name: Optional[str] = None,
        program_id: Optional[UUID] = None,
    ) -> List[Experiment]:
        params = {
            "skip": skip,
            "limit": limit,
            "name": name,
            "program_id": str(program_id) if program_id else None,
        }
        payload = self._http._request(
            "GET",
            "/api/v1/experiments/",
            params={k: v for k, v in params.items() if v is not None},
        )
        return [Experiment.model_validate(item) for item in (payload or [])]

    def get(self, experiment_id: UUID) -> Experiment:
        payload = self._http._request("GET", f"/api/v1/experiments/{experiment_id}")
        return Experiment.model_validate(payload)

    def create(self, experiment: ExperimentCreate) -> Experiment:
        payload = self._http._request(
            "POST",
            "/api/v1/experiments/",
            json=experiment.model_dump(mode="json", exclude_none=True),
        )
        return Experiment.model_validate(payload)

    def update(self, experiment_id: UUID, experiment: ExperimentUpdate) -> Experiment:
        payload = self._http._request(
            "PATCH",
            f"/api/v1/experiments/{experiment_id}",
            json=experiment.model_dump(mode="json", exclude_none=True),
        )
        return Experiment.model_validate(payload)

    def delete(self, experiment_id: UUID) -> None:
        self._http._request("DELETE", f"/api/v1/experiments/{experiment_id}")

    def annotate(self, experiment_id: UUID, text: str, annotation_type: Optional[str] = None) -> Dict[str, Any]:
        payload = self._http._request(
            "POST",
            f"/api/v1/experiments/{experiment_id}/annotations",
            json={"text": text, "annotation_type": annotation_type},
        )
        return payload or {}


