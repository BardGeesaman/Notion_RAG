"""
Datasets resource client.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union
from uuid import UUID

from amprenta_rag.api.schemas import Dataset, DatasetCreate, DatasetUpdate
from amprenta_rag.client.base import BaseHTTPClient
from amprenta_rag.models.domain import OmicsType


class DatasetsClient:
    def __init__(self, http: BaseHTTPClient):
        self._http = http

    def list(
        self,
        *,
        skip: int = 0,
        limit: int = 100,
        name: Optional[str] = None,
        omics_type: Optional[Union[OmicsType, str]] = None,
        program_id: Optional[UUID] = None,
        experiment_id: Optional[UUID] = None,
    ) -> List[Dataset]:
        params = {
            "skip": skip,
            "limit": limit,
            "name": name,
            "omics_type": (omics_type.value if isinstance(omics_type, OmicsType) else omics_type),
            "program_id": str(program_id) if program_id else None,
            "experiment_id": str(experiment_id) if experiment_id else None,
        }
        payload = self._http._request("GET", "/api/v1/datasets/", params={k: v for k, v in params.items() if v is not None})
        return [Dataset.model_validate(item) for item in (payload or [])]

    def get(self, dataset_id: UUID) -> Dataset:
        payload = self._http._request("GET", f"/api/v1/datasets/{dataset_id}")
        return Dataset.model_validate(payload)

    def create(self, dataset: DatasetCreate) -> Dataset:
        payload = self._http._request(
            "POST",
            "/api/v1/datasets/",
            json=dataset.model_dump(mode="json", exclude_none=True),
        )
        return Dataset.model_validate(payload)

    def update(self, dataset_id: UUID, dataset: DatasetUpdate) -> Dataset:
        payload = self._http._request(
            "PATCH",
            f"/api/v1/datasets/{dataset_id}",
            json=dataset.model_dump(mode="json", exclude_none=True),
        )
        return Dataset.model_validate(payload)

    def delete(self, dataset_id: UUID) -> None:
        self._http._request("DELETE", f"/api/v1/datasets/{dataset_id}")

    def annotate(self, dataset_id: UUID, text: str, annotation_type: Optional[str] = None) -> Dict[str, Any]:
        payload = self._http._request(
            "POST",
            f"/api/v1/datasets/{dataset_id}/annotations",
            json={"text": text, "annotation_type": annotation_type},
        )
        return payload or {}


