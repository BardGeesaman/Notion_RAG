"""
Features resource client.
"""

from __future__ import annotations

from typing import List, Optional, Union
from uuid import UUID

from amprenta_rag.api.schemas import Feature, FeatureCreate, FeatureUpdate
from amprenta_rag.client.base import BaseHTTPClient
from amprenta_rag.models.domain import FeatureType


class FeaturesClient:
    def __init__(self, http: BaseHTTPClient):
        self._http = http

    def list(
        self,
        *,
        skip: int = 0,
        limit: int = 100,
        name: Optional[str] = None,
        feature_type: Optional[Union[FeatureType, str]] = None,
        dataset_id: Optional[UUID] = None,
    ) -> List[Feature]:
        params = {
            "skip": skip,
            "limit": limit,
            "name": name,
            "feature_type": (feature_type.value if isinstance(feature_type, FeatureType) else feature_type),
            "dataset_id": str(dataset_id) if dataset_id else None,
        }
        payload = self._http._request(
            "GET",
            "/api/v1/features/",
            params={k: v for k, v in params.items() if v is not None},
        )
        return [Feature.model_validate(item) for item in (payload or [])]

    def get(self, feature_id: UUID) -> Feature:
        payload = self._http._request("GET", f"/api/v1/features/{feature_id}")
        return Feature.model_validate(payload)

    def create(self, feature: FeatureCreate) -> Feature:
        payload = self._http._request(
            "POST",
            "/api/v1/features/",
            json=feature.model_dump(mode="json", exclude_none=True),
        )
        return Feature.model_validate(payload)

    def update(self, feature_id: UUID, feature: FeatureUpdate) -> Feature:
        payload = self._http._request(
            "PATCH",
            f"/api/v1/features/{feature_id}",
            json=feature.model_dump(mode="json", exclude_none=True),
        )
        return Feature.model_validate(payload)

    def delete(self, feature_id: UUID) -> None:
        self._http._request("DELETE", f"/api/v1/features/{feature_id}")


