"""
Compounds resource client.

Note: the current API router (`amprenta_rag/api/routers/compounds.py`) exposes:
- GET /api/v1/compounds/
- GET /api/v1/compounds/{compound_id}
- GET /api/v1/compounds/{compound_id}/programs

CRUD-shaped methods are provided for interface consistency; create/update/delete
will raise APIError if the server does not implement those routes yet.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from amprenta_rag.api.schemas import CompoundResponse, ProgramLinkResponse
from amprenta_rag.client.base import BaseHTTPClient


class CompoundsClient:
    def __init__(self, http: BaseHTTPClient):
        self._http = http

    def list(self) -> List[CompoundResponse]:
        payload = self._http._request("GET", "/api/v1/compounds/")
        return [CompoundResponse.model_validate(item) for item in (payload or [])]

    def get(self, compound_id: str) -> CompoundResponse:
        payload = self._http._request("GET", f"/api/v1/compounds/{compound_id}")
        return CompoundResponse.model_validate(payload)

    def get_programs(self, compound_id: str) -> List[ProgramLinkResponse]:
        payload = self._http._request("GET", f"/api/v1/compounds/{compound_id}/programs")
        return [ProgramLinkResponse.model_validate(item) for item in (payload or [])]

    def create(self, compound: CompoundResponse) -> CompoundResponse:
        payload = self._http._request(
            "POST",
            "/api/v1/compounds/",
            json=compound.model_dump(mode="json", exclude_none=True),
        )
        return CompoundResponse.model_validate(payload)

    def update(self, compound_id: str, compound: CompoundResponse) -> CompoundResponse:
        payload = self._http._request(
            "PATCH",
            f"/api/v1/compounds/{compound_id}",
            json=compound.model_dump(mode="json", exclude_none=True),
        )
        return CompoundResponse.model_validate(payload)

    def delete(self, compound_id: str) -> None:
        self._http._request("DELETE", f"/api/v1/compounds/{compound_id}")

    def annotate(self, compound_id: str, text: str, annotation_type: Optional[str] = None) -> Dict[str, Any]:
        payload = self._http._request(
            "POST",
            f"/api/v1/compounds/{compound_id}/annotations",
            json={"text": text, "annotation_type": annotation_type},
        )
        return payload or {}


