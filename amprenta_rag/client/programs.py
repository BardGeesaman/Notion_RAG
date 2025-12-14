"""
Programs resource client.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from amprenta_rag.api.schemas import Program, ProgramCreate, ProgramUpdate
from amprenta_rag.client.base import BaseHTTPClient


class ProgramsClient:
    def __init__(self, http: BaseHTTPClient):
        self._http = http

    def list(
        self,
        *,
        skip: int = 0,
        limit: int = 100,
        name: Optional[str] = None,
    ) -> List[Program]:
        params = {
            "skip": skip,
            "limit": limit,
            "name": name,
        }
        payload = self._http._request(
            "GET",
            "/api/v1/programs/",
            params={k: v for k, v in params.items() if v is not None},
        )
        return [Program.model_validate(item) for item in (payload or [])]

    def get(self, program_id: UUID) -> Program:
        payload = self._http._request("GET", f"/api/v1/programs/{program_id}")
        return Program.model_validate(payload)

    def create(self, program: ProgramCreate) -> Program:
        payload = self._http._request(
            "POST",
            "/api/v1/programs/",
            json=program.model_dump(mode="json", exclude_none=True),
        )
        return Program.model_validate(payload)

    def update(self, program_id: UUID, program: ProgramUpdate) -> Program:
        payload = self._http._request(
            "PATCH",
            f"/api/v1/programs/{program_id}",
            json=program.model_dump(mode="json", exclude_none=True),
        )
        return Program.model_validate(payload)

    def delete(self, program_id: UUID) -> None:
        self._http._request("DELETE", f"/api/v1/programs/{program_id}")


