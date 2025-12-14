"""
Signatures resource client.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from amprenta_rag.api.schemas import Signature, SignatureCreate, SignatureUpdate
from amprenta_rag.client.base import BaseHTTPClient


class SignaturesClient:
    def __init__(self, http: BaseHTTPClient):
        self._http = http

    def list(
        self,
        *,
        skip: int = 0,
        limit: int = 100,
        name: Optional[str] = None,
        program_id: Optional[UUID] = None,
    ) -> List[Signature]:
        params = {
            "skip": skip,
            "limit": limit,
            "name": name,
            "program_id": str(program_id) if program_id else None,
        }
        payload = self._http._request(
            "GET",
            "/api/v1/signatures/",
            params={k: v for k, v in params.items() if v is not None},
        )
        return [Signature.model_validate(item) for item in (payload or [])]

    def get(self, signature_id: UUID) -> Signature:
        payload = self._http._request("GET", f"/api/v1/signatures/{signature_id}")
        return Signature.model_validate(payload)

    def create(self, signature: SignatureCreate) -> Signature:
        payload = self._http._request(
            "POST",
            "/api/v1/signatures/",
            json=signature.model_dump(mode="json", exclude_none=True),
        )
        return Signature.model_validate(payload)

    def update(self, signature_id: UUID, signature: SignatureUpdate) -> Signature:
        payload = self._http._request(
            "PATCH",
            f"/api/v1/signatures/{signature_id}",
            json=signature.model_dump(mode="json", exclude_none=True),
        )
        return Signature.model_validate(payload)

    def delete(self, signature_id: UUID) -> None:
        self._http._request("DELETE", f"/api/v1/signatures/{signature_id}")


