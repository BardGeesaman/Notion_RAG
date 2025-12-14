"""
Screening resource client.
"""

from __future__ import annotations

from typing import List

from amprenta_rag.api.schemas import CampaignResponse, HTSHitResponse
from amprenta_rag.client.base import BaseHTTPClient


class ScreeningClient:
    def __init__(self, http: BaseHTTPClient):
        self._http = http

    def list_campaigns(self) -> List[CampaignResponse]:
        """List all HTS campaigns."""
        payload = self._http._request("GET", "/api/v1/screening/campaigns")
        return [CampaignResponse.model_validate(item) for item in (payload or [])]

    def get_campaign(self, campaign_id: str) -> CampaignResponse:
        """Get a campaign by campaign_id."""
        payload = self._http._request("GET", f"/api/v1/screening/campaigns/{campaign_id}")
        return CampaignResponse.model_validate(payload)

    def get_campaign_hits(self, campaign_id: str) -> List[HTSHitResponse]:
        """Get hit results for a campaign."""
        payload = self._http._request("GET", f"/api/v1/screening/campaigns/{campaign_id}/hits")
        return [HTSHitResponse.model_validate(item) for item in (payload or [])]

