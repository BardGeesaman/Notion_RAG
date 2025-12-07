from __future__ import annotations

from fastapi import APIRouter, HTTPException

from amprenta_rag.api.services import screening as service

router = APIRouter()


@router.get("/campaigns", summary="List HTS campaigns")
def list_campaigns():
    return service.list_campaigns()


@router.get("/campaigns/{campaign_id}", summary="Get HTS campaign by ID")
def get_campaign(campaign_id: str):
    campaign = service.get_campaign(campaign_id)
    if not campaign:
        raise HTTPException(status_code=404, detail="Campaign not found")
    return campaign


@router.get("/campaigns/{campaign_id}/hits", summary="Get HTS hits for campaign")
def get_campaign_hits(campaign_id: str):
    return service.get_campaign_hits(campaign_id)

