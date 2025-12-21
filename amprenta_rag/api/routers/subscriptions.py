"""
Subscription API routes.

These endpoints manage user subscriptions (create, read, update, delete) and
support subscription checks (e.g., polling for updates). All routes are scoped
to the optional authenticated user when available.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Depends

from amprenta_rag.api.dependencies import get_optional_user_id

from amprenta_rag.api.schemas import (
    SubscriptionCreate,
    SubscriptionUpdate,
    SubscriptionResponse,
)
from amprenta_rag.api.services import subscription_service

router = APIRouter()


@router.get("/subscriptions", response_model=List[SubscriptionResponse])
def list_subscriptions(user_id: Optional[UUID] = Depends(get_optional_user_id)) -> List[SubscriptionResponse]:
    """List subscriptions visible to the current user (or all if unauthenticated)."""
    subs = subscription_service.list_subscriptions(user_id)
    return [SubscriptionResponse.model_validate(s) for s in subs]


@router.post("/subscriptions", response_model=SubscriptionResponse, status_code=201)
def create_subscription(payload: SubscriptionCreate, user_id: Optional[UUID] = Depends(get_optional_user_id)) -> SubscriptionResponse:
    """Create a new subscription for the current user."""
    sub = subscription_service.create_subscription(payload, user_id)
    return SubscriptionResponse.model_validate(sub)


@router.get("/subscriptions/{sub_id}", response_model=SubscriptionResponse)
def get_subscription(sub_id: UUID, user_id: Optional[UUID] = Depends(get_optional_user_id)) -> SubscriptionResponse:
    """Fetch a single subscription by ID if visible to the current user."""
    sub = subscription_service.get_subscription(sub_id, user_id)
    if not sub:
        raise HTTPException(status_code=404, detail="Subscription not found")
    return SubscriptionResponse.model_validate(sub)


@router.put("/subscriptions/{sub_id}", response_model=SubscriptionResponse)
def update_subscription(sub_id: UUID, payload: SubscriptionUpdate, user_id: Optional[UUID] = Depends(get_optional_user_id)) -> SubscriptionResponse:
    """Update an existing subscription owned by/visible to the current user."""
    sub = subscription_service.update_subscription(sub_id, payload, user_id)
    if not sub:
        raise HTTPException(status_code=404, detail="Subscription not found")
    return SubscriptionResponse.model_validate(sub)


@router.delete("/subscriptions/{sub_id}", status_code=204)
def delete_subscription(sub_id: UUID, user_id: Optional[UUID] = Depends(get_optional_user_id)) -> None:
    """Delete a subscription owned by/visible to the current user."""
    deleted = subscription_service.delete_subscription(sub_id, user_id)
    if not deleted:
        raise HTTPException(status_code=404, detail="Subscription not found")
    return None


@router.post("/subscriptions/{sub_id}/check", response_model=dict)
def check_subscription(sub_id: UUID, user_id: Optional[UUID] = Depends(get_optional_user_id)) -> dict:
    """Run a lightweight check for a subscription (e.g., new alerts/changes)."""
    result = subscription_service.check_subscription(sub_id, user_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Subscription not found")
    return result

