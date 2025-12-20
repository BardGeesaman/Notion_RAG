from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from amprenta_rag.api.schemas import SubscriptionCreate, SubscriptionUpdate
from amprenta_rag.database.models import RepositorySubscription
from amprenta_rag.database.session import db_session


def _to_dict(sub: RepositorySubscription) -> dict:
    return {
        "id": sub.id,
        "name": sub.name,
        "repository_source": sub.repository_source,
        "query_params": sub.query_params,
        "notify_email": sub.notify_email,
        "notify_in_app": sub.notify_in_app,
        "is_active": sub.is_active,
        "user_id": sub.user_id,
        "last_checked": sub.last_checked,
        "created_at": sub.created_at,
        "updated_at": sub.updated_at,
    }


def list_subscriptions(user_id: Optional[UUID]) -> List[dict]:
    with db_session() as db:
        query = db.query(RepositorySubscription).order_by(RepositorySubscription.created_at.desc())
        if user_id:
            query = query.filter(RepositorySubscription.user_id == user_id)
        subs = query.all()
        return [_to_dict(s) for s in subs]


def create_subscription(payload: SubscriptionCreate, user_id: Optional[UUID]) -> dict:
    with db_session() as db:
        sub = RepositorySubscription(
            name=payload.name,
            repository_source=payload.repository_source,
            query_params=payload.query_params,
            notify_email=payload.notify_email,
            notify_in_app=payload.notify_in_app,
            is_active=True,
            user_id=user_id,  # type: ignore[arg-type]
        )
        db.add(sub)
        db.flush()
        result = _to_dict(sub)
        db.commit()
        return result


def get_subscription(sub_id: UUID, user_id: Optional[UUID]) -> Optional[dict]:
    with db_session() as db:
        query = db.query(RepositorySubscription).filter(RepositorySubscription.id == sub_id)
        if user_id:
            query = query.filter(RepositorySubscription.user_id == user_id)
        sub = query.first()
        return _to_dict(sub) if sub else None


def update_subscription(sub_id: UUID, payload: SubscriptionUpdate, user_id: Optional[UUID]) -> Optional[dict]:
    with db_session() as db:
        query = db.query(RepositorySubscription).filter(RepositorySubscription.id == sub_id)
        if user_id:
            query = query.filter(RepositorySubscription.user_id == user_id)
        sub = query.first()
        if not sub:
            return None
        if payload.name is not None:
            sub.name = payload.name
        if payload.notify_email is not None:
            sub.notify_email = payload.notify_email
        if payload.notify_in_app is not None:
            sub.notify_in_app = payload.notify_in_app
        if payload.is_active is not None:
            sub.is_active = payload.is_active
        if payload.query_params is not None:
            sub.query_params = payload.query_params
        db.flush()
        result = _to_dict(sub)
        db.commit()
        return result


def delete_subscription(sub_id: UUID, user_id: Optional[UUID]) -> bool:
    with db_session() as db:
        query = db.query(RepositorySubscription).filter(RepositorySubscription.id == sub_id)
        if user_id:
            query = query.filter(RepositorySubscription.user_id == user_id)
        sub = query.first()
        if not sub:
            return False
        db.delete(sub)
        db.commit()
        return True


def check_subscription(sub_id: UUID, user_id: Optional[UUID]):
    """Placeholder manual check endpoint; returns empty matches for now."""
    sub = get_subscription(sub_id, user_id)
    if not sub:
        return None
    return {"new_matches": []}

