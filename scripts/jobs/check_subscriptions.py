"""
Background job to check repository subscriptions for new datasets.

Run: python scripts/jobs/check_subscriptions.py
"""

from __future__ import annotations

# --- Robust import for jobs: always importable regardless of directory ---
import sys
import os
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
# -------------------------------------------------------------------------

from datetime import datetime, timezone
from typing import List
import logging
from uuid import UUID

from sqlalchemy import or_

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import RepositorySubscription, Dataset, RepositoryNotification

logger = logging.getLogger(__name__)


def _sub_to_dict(sub: RepositorySubscription) -> dict:
    return {
        "id": str(sub.id),
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


def _dataset_to_dict(ds: Dataset) -> dict:
    return {
        "id": str(ds.id),
        "name": ds.name,
        "data_origin": ds.data_origin,
        "created_at": ds.created_at,
    }


def send_notifications(sub: dict, alerts: List[dict]):
    """Send notifications for new alerts (MVP logs only)."""
    if not alerts:
        return
    logger.info("[%s] %d new alerts created (in-app)", sub.get("name"), len(alerts))

    # Email placeholder
    if sub.get("notify_email"):
        logger.info("[%s] Email notification would be sent to user (TODO SMTP)", sub.get("name"))

    # In-app polling via alerts API/bell already handled
    if sub.get("notify_in_app"):
        logger.info("[%s] In-app notification ready (alerts API)", sub.get("name"))


def fetch_active_subscriptions() -> List[dict]:
    with db_session() as db:
        subs = (
            db.query(RepositorySubscription)
            .filter(RepositorySubscription.is_active.is_(True))
            .order_by(RepositorySubscription.created_at.desc())
            .all()
        )
        return [_sub_to_dict(s) for s in subs]


def find_new_datasets(sub: dict) -> List[dict]:
    with db_session() as db:
        query = db.query(Dataset)
        if sub.get("repository_source") and sub.get("repository_source") != "all":
            query = query.filter(Dataset.data_origin == sub["repository_source"])
        if sub.get("last_checked"):
            query = query.filter(Dataset.created_at > sub["last_checked"])

        # Keyword filtering from query_params
        query_params = sub.get("query_params") or {}
        keywords = query_params.get("keywords")
        if keywords:
            if isinstance(keywords, str):
                kw_list = [k.strip() for k in keywords.split(",") if k.strip()]
            else:
                kw_list = [k.strip() for k in keywords if k.strip()]
            if kw_list:
                conditions = []
                for kw in kw_list:
                    term = f"%{kw}%"
                    conditions.append(Dataset.name.ilike(term))
                    conditions.append(Dataset.description.ilike(term))
                query = query.filter(or_(*conditions))

        results = query.order_by(Dataset.created_at.desc()).all()
        return [_dataset_to_dict(d) for d in results]


def update_last_checked(sub_ids: List[str]):
    now = datetime.now(timezone.utc)
    with db_session() as db:
        db.query(RepositorySubscription).filter(RepositorySubscription.id.in_(sub_ids)).update(
            {RepositorySubscription.last_checked: now}, synchronize_session=False
        )
        db.commit()


def main():
    subs = fetch_active_subscriptions()
    touched: List[str] = []
    for sub in subs:
        try:
            new_datasets = find_new_datasets(sub)
            if new_datasets:
                logger.info(
                    "[Subscription] %s (%s) new matches: %d",
                    sub.get("name"),
                    sub.get("repository_source"),
                    len(new_datasets),
                )
                with db_session() as db:
                    matched_datasets: List[dict] = []
                    for ds in new_datasets:
                        logger.info("  - %s | %s | %s | %s", ds.get("id"), ds.get("name"), ds.get("data_origin"), ds.get("created_at"))
                        alert = Alert(
                            subscription_id=UUID(sub.get("id")),
                            dataset_id=UUID(ds.get("id")),
                        )
                        db.add(alert)
                        matched_datasets.append(ds)
                    db.commit()
                send_notifications(sub, matched_datasets)
            touched.append(sub.get("id"))
        except Exception as exc:  # noqa: BLE001
            logger.error("Error processing subscription %s: %s", sub.get("id"), exc)
    if touched:
        update_last_checked(touched)
    logger.info("Checked %d subscriptions.", len(subs))


if __name__ == "__main__":
    main()

