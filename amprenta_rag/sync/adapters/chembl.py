from __future__ import annotations

import asyncio
import hashlib
import json
from datetime import datetime, timezone
from typing import Any, AsyncIterator, Dict, Optional
from uuid import UUID

import httpx

from amprenta_rag.database.models import ChEMBLActivity, Compound
from amprenta_rag.sync.adapters.base import BaseSyncAdapter


class ChEMBLAdapter(BaseSyncAdapter):
    source = "chembl"

    API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
    PAGE_LIMIT = 1000
    PAGE_SLEEP_S = 0.1
    ACTIVITY_TYPES = {"IC50", "Ki", "Kd", "EC50"}

    async def fetch_records(self, since: datetime | None) -> AsyncIterator[dict]:
        offset = 0

        # Prefer server-side filtering if supported; also enforce client-side filtering.
        standard_type_in = ",".join(sorted(self.ACTIVITY_TYPES))

        while True:
            url = f"{self.API_BASE}/activity.json"
            params = {"limit": self.PAGE_LIMIT, "offset": offset, "standard_type__in": standard_type_in}

            payload = await self._fetch_page(url=url, params=params)
            activities = payload.get("activities") or []

            for a in activities:
                stype = (a.get("standard_type") or "").strip()
                if stype not in self.ACTIVITY_TYPES:
                    continue

                if since is not None and not self._passes_since_filter(a, since):
                    continue

                activity_id = a.get("activity_id") or a.get("id")
                if activity_id is None:
                    # Can't build external_id; skip.
                    continue

                yield {
                    "external_id": f"CHEMBL_ACT_{activity_id}",
                    "molecule_chembl_id": a.get("molecule_chembl_id"),
                    "assay_chembl_id": a.get("assay_chembl_id"),
                    "standard_type": stype,
                    "standard_value": a.get("standard_value"),
                    "standard_units": a.get("standard_units"),
                    "standard_relation": a.get("standard_relation"),
                    "target_chembl_id": a.get("target_chembl_id"),
                    "target_pref_name": a.get("target_pref_name"),
                    "target_type": a.get("target_type"),
                    "assay_type": a.get("assay_type"),
                }

            # Pagination: prefer next link if present, else offset-based until no results.
            next_url = (payload.get("page_meta") or {}).get("next")
            if next_url:
                offset += self.PAGE_LIMIT
            else:
                if len(activities) < self.PAGE_LIMIT:
                    break
                offset += self.PAGE_LIMIT

            await asyncio.sleep(self.PAGE_SLEEP_S)

    def compute_checksum(self, record: dict) -> str:
        payload = json.dumps(record, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        return hashlib.md5(payload).hexdigest()  # noqa: S324

    def map_to_entity(self, record: dict, db_session) -> tuple[str, UUID | None]:
        molecule_chembl_id = (record.get("molecule_chembl_id") or "").strip()
        if not molecule_chembl_id:
            return ("activity", None)

        compound: Optional[Compound] = (
            db_session.query(Compound).filter(Compound.chembl_id == molecule_chembl_id).first()
        )
        if compound is None:
            return ("activity", None)

        chembl_assay_id = (record.get("assay_chembl_id") or "").strip()
        activity_type = (record.get("standard_type") or "").strip() or None
        relation = (record.get("standard_relation") or "").strip() or None
        units = (record.get("standard_units") or "").strip() or None

        # Best-effort de-dup: match exact assay + activity_type + value + relation + units.
        value = record.get("standard_value")
        try:
            value_f = float(value) if value is not None else None
        except Exception:  # noqa: BLE001
            value_f = None

        existing: Optional[ChEMBLActivity] = (
            db_session.query(ChEMBLActivity)
            .filter(
                ChEMBLActivity.compound_id == compound.id,
                ChEMBLActivity.chembl_molecule_id == molecule_chembl_id,
                ChEMBLActivity.chembl_assay_id == chembl_assay_id,
                ChEMBLActivity.activity_type == activity_type,
                ChEMBLActivity.relation == relation,
                ChEMBLActivity.units == units,
                ChEMBLActivity.value == value_f,
            )
            .first()
        )

        now = datetime.now(timezone.utc)
        if existing is None:
            existing = ChEMBLActivity(
                compound_id=compound.id,
                chembl_molecule_id=molecule_chembl_id,
                chembl_assay_id=chembl_assay_id,
                activity_type=activity_type,
                value=value_f,
                units=units,
                relation=relation,
                target_chembl_id=(record.get("target_chembl_id") or "").strip() or None,
                target_name=(record.get("target_pref_name") or "").strip() or None,
                assay_type=(record.get("assay_type") or "").strip() or None,
                synced_at=now,
            )
            db_session.add(existing)
            db_session.flush()
        else:
            existing.target_chembl_id = (record.get("target_chembl_id") or "").strip() or None
            existing.target_name = (record.get("target_pref_name") or "").strip() or None
            existing.assay_type = (record.get("assay_type") or "").strip() or None
            existing.synced_at = now
            db_session.flush()

        return ("activity", existing.id)

    async def _fetch_page(self, url: str, params: Optional[Dict[str, Any]] = None) -> dict:
        max_retries = 6
        backoff_s = 0.5

        for attempt in range(max_retries):
            try:
                async with httpx.AsyncClient(timeout=30.0, headers={"Accept": "application/json"}) as client:
                    resp = await client.get(url, params=params)

                # Rate limit handling
                if resp.status_code == 429:
                    retry_after = resp.headers.get("Retry-After")
                    if retry_after:
                        try:
                            sleep_s = float(retry_after)
                        except Exception:  # noqa: BLE001
                            sleep_s = backoff_s
                    else:
                        sleep_s = backoff_s

                    await asyncio.sleep(sleep_s)
                    backoff_s *= 2
                    continue

                # Retry transient server errors
                if 500 <= resp.status_code < 600:
                    await asyncio.sleep(backoff_s)
                    backoff_s *= 2
                    continue

                resp.raise_for_status()
                return resp.json()
            except httpx.HTTPError:
                if attempt >= max_retries - 1:
                    raise
                await asyncio.sleep(backoff_s)
                backoff_s *= 2

        raise RuntimeError("Unreachable: fetch retries exhausted")

    @staticmethod
    def _passes_since_filter(raw_activity: dict, since: datetime) -> bool:
        """Best-effort client-side incremental filter.

        ChEMBL activities do not reliably expose a single updated timestamp in the JSON payload.
        We attempt a few common keys and fall back to 'include' when none exist.
        """

        keys = ("updated_on", "modified_on", "last_updated", "created_on")
        for k in keys:
            v = raw_activity.get(k)
            if not v:
                continue
            try:
                # Support both date-only and full ISO timestamps.
                dt = datetime.fromisoformat(str(v).replace("Z", "+00:00"))
                if dt.tzinfo is None:
                    dt = dt.replace(tzinfo=timezone.utc)
                return dt >= since
            except Exception:  # noqa: BLE001
                continue
        return True


