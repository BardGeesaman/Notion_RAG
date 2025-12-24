from __future__ import annotations

import asyncio
import hashlib
import json
from datetime import datetime, timezone
from typing import Any, AsyncIterator, Callable, Dict, Iterable, List, Optional, Set, Tuple
from uuid import UUID

import httpx

from amprenta_rag.database.models import Compound, PubChemBioassay
from amprenta_rag.sync.adapters.base import BaseSyncAdapter


class PubChemAdapter(BaseSyncAdapter):
    source = "pubchem"

    API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    CID_BATCH = 100
    REQUEST_SLEEP_S = 0.2

    def __init__(self, db_session: Callable[[], Any]):
        # db_session is expected to be the context manager from amprenta_rag.database.session
        self._db_session = db_session

    async def fetch_records(self, since: datetime | None) -> AsyncIterator[dict]:
        # since is currently unused: PubChem PUG-REST responses do not reliably provide a per-record updated timestamp.
        _ = since

        cids = self._get_local_cids()
        if not cids:
            return

        for cid_batch in _chunked(cids, self.CID_BATCH):
            # Step 2: assaysummary per CID batch -> AIDs
            cids_str = ",".join(str(c) for c in cid_batch)
            url = f"{self.API_BASE}/compound/cid/{cids_str}/assaysummary/JSON"
            payload = await self._fetch_json(url)
            aids = _extract_aids(payload)

            # Step 3: concise per AID -> results for many CIDs; we filter back to our CID set.
            cid_set = set(cid_batch)
            for aid in sorted(aids):
                url2 = f"{self.API_BASE}/assay/aid/{aid}/concise/JSON"
                concise = await self._fetch_json(url2)
                assay_name = _extract_assay_name(concise)
                for rec in _extract_concise_activity_records(concise):
                    cid = rec.get("pubchem_cid")
                    if cid is None or cid not in cid_set:
                        continue

                    yield {
                        "external_id": f"PUBCHEM_{cid}_{aid}",
                        "pubchem_cid": cid,
                        "aid": int(aid),
                        "activity_outcome": rec.get("activity_outcome"),
                        "activity_score": rec.get("activity_score"),
                        "assay_name": assay_name,
                    }

                await asyncio.sleep(self.REQUEST_SLEEP_S)

            await asyncio.sleep(self.REQUEST_SLEEP_S)

    def compute_checksum(self, record: dict) -> str:
        payload = json.dumps(record, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        return hashlib.md5(payload).hexdigest()  # noqa: S324

    def map_to_entity(self, record: dict, db_session) -> tuple[str, UUID | None]:
        cid = record.get("pubchem_cid")
        aid = record.get("aid")
        if cid is None or aid is None:
            return ("bioassay", None)

        try:
            cid_i = int(cid)
            aid_i = int(aid)
        except Exception:  # noqa: BLE001
            return ("bioassay", None)

        compound: Optional[Compound] = (
            db_session.query(Compound).filter(Compound.pubchem_cid == cid_i).first()
        )
        if compound is None:
            return ("bioassay", None)

        existing: Optional[PubChemBioassay] = (
            db_session.query(PubChemBioassay)
            .filter(
                PubChemBioassay.compound_id == compound.id,
                PubChemBioassay.pubchem_cid == cid_i,
                PubChemBioassay.aid == aid_i,
            )
            .first()
        )

        now = datetime.now(timezone.utc)
        score = record.get("activity_score")
        try:
            score_f = float(score) if score is not None else None
        except Exception:  # noqa: BLE001
            score_f = None

        if existing is None:
            existing = PubChemBioassay(
                compound_id=compound.id,
                pubchem_cid=cid_i,
                aid=aid_i,
                activity_outcome=(record.get("activity_outcome") or None),
                activity_score=score_f,
                assay_name=(record.get("assay_name") or None),
                synced_at=now,
            )
            db_session.add(existing)
            db_session.flush()
        else:
            existing.activity_outcome = (record.get("activity_outcome") or None)
            existing.activity_score = score_f
            existing.assay_name = (record.get("assay_name") or None)
            existing.synced_at = now
            db_session.flush()

        return ("bioassay", existing.id)

    async def _fetch_json(self, url: str) -> dict:
        max_retries = 6
        backoff_s = 0.5

        for attempt in range(max_retries):
            try:
                async with httpx.AsyncClient(timeout=45.0, headers={"Accept": "application/json"}) as client:
                    resp = await client.get(url)

                if resp.status_code in (429, 503):
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

    def _get_local_cids(self) -> List[int]:
        with self._db_session() as db:  # type: ignore[assignment]
            rows = (
                db.query(Compound.pubchem_cid)
                .filter(Compound.pubchem_cid.isnot(None))
                .distinct()
                .all()
            )
        out: List[int] = []
        for (cid,) in rows:
            try:
                if cid is not None:
                    out.append(int(cid))
            except Exception:  # noqa: BLE001
                continue
        return sorted(set(out))


def _chunked(items: List[int], size: int) -> Iterable[List[int]]:
    for i in range(0, len(items), max(1, int(size))):
        yield items[i : i + size]


def _extract_aids(payload: Any) -> Set[int]:
    # PubChem responses vary by endpoint; be robust by scanning common keys.
    return _collect_ints_by_key(payload, {"AID", "aid"})


def _extract_assay_name(payload: Any) -> Optional[str]:
    # Best-effort: find first non-empty string for key "Name"/"name" under assay description objects.
    # If that fails, return None.
    candidates = _collect_strings_by_key(payload, {"Name", "name"})
    for s in candidates:
        s2 = s.strip()
        if len(s2) >= 3:
            return s2
    return None


def _extract_concise_activity_records(payload: Any) -> List[Dict[str, Any]]:
    # Scan for dicts that look like an activity row: CID + outcome/score.
    out: List[Dict[str, Any]] = []

    def visit(obj: Any) -> None:
        if isinstance(obj, dict):
            lower_keys = {str(k).lower() for k in obj.keys()}
            if "cid" in lower_keys and ("activityoutcome" in lower_keys or "activity_outcome" in lower_keys):
                cid = _get_first(obj, ("CID", "cid"))
                outcome = _get_first(obj, ("ActivityOutcome", "activityoutcome", "activity_outcome"))
                score = _get_first(obj, ("ActivityScore", "activityscore", "activity_score", "score"))
                try:
                    cid_i = int(cid) if cid is not None else None
                except Exception:  # noqa: BLE001
                    cid_i = None
                rec: Dict[str, Any] = {
                    "pubchem_cid": cid_i,
                    "activity_outcome": str(outcome) if outcome is not None else None,
                    "activity_score": score,
                }
                out.append(rec)

            for v in obj.values():
                visit(v)
        elif isinstance(obj, list):
            for v in obj:
                visit(v)

    visit(payload)

    # De-dup by (cid, outcome, score)
    seen: Set[Tuple[Any, Any, Any]] = set()
    deduped: List[Dict[str, Any]] = []
    for r in out:
        key = (r.get("pubchem_cid"), r.get("activity_outcome"), r.get("activity_score"))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(r)
    return deduped


def _collect_ints_by_key(obj: Any, keys: Set[str]) -> Set[int]:
    out: Set[int] = set()
    keys_l = {k.lower() for k in keys}

    def visit(x: Any) -> None:
        if isinstance(x, dict):
            for k, v in x.items():
                if str(k).lower() in keys_l:
                    if isinstance(v, list):
                        for vv in v:
                            try:
                                out.add(int(vv))
                            except Exception:  # noqa: BLE001
                                continue
                    else:
                        try:
                            out.add(int(v))
                        except Exception:  # noqa: BLE001
                            pass
                visit(v)
        elif isinstance(x, list):
            for v in x:
                visit(v)

    visit(obj)
    return out


def _collect_strings_by_key(obj: Any, keys: Set[str]) -> List[str]:
    out: List[str] = []
    keys_l = {k.lower() for k in keys}

    def visit(x: Any) -> None:
        if isinstance(x, dict):
            for k, v in x.items():
                if str(k).lower() in keys_l and isinstance(v, str):
                    out.append(v)
                visit(v)
        elif isinstance(x, list):
            for v in x:
                visit(v)

    visit(obj)
    return out


def _get_first(d: Dict[str, Any], keys: Tuple[str, ...]) -> Any:
    for k in keys:
        if k in d:
            return d.get(k)
        lk = k.lower()
        for kk in d.keys():
            if str(kk).lower() == lk:
                return d.get(kk)
    return None


