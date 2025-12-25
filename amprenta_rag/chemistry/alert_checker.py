"""Unified structural alert checker service (PAINS/Brenk/Lilly)."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Optional

from amprenta_rag.chemistry.alerts import (
    AlertResult,
    BRENKFilter,
    LillyFilter,
    PAINSFilter,
    RDKIT_AVAILABLE,
)


def get_traffic_light(alerts: List[AlertResult]) -> str:
    """Compute traffic light from alerts.

    - RED: any PAINS match (high severity)
    - YELLOW: Brenk/Lilly only
    - GREEN: no alerts
    """
    if any(a.alert_type.upper() == "PAINS" for a in alerts):
        return "RED"
    if alerts:
        return "YELLOW"
    return "GREEN"


def list_available_filters() -> List[Dict[str, Any]]:
    """List supported filters and (best-effort) pattern counts."""
    pains = PAINSFilter()
    brenk = BRENKFilter()
    lilly = LillyFilter()
    return [
        {
            "name": "pains",
            "pattern_count": getattr(pains, "pattern_count", 0),
            "description": "RDKit PAINS FilterCatalog (PAINS_A/B/C)",
        },
        {
            "name": "brenk",
            "pattern_count": getattr(brenk, "pattern_count", 0),
            "description": "RDKit Brenk unwanted substructures FilterCatalog",
        },
        {
            "name": "lilly",
            "pattern_count": getattr(lilly, "pattern_count", 0),
            "description": "Curated Lilly-style SMARTS liabilities (MVP subset)",
        },
    ]


class StructuralAlertChecker:
    """Run selected structural alert filters on compounds."""

    def __init__(self, filters: Optional[List[str]] = None):
        wanted = [s.strip().lower() for s in (filters or ["pains", "brenk", "lilly"])]
        wanted_set = set(wanted)
        self.filters: Dict[str, Any] = {}

        if "pains" in wanted_set:
            self.filters["pains"] = PAINSFilter()
        if "brenk" in wanted_set:
            self.filters["brenk"] = BRENKFilter()
        if "lilly" in wanted_set:
            self.filters["lilly"] = LillyFilter()

    def check_compound(self, smiles: str) -> Dict[str, Any]:
        smi = str(smiles or "").strip()
        if not smi:
            return {
                "smiles": smi,
                "is_clean": False,
                "traffic_light": "RED",
                "alert_count": 0,
                "alerts": [],
                "summary": {"pains_count": 0, "brenk_count": 0, "lilly_count": 0},
                "error": "Empty SMILES",
            }

        # Validate SMILES when RDKit is present; otherwise allow graceful empty.
        if RDKIT_AVAILABLE:
            from rdkit import Chem  # type: ignore

            if Chem.MolFromSmiles(smi) is None:
                return {
                    "smiles": smi,
                    "is_clean": False,
                    "traffic_light": "RED",
                    "alert_count": 0,
                    "alerts": [],
                    "summary": {"pains_count": 0, "brenk_count": 0, "lilly_count": 0},
                    "error": "Invalid SMILES",
                }

        alerts: List[AlertResult] = []
        pains_count = 0
        brenk_count = 0
        lilly_count = 0

        if "pains" in self.filters:
            res = self.filters["pains"].check(smi)
            pains_count = len(res)
            alerts.extend(res)
        if "brenk" in self.filters:
            res = self.filters["brenk"].check(smi)
            brenk_count = len(res)
            alerts.extend(res)
        if "lilly" in self.filters:
            res = self.filters["lilly"].check(smi)
            lilly_count = len(res)
            alerts.extend(res)

        tl = get_traffic_light(alerts)
        return {
            "smiles": smi,
            "is_clean": len(alerts) == 0,
            "traffic_light": tl,
            "alert_count": len(alerts),
            "alerts": [a.to_dict() for a in alerts],
            "summary": {"pains_count": pains_count, "brenk_count": brenk_count, "lilly_count": lilly_count},
        }

    def check_batch(self, smiles_list: List[str], n_jobs: int = 4) -> List[Dict[str, Any]]:
        n_jobs = max(1, int(n_jobs))

        def _one(smi: str) -> Dict[str, Any]:
            try:
                return self.check_compound(smi)
            except Exception as e:  # noqa: BLE001
                return {"smiles": str(smi or ""), "error": str(e), "alerts": [], "alert_count": 0}

        with ThreadPoolExecutor(max_workers=n_jobs) as ex:
            return list(ex.map(_one, list(smiles_list or [])))


__all__ = ["StructuralAlertChecker", "get_traffic_light", "list_available_filters"]


