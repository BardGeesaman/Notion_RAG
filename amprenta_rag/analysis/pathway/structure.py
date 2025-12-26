"""KEGG pathway structure fetch + cache (KGML parsing).

KEGG-only MVP: fetch KGML, parse nodes/edges, normalize coordinates, and cache on disk.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
from xml.etree import ElementTree as ET

import requests

from amprenta_rag.logging_utils import get_logger


logger = get_logger(__name__)

KEGG_BASE_URL = "https://rest.kegg.jp"
CACHE_DIR = Path("data/pathway_cache")
CACHE_TTL_SECONDS = 7 * 24 * 60 * 60  # 7 days
KEGG_REQUEST_DELAY_SECONDS = 0.5

_LAST_KEGG_REQUEST_TS: float = 0.0


@dataclass(frozen=True)
class PathwayNode:
    id: str
    name: str
    type: str
    x: float
    y: float
    kegg_ids: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class PathwayEdge:
    source: str
    target: str
    type: str
    subtype: Optional[str] = None


@dataclass(frozen=True)
class PathwayStructure:
    pathway_id: str
    name: str
    nodes: List[PathwayNode] = field(default_factory=list)
    edges: List[PathwayEdge] = field(default_factory=list)
    organism: str = ""


def _rate_limit_kegg() -> None:
    global _LAST_KEGG_REQUEST_TS
    now = time.monotonic()
    wait = KEGG_REQUEST_DELAY_SECONDS - (now - _LAST_KEGG_REQUEST_TS)
    if wait > 0:
        time.sleep(wait)
    _LAST_KEGG_REQUEST_TS = time.monotonic()


def _node_type(entry_type: str | None) -> Optional[str]:
    t = (entry_type or "").strip().lower()
    if t in {"reaction", "brite"}:
        return None
    if t in {"gene", "ortholog", "enzyme"}:
        return "gene"
    if t == "compound":
        return "compound"
    if t == "map":
        return "sub-pathway"
    if t == "group":
        # KGML "group" is often a complex; treat as gene-like if it has graphics.
        return "gene"
    return t or None


def _parse_kegg_ids(entry_name: str | None) -> List[str]:
    if not entry_name:
        return []
    # KGML "name" can contain space-separated KEGG IDs, e.g. "hsa:7157 hsa:7422"
    return [p.strip() for p in str(entry_name).split() if p.strip()]


def _empty_structure(pathway_id: str, *, name: Optional[str] = None, organism: Optional[str] = None) -> PathwayStructure:
    org = organism or (pathway_id[:3] if len(pathway_id) >= 3 else "")
    return PathwayStructure(pathway_id=pathway_id, name=name or pathway_id, nodes=[], edges=[], organism=org)


def fetch_kegg_pathway_structure(pathway_id: str) -> PathwayStructure:
    """
    Fetch KEGG KGML and parse into a PathwayStructure.

    Rate-limited to 0.5s between requests.
    """
    pid = (pathway_id or "").strip()
    if not pid:
        return _empty_structure("unknown")

    url = f"{KEGG_BASE_URL}/get/{pid}/kgml"
    try:
        _rate_limit_kegg()
        resp = requests.get(url, timeout=20)
        if resp.status_code != 200 or not resp.text:
            logger.warning("[PATHWAY][STRUCTURE] KGML fetch failed for %s: status=%s", pid, resp.status_code)
            return _empty_structure(pid)

        try:
            root = ET.fromstring(resp.text)
        except Exception as e:  # noqa: BLE001
            logger.warning("[PATHWAY][STRUCTURE] KGML parse error for %s: %r", pid, e)
            return _empty_structure(pid)

        title = root.get("title") or pid
        org = root.get("org") or (pid[:3] if len(pid) >= 3 else "")

        nodes: List[PathwayNode] = []
        edges: List[PathwayEdge] = []

        # Parse nodes from <entry>
        for entry in root.findall("entry"):
            eid = entry.get("id")
            if not eid:
                continue
            mapped = _node_type(entry.get("type"))
            if mapped is None:
                continue

            kegg_ids = _parse_kegg_ids(entry.get("name"))

            # Graphics carries layout; skip if missing.
            g = entry.find("graphics")
            if g is None:
                continue
            try:
                x = float(g.get("x") or 0.0)
                y = float(g.get("y") or 0.0)
            except Exception:
                x, y = 0.0, 0.0
            name = g.get("name") or entry.get("name") or eid

            nodes.append(PathwayNode(id=str(eid), name=str(name), type=mapped, x=x, y=y, kegg_ids=kegg_ids))

        # Normalize coordinates
        max_x = max((n.x for n in nodes), default=1.0)
        max_y = max((n.y for n in nodes), default=1.0)
        max_x = max(max_x, 1.0)
        max_y = max(max_y, 1.0)
        nodes = [
            PathwayNode(id=n.id, name=n.name, type=n.type, x=float(n.x / max_x), y=float(n.y / max_y), kegg_ids=list(n.kegg_ids))
            for n in nodes
        ]

        # Parse edges from <relation>
        for rel in root.findall("relation"):
            s = rel.get("entry1")
            t = rel.get("entry2")
            rtype = rel.get("type") or "relation"
            if not s or not t:
                continue
            subtypes = rel.findall("subtype")
            if not subtypes:
                edges.append(PathwayEdge(source=str(s), target=str(t), type=str(rtype), subtype=None))
            else:
                for st in subtypes:
                    st_name = st.get("name") or st.get("value")
                    edges.append(PathwayEdge(source=str(s), target=str(t), type=str(rtype), subtype=str(st_name) if st_name else None))

        # Parse edges from <reaction> (substrate -> product)
        for rxn in root.findall("reaction"):
            rtype = rxn.get("type") or "reaction"
            subtype = rxn.get("name") or rtype
            subs = [s.get("id") for s in rxn.findall("substrate") if s.get("id")]
            prods = [p.get("id") for p in rxn.findall("product") if p.get("id")]
            for s in subs:
                for p in prods:
                    edges.append(PathwayEdge(source=str(s), target=str(p), type="reaction", subtype=str(subtype)))

        return PathwayStructure(pathway_id=pid, name=str(title), nodes=nodes, edges=edges, organism=str(org))

    except Exception as e:  # noqa: BLE001
        logger.warning("[PATHWAY][STRUCTURE] Failed to fetch structure for %s: %r", pid, e)
        return _empty_structure(pid)


def load_cached_structure(pathway_id: str) -> Optional[PathwayStructure]:
    """
    Load cached structure from data/pathway_cache/{pathway_id}.json if present and fresh.

    Cache expires after 7 days based on file mtime.
    """
    pid = (pathway_id or "").strip()
    if not pid:
        return None
    path = CACHE_DIR / f"{pid}.json"
    try:
        if not path.exists():
            return None
        age = time.time() - path.stat().st_mtime
        if age > CACHE_TTL_SECONDS:
            return None
        raw = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            return None
        nodes = [PathwayNode(**n) for n in (raw.get("nodes") or []) if isinstance(n, dict)]
        edges = [PathwayEdge(**e) for e in (raw.get("edges") or []) if isinstance(e, dict)]
        return PathwayStructure(
            pathway_id=str(raw.get("pathway_id") or pid),
            name=str(raw.get("name") or pid),
            nodes=nodes,
            edges=edges,
            organism=str(raw.get("organism") or (pid[:3] if len(pid) >= 3 else "")),
        )
    except Exception:
        return None


def cache_pathway_structure(pathway_id: str, structure: PathwayStructure) -> None:
    """Persist structure to data/pathway_cache/{pathway_id}.json."""
    pid = (pathway_id or "").strip()
    if not pid:
        return
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        path = CACHE_DIR / f"{pid}.json"
        payload: Dict[str, object] = asdict(structure)  # nested dataclasses -> dict
        path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    except Exception as e:  # noqa: BLE001
        logger.debug("[PATHWAY][STRUCTURE] Cache write failed for %s: %r", pid, e)


def get_pathway_structure(pathway_id: str) -> PathwayStructure:
    """Try cache first; otherwise fetch from KEGG and cache."""
    cached = load_cached_structure(pathway_id)
    if cached is not None:
        return cached
    st = fetch_kegg_pathway_structure(pathway_id)
    cache_pathway_structure(pathway_id, st)
    return st


def edge_style(edge: PathwayEdge) -> Dict[str, str]:
    """
    Map edge subtype to visualization style.

    - activation/expression -> green solid
    - inhibition/repression -> red dashed
    - other -> gray solid
    """
    st = (edge.subtype or "").lower()
    if any(k in st for k in ("activation", "expression")):
        return {"style": "solid", "color": "#2E8B57"}
    if any(k in st for k in ("inhibition", "repression")):
        return {"style": "dashed", "color": "#C73E1D"}
    return {"style": "solid", "color": "#9E9E9E"}


__all__ = [
    "PathwayNode",
    "PathwayEdge",
    "PathwayStructure",
    "fetch_kegg_pathway_structure",
    "load_cached_structure",
    "cache_pathway_structure",
    "get_pathway_structure",
    "edge_style",
]


