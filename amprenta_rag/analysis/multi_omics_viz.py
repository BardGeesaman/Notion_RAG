from __future__ import annotations

from collections import Counter, defaultdict
from typing import Dict, List, Optional, Set, Tuple
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset, Feature, dataset_feature_assoc


OMICS_COLORS: Dict[str, str] = {
    "transcriptomics": "#E74C3C",
    "proteomics": "#3498DB",
    "metabolomics": "#2ECC71",
    "lipidomics": "#9B59B6",
    "other": "#95A5A6",
}

_FEATURE_NODE_COLOR = "#7F8C8D"
_DATASET_NODE_COLOR = "#95A5A6"

# Guardrails for visualization payload sizes.
_MAX_SHARED_FEATURES_ALLUVIAL = 50
_MAX_FEATURES_UPSET_MATRIX = 500


def _normalize_feature_name(name: str) -> str:
    """
    Normalize feature names for cross-dataset matching.

    - lowercase
    - strip whitespace
    - remove common suffixes: "_HUMAN", "_MOUSE"
    """
    s = (name or "").strip().lower()
    for suf in ("_human", "_mouse"):
        if s.endswith(suf):
            s = s[: -len(suf)]
            s = s.strip()
    return s


def _feature_unification_keys(feature: Feature) -> List[str]:
    """
    Keys used to unify features across datasets.

    Phase 1: feature.normalized_name (if present) and feature.name
    Phase 2: external_ids["hgnc_symbol"], external_ids["uniprot_id"] (if present)
    """
    keys: List[str] = []
    base = getattr(feature, "normalized_name", None) or getattr(feature, "name", None) or ""
    if base:
        keys.append(_normalize_feature_name(str(base)))
    nm = getattr(feature, "name", None)
    if nm:
        keys.append(_normalize_feature_name(str(nm)))

    ext = getattr(feature, "external_ids", None)
    if isinstance(ext, dict):
        hgnc = ext.get("hgnc_symbol")
        if hgnc:
            keys.append(_normalize_feature_name(str(hgnc)))
        unp = ext.get("uniprot_id")
        if unp:
            keys.append(_normalize_feature_name(str(unp)))

    # Dedup preserving order
    out: List[str] = []
    seen: Set[str] = set()
    for k in keys:
        if k and k not in seen:
            out.append(k)
            seen.add(k)
    return out


def _build_unified_feature_index(features: List[Feature]) -> Dict[str, List[Feature]]:
    """
    Group Feature rows by a normalized key so features from different datasets
    can be treated as "the same" entity for overlap visualization.
    """
    idx: Dict[str, List[Feature]] = defaultdict(list)
    for f in features or []:
        for k in _feature_unification_keys(f):
            idx[k].append(f)
    return dict(idx)


def _omics_bucket(dataset_omics_type: Optional[str]) -> str:
    t = (dataset_omics_type or "").strip().lower()
    if t in OMICS_COLORS:
        return t
    return "other"


def _dataset_node_id(dataset_id: UUID) -> str:
    return f"dataset:{dataset_id}"


def _omics_node_id(omics_type: str) -> str:
    return f"omics:{omics_type}"


def _feature_node_id(feature_key: str) -> str:
    return f"feature:{feature_key}"


def _feature_label_for_group(features: List[Feature]) -> str:
    """
    Choose a stable, human-friendly label for a unified feature group.
    Prefer HGNC for gene-like features if available, otherwise Feature.name.
    """
    if not features:
        return ""
    # Prefer HGNC if present
    for f in features:
        ext = getattr(f, "external_ids", None)
        if isinstance(ext, dict) and ext.get("hgnc_symbol"):
            return str(ext["hgnc_symbol"])
    # Fall back to shortest name (often cleaner)
    names = [str(getattr(f, "name", "") or "") for f in features if getattr(f, "name", None)]
    if names:
        return sorted(names, key=lambda x: (len(x), x.lower()))[0]
    return str(getattr(features[0], "id", "feature"))


def compute_alluvial_data(dataset_ids: List[UUID], db: Session) -> Dict:
    """
    Build a simple Sankey/alluvial payload describing:
      dataset -> omics_type -> shared feature (top 50)

    Returns:
      {"nodes": [...], "links": [...]}
    """
    if not dataset_ids:
        return {"nodes": [], "links": []}

    datasets: List[Dataset] = (
        db.query(Dataset).filter(Dataset.id.in_(list(dataset_ids))).all()
    )
    ds_by_id: Dict[UUID, Dataset] = {d.id: d for d in datasets if getattr(d, "id", None)}

    # Fetch feature rows per dataset via association table.
    rows: List[Tuple[UUID, Feature]] = (
        db.query(dataset_feature_assoc.c.dataset_id, Feature)
        .join(Feature, Feature.id == dataset_feature_assoc.c.feature_id)
        .filter(dataset_feature_assoc.c.dataset_id.in_(list(ds_by_id.keys())))
        .all()
    )

    dataset_features: Dict[UUID, List[Feature]] = defaultdict(list)
    all_features: List[Feature] = []
    for dsid, feat in rows:
        dataset_features[dsid].append(feat)
        all_features.append(feat)

    unified = _build_unified_feature_index(all_features)

    # Compute which datasets each unified key appears in (dedup per dataset).
    key_to_datasets: Dict[str, Set[UUID]] = defaultdict(set)
    dataset_to_keys: Dict[UUID, Set[str]] = defaultdict(set)
    for dsid, feats in dataset_features.items():
        for f in feats:
            for k in _feature_unification_keys(f):
                dataset_to_keys[dsid].add(k)
                key_to_datasets[k].add(dsid)

    # Shared features only (present in >= 2 datasets), pick top 50 by dataset count then label.
    shared_keys = [k for k, dsids in key_to_datasets.items() if len(dsids) >= 2]
    shared_keys.sort(
        key=lambda k: (
            -len(key_to_datasets.get(k, set())),
            _feature_label_for_group(unified.get(k, [])),
            k,
        )
    )
    shared_keys = shared_keys[:_MAX_SHARED_FEATURES_ALLUVIAL]

    # Nodes: datasets + omics types (present) + features (shared top N)
    nodes: List[Dict] = []
    node_index: Dict[str, int] = {}

    def _add_node(node_id: str, label: str, color: str) -> int:
        if node_id in node_index:
            return node_index[node_id]
        node_index[node_id] = len(nodes)
        nodes.append({"id": node_id, "label": label, "color": color})
        return node_index[node_id]

    # Dataset nodes (stable order by name then UUID)
    for ds in sorted(ds_by_id.values(), key=lambda d: ((d.name or "").lower(), str(d.id))):
        _add_node(_dataset_node_id(ds.id), ds.name or str(ds.id), _DATASET_NODE_COLOR)

    # Omics nodes (only those present)
    omics_present = sorted({_omics_bucket(ds.omics_type) for ds in ds_by_id.values()})
    for ot in omics_present:
        _add_node(_omics_node_id(ot), ot, OMICS_COLORS.get(ot, OMICS_COLORS["other"]))

    # Feature nodes
    for k in shared_keys:
        label = _feature_label_for_group(unified.get(k, [])) or k
        _add_node(_feature_node_id(k), label, _FEATURE_NODE_COLOR)

    # Links:
    #   dataset -> omics_type (count of unified keys in that dataset)
    #   omics_type -> feature (count of datasets of that omics containing the feature)
    links: List[Dict] = []

    # dataset -> omics_type
    for ds in ds_by_id.values():
        dsid = ds.id
        ot = _omics_bucket(ds.omics_type)
        value = len(dataset_to_keys.get(dsid, set()))
        if value <= 0:
            continue
        links.append(
            {
                "source": node_index[_dataset_node_id(dsid)],
                "target": node_index[_omics_node_id(ot)],
                "value": int(value),
            }
        )

    # omics_type -> feature (aggregate across datasets)
    omics_feature_counts: Dict[Tuple[str, str], int] = Counter()
    for ds in ds_by_id.values():
        ot = _omics_bucket(ds.omics_type)
        keys = dataset_to_keys.get(ds.id, set())
        for k in shared_keys:
            if k in keys:
                omics_feature_counts[(ot, k)] += 1

    for (ot, k), cnt in sorted(
        omics_feature_counts.items(),
        key=lambda x: (-x[1], x[0][0], _feature_label_for_group(unified.get(x[0][1], [])), x[0][1]),
    ):
        links.append(
            {
                "source": node_index[_omics_node_id(ot)],
                "target": node_index[_feature_node_id(k)],
                "value": int(cnt),
            }
        )

    return {"nodes": nodes, "links": links}


def compute_upset_data(dataset_ids: List[UUID], db: Session) -> Dict:
    """
    Compute a simple UpSet payload across datasets based on unified features.

    Returns:
      {
        "sets": [{"id","label","omics_type","color"}, ...],
        "intersections": [{"set_ids":[...], "count": int, "key": str}, ...],
        "matrix": [{"key": str, "label": str, "presence": [0/1,...], "dataset_count": int}, ...]
      }
    """
    if not dataset_ids:
        return {"sets": [], "intersections": [], "matrix": []}

    datasets: List[Dataset] = (
        db.query(Dataset).filter(Dataset.id.in_(list(dataset_ids))).all()
    )
    # Stable ordering for columns
    datasets = sorted(datasets, key=lambda d: ((d.name or "").lower(), str(d.id)))
    ds_ids: List[UUID] = [d.id for d in datasets]

    rows: List[Tuple[UUID, Feature]] = (
        db.query(dataset_feature_assoc.c.dataset_id, Feature)
        .join(Feature, Feature.id == dataset_feature_assoc.c.feature_id)
        .filter(dataset_feature_assoc.c.dataset_id.in_(list(ds_ids)))
        .all()
    )

    dataset_to_keys: Dict[UUID, Set[str]] = defaultdict(set)
    all_features: List[Feature] = []
    for dsid, feat in rows:
        all_features.append(feat)
        for k in _feature_unification_keys(feat):
            dataset_to_keys[dsid].add(k)

    unified = _build_unified_feature_index(all_features)

    # Union of all keys
    all_keys: Set[str] = set()
    for s in dataset_to_keys.values():
        all_keys |= set(s)

    # Cap matrix size for safety: keep most frequent (by dataset prevalence), then alphabetical.
    key_prevalence = {k: sum(1 for dsid in ds_ids if k in dataset_to_keys.get(dsid, set())) for k in all_keys}
    keys_sorted = sorted(
        all_keys,
        key=lambda k: (-key_prevalence.get(k, 0), _feature_label_for_group(unified.get(k, [])), k),
    )
    keys_sorted = keys_sorted[:_MAX_FEATURES_UPSET_MATRIX]

    # Sets
    sets_out: List[Dict] = []
    for d in datasets:
        ot = _omics_bucket(d.omics_type)
        sets_out.append(
            {
                "id": str(d.id),
                "name": d.name or str(d.id),  # Changed from 'label' to 'name' to match schema
                "omics_type": ot,
                "color": OMICS_COLORS.get(ot, OMICS_COLORS["other"]),
                "size": len(dataset_to_keys.get(d.id, set())),  # Add size field
            }
        )

    # Matrix + intersection counts
    intersections_counter: Counter[int] = Counter()
    matrix_out: List[Dict] = []

    for k in keys_sorted:
        presence: List[int] = []
        mask = 0
        for i, dsid in enumerate(ds_ids):
            has = 1 if k in dataset_to_keys.get(dsid, set()) else 0
            presence.append(has)
            if has:
                mask |= 1 << i
        dataset_count = sum(presence)
        if dataset_count == 0:
            continue
        intersections_counter[mask] += 1

        label = _feature_label_for_group(unified.get(k, [])) or k
        matrix_out.append(
            {
                "key": k,
                "label": label,
                "presence": presence,
                "dataset_count": int(dataset_count),
            }
        )

    # Intersections output (sorted by size desc, then mask)
    intersections_out: List[Dict] = []
    for mask, cnt in sorted(intersections_counter.items(), key=lambda x: (-x[1], x[0])):
        set_ids = [str(ds_ids[i]) for i in range(len(ds_ids)) if (mask >> i) & 1]
        bitkey = "".join("1" if (mask >> i) & 1 else "0" for i in range(len(ds_ids)))
        intersections_out.append({"set_ids": set_ids, "count": int(cnt), "key": bitkey})

    return {"sets": sets_out, "intersections": intersections_out, "matrix": matrix_out}


__all__ = [
    "OMICS_COLORS",
    "_normalize_feature_name",
    "_build_unified_feature_index",
    "compute_alluvial_data",
    "compute_upset_data",
]


