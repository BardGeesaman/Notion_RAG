"""Connectivity scoring utilities."""

from __future__ import annotations

import math
from typing import Dict


def compute_connectivity_score(query_vec: Dict[int, float], lincs_vec: Dict[int, float]) -> float:
    """Compute weighted cosine similarity between two vectors keyed by Entrez IDs.

    Weighting scheme:
    - Use per-gene weight = abs(query_vec[g]) (emphasizes strongly-regulated query genes).

    Returns a score in [-1, 1]. If vectors have no overlap or are degenerate, returns 0.0.
    """
    if not query_vec or not lincs_vec:
        return 0.0

    keys = set(query_vec.keys()) & set(lincs_vec.keys())
    if not keys:
        return 0.0

    num = 0.0
    qn = 0.0
    ln = 0.0
    for k in keys:
        q = float(query_vec.get(k, 0.0))
        lv = float(lincs_vec.get(k, 0.0))
        w = abs(q)
        num += w * q * lv
        qn += w * q * q
        ln += w * lv * lv

    if qn <= 0.0 or ln <= 0.0:
        return 0.0

    score = num / (math.sqrt(qn) * math.sqrt(ln))
    # clamp for numeric safety
    if score > 1.0:
        return 1.0
    if score < -1.0:
        return -1.0
    return float(score)


__all__ = ["compute_connectivity_score"]


