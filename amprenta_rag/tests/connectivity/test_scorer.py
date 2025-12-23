from __future__ import annotations

from amprenta_rag.connectivity.scorer import compute_connectivity_score


def test_compute_connectivity_score_basic():
    q = {1: 1.0, 2: -1.0}
    # identical direction -> positive
    l1 = {1: 2.0, 2: -2.0}
    s1 = compute_connectivity_score(q, l1)
    assert s1 > 0.9

    # opposite direction -> negative
    l2 = {1: -2.0, 2: 2.0}
    s2 = compute_connectivity_score(q, l2)
    assert s2 < -0.9


def test_compute_connectivity_score_no_overlap():
    assert compute_connectivity_score({1: 1.0}, {2: 1.0}) == 0.0


