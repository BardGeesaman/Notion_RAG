from __future__ import annotations

import pytest


pytestmark = pytest.mark.requires_server


def _skip_if_backend_unavailable(resp) -> None:
    if resp.status_code >= 500:
        # Common in dev/CI when DB isn't reachable or migrations aren't applied yet.
        pytest.skip(f"Phenotype API backend unavailable (status={resp.status_code}): {resp.text[:500]}")


def test_get_phenotype_genes_returns_list(fastapi_server: str) -> None:
    """GET /api/phenotypes/{hpo_id}/genes returns a list (skip if DB has no HPO data)."""
    import httpx

    hpo_id = "HP:0000001"

    # Probe via expand-query to decide whether there is any phenotype-gene data loaded.
    with httpx.Client(timeout=10) as client:
        probe = client.post(
            f"{fastapi_server}/api/phenotypes/expand-query",
            json={"query": f"Patient has {hpo_id}"},
        )
        _skip_if_backend_unavailable(probe)
        probe.raise_for_status()
        probe_json = probe.json()
        if int(probe_json.get("gene_count", 0)) == 0:
            pytest.skip("No phenotype-gene associations loaded in this environment.")

        resp = client.get(f"{fastapi_server}/api/phenotypes/{hpo_id}/genes")
        _skip_if_backend_unavailable(resp)
        resp.raise_for_status()
        genes = resp.json()

    assert isinstance(genes, list)


def test_expand_query_extracts_hpo_ids(fastapi_server: str) -> None:
    """POST /api/phenotypes/expand-query extracts HPO IDs from text."""
    import httpx

    with httpx.Client(timeout=10) as client:
        resp = client.post(
            f"{fastapi_server}/api/phenotypes/expand-query",
            json={"query": "Patient has HP:0000002 and HP:0000001"},
        )
        _skip_if_backend_unavailable(resp)
        resp.raise_for_status()
        out = resp.json()

    assert out["hpo_ids"] == ["HP:0000001", "HP:0000002"]
    assert isinstance(out["genes"], list)
    assert isinstance(out["gene_count"], int)


