from __future__ import annotations

import json
from pathlib import Path

from amprenta_rag.connectivity.lincs_parser import parse_gct


def test_parse_gct_with_sidecar_meta(tmp_path: Path):
    gct = tmp_path / "tiny.gct"
    gct.write_text(
        "\n".join(
            [
                "#1.2",
                "2\t2",
                "Name\tDescription\tsigA\tsigB",
                "101\tGENE1\t1.0\t-1.5",
                "102\tGENE2\t0.5\t2.0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    meta = tmp_path / "tiny.gct.meta.json"
    meta.write_text(
        json.dumps(
            {
                "sigA": {"pert_iname": "drugA", "pert_id": "BRD-1", "cell_id": "A375"},
                "sigB": {"pert_iname": "drugB", "pert_id": "BRD-2", "cell_id": "PC3"},
            }
        ),
        encoding="utf-8",
    )

    rows = list(parse_gct(str(gct)))
    assert len(rows) == 2
    a = next(r for r in rows if r.sig_id == "sigA")
    b = next(r for r in rows if r.sig_id == "sigB")
    assert a.pert_iname == "drugA"
    assert a.cell_id == "A375"
    assert a.gene_expression == {101: 1.0, 102: 0.5}
    assert b.gene_expression == {101: -1.5, 102: 2.0}


