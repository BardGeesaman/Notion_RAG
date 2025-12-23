from __future__ import annotations

from pathlib import Path

import pytest

from amprenta_rag.crispr.result_parser import parse_gene_summary


def test_parse_gene_summary_basic(tmp_path: Path):
    p = tmp_path / "mageck.gene_summary.txt"
    p.write_text(
        "\n".join(
            [
                "id\tneg|lfc\tpos|lfc\tneg|p-value\tpos|p-value\tfdr",
                "TP53\t-1.2\t0.3\t0.001\t0.5\t0.01",
                "BRCA1\t-0.2\t1.1\t0.2\t0.01\t0.04",
                "",
            ]
        ),
        encoding="utf-8",
    )
    out = parse_gene_summary(str(p))
    assert len(out) == 2
    assert out[0]["gene"] == "TP53"
    assert out[0]["neg_lfc"] == pytest.approx(-1.2)
    assert out[0]["pos_lfc"] == pytest.approx(0.3)
    assert out[0]["neg_p"] == pytest.approx(0.001)
    assert out[0]["pos_p"] == pytest.approx(0.5)
    assert out[0]["fdr"] == pytest.approx(0.01)


