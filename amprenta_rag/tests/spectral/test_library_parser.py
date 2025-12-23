from __future__ import annotations

from pathlib import Path

from amprenta_rag.spectral.library_parser import parse_mgf


def test_parse_mgf_basic(tmp_path: Path):
    mgf = tmp_path / "x.mgf"
    mgf.write_text(
        "\n".join(
            [
                "BEGIN IONS",
                "TITLE=PC 34:1",
                "PEPMASS=760.585 1234",
                "CHARGE=1+",
                "100.0 10",
                "150.0 20",
                "END IONS",
                "",
            ]
        ),
        encoding="utf-8",
    )
    out = parse_mgf(str(mgf))
    assert len(out) == 1
    s = out[0]
    assert abs(s["precursor_mz"] - 760.585) < 1e-6
    assert s["name"] == "PC 34:1"
    assert len(s["peaks"]) == 2


