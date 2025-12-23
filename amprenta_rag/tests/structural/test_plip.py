from __future__ import annotations

from pathlib import Path

from amprenta_rag.structural.plip_parser import parse_plip_xml
from amprenta_rag.structural.plip_runner import run_plip


def test_run_plip_prefers_docker(monkeypatch, tmp_path: Path):
    complex_pdb = tmp_path / "complex.pdb"
    complex_pdb.write_text("ATOM      1  N   ALA A   1      0.0   0.0   0.0\nEND\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    calls = []

    def _fake_which(name: str):
        return "/usr/bin/docker" if name == "docker" else None

    def _fake_run(cmd, capture_output, text, check):
        calls.append(cmd)

        class R:
            returncode = 0
            stdout = ""
            stderr = ""

        return R()

    import amprenta_rag.structural.plip_runner as mod

    monkeypatch.setattr(mod.shutil, "which", _fake_which)
    monkeypatch.setattr(mod.subprocess, "run", _fake_run)

    ok = run_plip(str(complex_pdb), str(out_dir))
    assert ok is True
    assert calls and calls[0][0] == "docker"


def test_parse_plip_xml_counts(tmp_path: Path):
    xml = tmp_path / "report.xml"
    xml.write_text(
        """<?xml version="1.0"?>
<report>
  <hydrogen_bonds>
    <hydrogen_bond>
      <ligatom>C1</ligatom>
      <reschain>A</reschain>
      <restype>GLY</restype>
      <resnr>12</resnr>
      <dist>2.9</dist>
      <angle>160.0</angle>
    </hydrogen_bond>
  </hydrogen_bonds>
  <hydrophobic_interactions>
    <hydrophobic_interaction>
      <ligatom>C2</ligatom>
      <reschain>B</reschain>
      <restype>LEU</restype>
      <resnr>5</resnr>
      <dist>3.8</dist>
    </hydrophobic_interaction>
  </hydrophobic_interactions>
</report>
""",
        encoding="utf-8",
    )

    res = parse_plip_xml(str(xml))
    assert res.num_hbonds == 1
    assert res.num_hydrophobic == 1
    assert res.total_interactions == 2
    assert any(i.interaction_type == "hbond" for i in res.interactions)


