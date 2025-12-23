from __future__ import annotations

from pathlib import Path

from amprenta_rag.structural.vina_parser import parse_vina_output
from amprenta_rag.structural.vina_runner import run_vina


def test_parse_vina_output_extracts_result(tmp_path: Path):
    out = tmp_path / "out.pdbqt"
    log = tmp_path / "vina.log"
    out.write_text(
        "REMARK VINA RESULT:     -7.6      0.000      0.000\nATOM ...\n",
        encoding="utf-8",
    )
    log.write_text("Vina log\n", encoding="utf-8")

    res = parse_vina_output(str(out), str(log))
    assert res.binding_affinity == -7.6
    assert res.rmsd_lb == 0.0
    assert res.rmsd_ub == 0.0
    assert "REMARK VINA RESULT" in res.pdbqt_content


def test_run_vina_prefers_docker(monkeypatch, tmp_path: Path):
    rec = tmp_path / "rec.pdbqt"
    lig = tmp_path / "lig.pdbqt"
    rec.write_text("RECEPTOR\n", encoding="utf-8")
    lig.write_text("LIGAND\n", encoding="utf-8")
    out_dir = tmp_path / "run"

    calls = []

    def _fake_which(name: str):
        return "/usr/bin/docker" if name == "docker" else None

    def _fake_run(cmd, capture_output, text, check, cwd=None):
        calls.append(cmd)
        # simulate docker writing output
        (out_dir / "out.pdbqt").write_text("REMARK VINA RESULT: -1.0 0 0\n", encoding="utf-8")
        (out_dir / "vina.log").write_text("ok\n", encoding="utf-8")

        class R:
            returncode = 0
            stdout = ""
            stderr = ""

        return R()

    import amprenta_rag.structural.vina_runner as mod

    monkeypatch.setattr(mod.shutil, "which", _fake_which)
    monkeypatch.setattr(mod.subprocess, "run", _fake_run)

    ok = run_vina(
        receptor_pdbqt=str(rec),
        ligand_pdbqt=str(lig),
        center=(0.0, 0.0, 0.0),
        size=(20.0, 20.0, 20.0),
        output_dir=str(out_dir),
    )
    assert ok is True
    assert calls and calls[0][0] == "docker"


