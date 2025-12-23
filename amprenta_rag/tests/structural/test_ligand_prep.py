from __future__ import annotations

from pathlib import Path

import pytest

from amprenta_rag.structural.ligand_prep import prepare_ligand


def test_prepare_ligand_invalid_smiles(tmp_path: Path):
    out = tmp_path / "x.pdbqt"
    assert prepare_ligand("not_a_smiles", str(out)) is False


def test_prepare_ligand_obabel_fallback(monkeypatch, tmp_path: Path):
    pytest.importorskip("rdkit")
    # Force Meeko path to fail (either not installed or intentionally skipped).
    # Mock OpenBabel CLI presence + subprocess to avoid needing actual obabel.
    out = tmp_path / "lig.pdbqt"

    import amprenta_rag.structural.ligand_prep as mod

    monkeypatch.setattr(mod.shutil, "which", lambda name: "/usr/bin/obabel" if name == "obabel" else None)

    def _fake_run(cmd, capture_output, text, check):
        # Simulate writing the output file
        Path(cmd[-1]).write_text("REMARK PDBQT\n", encoding="utf-8")

        class R:
            returncode = 0
            stdout = ""
            stderr = ""

        return R()

    monkeypatch.setattr(mod.subprocess, "run", _fake_run)

    ok = prepare_ligand("CCO", str(out))
    assert ok is True
    assert out.exists()


def test_prepare_ligand_meeko_if_available(tmp_path: Path):
    pytest.importorskip("rdkit")
    pytest.importorskip("meeko")
    out = tmp_path / "meeko.pdbqt"
    assert prepare_ligand("CCO", str(out)) is True
    assert out.exists()


