from __future__ import annotations

from pathlib import Path

import pytest

from amprenta_rag.structural.prep import prepare_structure


MINI_PDB = """\
ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.560  13.207  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      13.000  14.650  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      12.560  15.650  10.000  1.00 20.00           O
TER
END
"""


def test_prepare_structure_missing_deps_raises(tmp_path: Path, monkeypatch):
    # Force ImportError path without requiring actual deps.
    import builtins

    # Temporarily make imports fail
    orig_import = builtins.__import__

    def _fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name in ("pdbfixer", "openmm", "openmm.app"):
            raise ImportError("missing")
        return orig_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", _fake_import)

    in_p = tmp_path / "in.pdb"
    out_p = tmp_path / "out.pdb"
    in_p.write_text(MINI_PDB, encoding="utf-8")

    with pytest.raises(ImportError, match="pdbfixer.*openmm"):
        prepare_structure(str(in_p), str(out_p), chain_id="A")


def test_prepare_structure_writes_output_if_deps_present(tmp_path: Path):
    pdbfixer = pytest.importorskip("pdbfixer")
    openmm = pytest.importorskip("openmm")
    assert pdbfixer and openmm  # appease linters

    in_p = tmp_path / "in.pdb"
    out_p = tmp_path / "out.pdb"
    in_p.write_text(MINI_PDB, encoding="utf-8")
    prepare_structure(str(in_p), str(out_p), chain_id="A")
    assert out_p.exists()
    assert out_p.read_text(encoding="utf-8")


