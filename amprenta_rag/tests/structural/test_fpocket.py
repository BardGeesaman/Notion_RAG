from __future__ import annotations

from pathlib import Path

from amprenta_rag.structural.fpocket_parser import parse_fpocket_output
from amprenta_rag.structural.fpocket_runner import run_fpocket


def test_run_fpocket_prefers_docker(monkeypatch, tmp_path: Path):
    in_p = tmp_path / "in.pdb"
    in_p.write_text("END\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    calls = []

    def _fake_which(name: str):
        return "/usr/bin/docker" if name == "docker" else None

    def _fake_run(cmd, capture_output, text, check, cwd=None):
        calls.append(cmd)
        class R:
            returncode = 0
            stdout = ""
            stderr = ""
        return R()

    import amprenta_rag.structural.fpocket_runner as mod

    monkeypatch.setattr(mod.shutil, "which", _fake_which)
    monkeypatch.setattr(mod.subprocess, "run", _fake_run)

    ok = run_fpocket(str(in_p), str(out_dir))
    assert ok is True
    assert calls and calls[0][0] == "docker"


def test_parse_fpocket_output_parses_and_sorts(tmp_path: Path):
    # Create minimal info file with two pockets
    info = tmp_path / "test_info.txt"
    info.write_text(
        "Pocket 1\nScore : 10.0\nVolume : 100.0\nPocket 2\nScore : 20.0\nVolume : 50.0\n",
        encoding="utf-8",
    )
    # Create pocket PDB for pocket 2 with REMARK residue hints (optional parser path).
    (tmp_path / "pocket2.pdb").write_text(
        "REMARK RESIDUES A:GLY:12 B:SER:7\nEND\n",
        encoding="utf-8",
    )
    # Create minimal pqr with two TER blocks so we get 2 centers.
    pqr = tmp_path / "test_pockets.pqr"
    pqr.write_text(
        "ATOM 1 X X 0 0 0  1.0 2.0 3.0\nTER\nATOM 1 X X 0 0 0  4.0 5.0 6.0\nTER\n",
        encoding="utf-8",
    )

    pockets = parse_fpocket_output(str(tmp_path))
    assert len(pockets) == 2
    # Sorted by score desc: pocket 2 first
    assert pockets[0].pocket_rank == 2
    assert pockets[0].score == 20.0
    assert pockets[0].residues == ["A:GLY:12", "B:SER:7"]
    assert pockets[1].pocket_rank == 1


def test_parse_fpocket_output_empty(tmp_path: Path):
    assert parse_fpocket_output(str(tmp_path)) == []


