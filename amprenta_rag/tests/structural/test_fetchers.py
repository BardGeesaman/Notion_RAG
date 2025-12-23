from __future__ import annotations

from pathlib import Path
from uuid import uuid4

from amprenta_rag.structural.alphafold_fetcher import fetch_from_alphafold
from amprenta_rag.structural.pdb_fetcher import fetch_from_pdb
from amprenta_rag.structural.pdb_parser import parse_pdb_metadata
from amprenta_rag.structural.storage import save_structure_file


class _Resp:
    def __init__(self, content: bytes, status_code: int = 200):
        self.content = content
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")


def test_fetch_from_pdb_makes_expected_url(monkeypatch):
    seen = {}

    def _fake_get(url: str, timeout: int):
        seen["url"] = url
        seen["timeout"] = timeout
        return _Resp(b"PDBDATA")

    import amprenta_rag.structural.pdb_fetcher as mod

    monkeypatch.setattr(mod.requests, "get", _fake_get)
    out = fetch_from_pdb("1ABC")
    assert out == b"PDBDATA"
    assert seen["url"] == "https://files.rcsb.org/download/1ABC.pdb"


def test_fetch_from_alphafold_makes_expected_url(monkeypatch):
    seen = {}

    def _fake_get(url: str, timeout: int):
        seen["url"] = url
        return _Resp(b"AFDATA")

    import amprenta_rag.structural.alphafold_fetcher as mod

    monkeypatch.setattr(mod.requests, "get", _fake_get)
    out = fetch_from_alphafold("P12345")
    assert out == b"AFDATA"
    assert seen["url"] == "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.pdb"


def test_save_structure_file_writes_to_expected_path(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("STRUCTURE_STORE_BASE", str(tmp_path))
    sid = uuid4()
    content = b"HELLO"
    path = save_structure_file(sid, content, "pdb")
    p = Path(path)
    assert p.exists()
    assert p.read_bytes() == content
    assert p.parent.name == str(sid)


def test_parse_pdb_metadata_extracts_method_resolution_and_seqres():
    pdb = (
        "EXPDTA    X-RAY DIFFRACTION\n"
        "REMARK   2 RESOLUTION.    2.00 ANGSTROMS.\n"
        "SEQRES   1 A    3  MET GLY ALA\n"
        "END\n"
    ).encode("utf-8")
    out = parse_pdb_metadata(pdb)
    assert out["method"] == "X-RAY DIFFRACTION"
    assert out["resolution"] == 2.0
    assert "A" in (out["sequences"] or {})


