from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_conformers_endpoint_returns_pdb_strings(monkeypatch):
    from amprenta_rag.api.routers import viz3d as viz3d_router

    monkeypatch.setattr(viz3d_router, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(viz3d_router, "generate_conformers", lambda smiles, n_conformers=5, method="ETKDG": [object()] * int(n_conformers))
    monkeypatch.setattr(viz3d_router, "optimize_conformer", lambda mol, force_field="MMFF": mol)
    monkeypatch.setattr(viz3d_router, "conformer_to_pdb", lambda mol, conf_id=0: "PDBSTR")
    monkeypatch.setattr(viz3d_router, "conformer_energies", lambda mol, prefer="MMFF": [1.23])

    resp = client.post("/api/viz3d/conformers", json={"smiles": "CCO", "n_conformers": 3, "optimize": True})
    assert resp.status_code == 200
    data = resp.json()
    assert data["pdb_strings"] == ["PDBSTR", "PDBSTR", "PDBSTR"]
    assert data["energies"] == [1.23, 1.23, 1.23]


def test_conformers_endpoint_invalid_smiles_returns_400(monkeypatch):
    from amprenta_rag.api.routers import viz3d as viz3d_router

    monkeypatch.setattr(viz3d_router, "RDKIT_AVAILABLE", True)

    def _boom(smiles, n_conformers=5, method="ETKDG"):  # noqa: ANN001
        raise ValueError("Invalid SMILES")

    monkeypatch.setattr(viz3d_router, "generate_conformers", _boom)
    resp = client.post("/api/viz3d/conformers", json={"smiles": "BAD", "n_conformers": 1, "optimize": True})
    assert resp.status_code == 400


def test_overlay_endpoint_returns_aligned_pdbs(monkeypatch):
    from amprenta_rag.api.routers import viz3d as viz3d_router

    monkeypatch.setattr(viz3d_router, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(viz3d_router, "generate_conformers", lambda smiles, n_conformers=1, method="ETKDG": [object()])
    monkeypatch.setattr(viz3d_router, "optimize_conformer", lambda mol, force_field="MMFF": mol)
    monkeypatch.setattr(viz3d_router, "align_molecules_to_reference", lambda mols, reference_idx=0: None)

    it = iter(["PDB_A", "PDB_B", "PDB_C"])
    monkeypatch.setattr(viz3d_router, "conformer_to_pdb", lambda mol, conf_id=0: next(it))

    resp = client.post("/api/viz3d/overlay", json={"smiles_list": ["A", "B", "C"], "reference_idx": 1})
    assert resp.status_code == 200
    data = resp.json()
    assert data["aligned_pdb_strings"] == ["PDB_A", "PDB_B", "PDB_C"]


def test_protein_endpoint_returns_pdb_string(monkeypatch, tmp_path):
    from amprenta_rag.api.routers import viz3d as viz3d_router

    sid = uuid4()
    fp = tmp_path / "x.pdb"
    fp.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n", encoding="utf-8")

    fake_stc = SimpleNamespace(
        id=sid,
        pdb_id="1ABC",
        alphafold_uniprot_id=None,
        source="pdb",
    )
    fake_sf = SimpleNamespace(
        id=uuid4(),
        structure_id=sid,
        file_type="pdb",
        file_path=str(fp),
        file_size_bytes=fp.stat().st_size,
        created_at=None,
    )

    class FakeQuery:
        def __init__(self, obj):
            self.obj = obj

        def filter(self, *args, **kwargs):  # noqa: ANN001
            return self

        def order_by(self, *args, **kwargs):  # noqa: ANN001
            return self

        def first(self):
            return self.obj

    class FakeDB:
        def __init__(self):
            self.calls = 0

        def query(self, *args, **kwargs):  # noqa: ANN001
            self.calls += 1
            if self.calls == 1:
                return FakeQuery(fake_stc)
            return FakeQuery(fake_sf)

    @contextmanager
    def fake_db_session():
        yield FakeDB()

    monkeypatch.setattr(viz3d_router, "db_session", fake_db_session)

    resp = client.get(f"/api/viz3d/protein/{sid}", params={"file_type": "pdb"})
    assert resp.status_code == 200
    data = resp.json()
    assert "ATOM" in data["pdb_string"]
    assert data["metadata"]["structure_id"] == str(sid)


