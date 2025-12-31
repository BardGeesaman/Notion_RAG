from __future__ import annotations

import asyncio
from contextlib import contextmanager
from types import SimpleNamespace
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient
from httpx import ASGITransport, AsyncClient

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


class TestAsyncViz3DAPI:
    """Test async execution of viz3d API endpoints."""

    @pytest.mark.asyncio
    async def test_conformers_async(self):
        """Test conformers endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.viz3d._sync_get_conformers') as mock_conformers:
            # Mock the conformers function
            from amprenta_rag.api.schemas import ConformerResponse
            mock_response = ConformerResponse(
                pdb_strings=["ASYNC PDB 1", "ASYNC PDB 2"],
                energies=[1.5, 2.3]
            )
            mock_conformers.return_value = mock_response
            
            transport = ASGITransport(app=app)
            async with AsyncClient(transport=transport, base_url="http://test") as client:
                response = await client.post(
                    "/api/viz3d/conformers",
                    json={
                        "smiles": "CCO",
                        "n_conformers": 2,
                        "optimize": True
                    }
                )
                
                assert response.status_code == 200
                data = response.json()
                assert data["pdb_strings"] == ["ASYNC PDB 1", "ASYNC PDB 2"]
                assert data["energies"] == [1.5, 2.3]
                
                # Verify async execution
                mock_conformers.assert_called_once()
                call_args = mock_conformers.call_args[0][0]
                assert call_args.smiles == "CCO"
                assert call_args.n_conformers == 2
                assert call_args.optimize == True

    @pytest.mark.asyncio
    async def test_overlay_async(self):
        """Test overlay endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.viz3d._sync_overlay') as mock_overlay:
            # Mock the overlay function
            from amprenta_rag.api.schemas import OverlayResponse
            mock_response = OverlayResponse(
                aligned_pdb_strings=["ASYNC ALIGNED 1", "ASYNC ALIGNED 2"]
            )
            mock_overlay.return_value = mock_response
            
            transport = ASGITransport(app=app)
            async with AsyncClient(transport=transport, base_url="http://test") as client:
                response = await client.post(
                    "/api/viz3d/overlay",
                    json={
                        "smiles_list": ["CCO", "CCC"],
                        "reference_idx": 0
                    }
                )
                
                assert response.status_code == 200
                data = response.json()
                assert data["aligned_pdb_strings"] == ["ASYNC ALIGNED 1", "ASYNC ALIGNED 2"]
                
                # Verify async execution
                mock_overlay.assert_called_once()
                call_args = mock_overlay.call_args[0][0]
                assert call_args.smiles_list == ["CCO", "CCC"]
                assert call_args.reference_idx == 0

    @pytest.mark.asyncio
    async def test_protein_async(self):
        """Test protein PDB endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.viz3d._sync_get_protein_pdb') as mock_protein:
            # Mock the protein function
            structure_id = uuid4()
            mock_response = {
                "pdb_string": "ASYNC PROTEIN PDB CONTENT",
                "metadata": {
                    "structure_id": str(structure_id),
                    "pdb_id": "1ABC",
                    "source": "pdb",
                    "file_type": "pdb"
                }
            }
            mock_protein.return_value = mock_response
            
            transport = ASGITransport(app=app)
            async with AsyncClient(transport=transport, base_url="http://test") as client:
                response = await client.get(
                    f"/api/viz3d/protein/{structure_id}?file_type=pdb"
                )
                
                assert response.status_code == 200
                data = response.json()
                assert data["pdb_string"] == "ASYNC PROTEIN PDB CONTENT"
                assert data["metadata"]["pdb_id"] == "1ABC"
                
                # Verify async execution
                mock_protein.assert_called_once_with(structure_id, "pdb")

    @pytest.mark.asyncio
    async def test_concurrent_requests(self):
        """Test multiple simultaneous requests with asyncio.gather()."""
        with patch('amprenta_rag.api.routers.viz3d._sync_get_conformers') as mock_conformers:
            with patch('amprenta_rag.api.routers.viz3d._sync_overlay') as mock_overlay:
                with patch('amprenta_rag.api.routers.viz3d._sync_get_protein_pdb') as mock_protein:
                    # Mock all functions to return different results
                    from amprenta_rag.api.schemas import ConformerResponse, OverlayResponse
                    
                    mock_conformers.return_value = ConformerResponse(
                        pdb_strings=["CONCURRENT CONF"],
                        energies=[1.0]
                    )
                    mock_overlay.return_value = OverlayResponse(
                        aligned_pdb_strings=["CONCURRENT OVERLAY"]
                    )
                    mock_protein.return_value = {
                        "pdb_string": "CONCURRENT PROTEIN",
                        "metadata": {"structure_id": str(uuid4())}
                    }
                    
                    transport = ASGITransport(app=app)
                    async with AsyncClient(transport=transport, base_url="http://test") as client:
                        # Define async request functions
                        async def get_conformers():
                            return await client.post(
                                "/api/viz3d/conformers",
                                json={"smiles": "CCO", "n_conformers": 1, "optimize": False}
                            )
                        
                        async def get_overlay():
                            return await client.post(
                                "/api/viz3d/overlay",
                                json={"smiles_list": ["CCO", "CCC"], "reference_idx": 0}
                            )
                        
                        async def get_protein():
                            return await client.get(f"/api/viz3d/protein/{uuid4()}")
                        
                        # Make 3 concurrent requests
                        tasks = [
                            get_conformers(),
                            get_overlay(),
                            get_protein()
                        ]
                        
                        results = await asyncio.gather(*tasks)
                        
                        # All requests should succeed
                        assert len(results) == 3
                        assert all(r.status_code == 200 for r in results)
                        
                        # Verify responses
                        conformers_data = results[0].json()
                        overlay_data = results[1].json()
                        protein_data = results[2].json()
                        
                        assert conformers_data["pdb_strings"] == ["CONCURRENT CONF"]
                        assert overlay_data["aligned_pdb_strings"] == ["CONCURRENT OVERLAY"]
                        assert protein_data["pdb_string"] == "CONCURRENT PROTEIN"
                        
                        # Verify all calls were made
                        assert mock_conformers.call_count == 1
                        assert mock_overlay.call_count == 1
                        assert mock_protein.call_count == 1


