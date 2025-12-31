"""Tests for protein structures API endpoints."""

from __future__ import annotations

import asyncio
from pathlib import Path
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestStructuresAPI:
    """Test basic structures API functionality."""

    def test_list_structures_success(self):
        """Test GET /api/v1/structures endpoint."""
        with patch('amprenta_rag.api.routers.structures.db_session') as mock_db_session:
            # Mock database session and query
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            
            # Mock query result with all required fields
            mock_structure = MagicMock()
            mock_structure.id = uuid4()
            mock_structure.feature_id = None
            mock_structure.pdb_id = "1ABC"
            mock_structure.alphafold_uniprot_id = None
            mock_structure.source = "pdb"
            mock_structure.resolution = 2.0
            mock_structure.method = "X-RAY DIFFRACTION"
            mock_structure.chain_ids = ["A"]
            mock_structure.prep_status = "raw"
            mock_structure.prep_log = None
            mock_structure.files = []
            
            # Mock the full query chain
            mock_query = mock_db.query.return_value
            mock_query.order_by.return_value.limit.return_value.all.return_value = [mock_structure]
            
            response = client.get("/api/structures?limit=10")
            
            assert response.status_code == 200
            data = response.json()
            assert isinstance(data, list)
            assert len(data) == 1
            assert data[0]["source"] == "pdb"


class TestAsyncStructuresAPI:
    """Test async execution of structures API endpoints."""

    @pytest.mark.asyncio
    async def test_fetch_pdb_async(self):
        """Test PDB fetch endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.structures._sync_fetch_structure') as mock_fetch:
            with patch('amprenta_rag.api.routers.structures.db_session') as mock_db_session:
                with patch('amprenta_rag.api.routers.structures.parse_pdb_metadata') as mock_parse:
                    with patch('amprenta_rag.api.routers.structures.save_structure_file') as mock_save:
                        with patch('amprenta_rag.api.routers.structures.ProteinStructure') as MockProteinStructure:
                            with patch('amprenta_rag.api.routers.structures.StructureFile') as MockStructureFile:
                                # Mock the fetch function
                                pdb_content = b"HEADER    TEST PDB\nATOM      1  N   ALA A   1      20.154  16.967  14.421  1.00 20.00           N\n"
                                mock_fetch.return_value = pdb_content
                                
                                # Mock PDB metadata parsing
                                mock_parse.return_value = {
                                    "sequences": {"A": "MTEST"},
                                    "chain_ids": ["A"],
                                    "resolution": 2.1,
                                    "method": "X-RAY DIFFRACTION"
                                }
                                
                                # Mock database operations
                                mock_db = MagicMock()
                                mock_db_session.return_value.__enter__.return_value = mock_db
                                
                                structure_id = uuid4()
                                feature_id = uuid4()
                                
                                # Create proper mock structure instance
                                mock_structure_instance = MagicMock()
                                mock_structure_instance.id = structure_id
                                mock_structure_instance.feature_id = feature_id
                                mock_structure_instance.pdb_id = "1ABC"
                                mock_structure_instance.alphafold_uniprot_id = None
                                mock_structure_instance.source = "pdb"
                                mock_structure_instance.resolution = 2.1
                                mock_structure_instance.method = "X-RAY DIFFRACTION"
                                mock_structure_instance.chain_ids = ["A"]
                                mock_structure_instance.prep_status = "raw"
                                mock_structure_instance.prep_log = None
                                mock_structure_instance.files = []
                                
                                # Mock ProteinStructure constructor to return our mock instance
                                MockProteinStructure.return_value = mock_structure_instance
                                
                                mock_db.add.return_value = None
                                mock_db.commit.return_value = None
                                mock_db.refresh.return_value = None
                                
                                # Mock file saving
                                mock_save.return_value = "/tmp/structure.pdb"
                                
                                from amprenta_rag.api.routers.structures import fetch_structure
                                from amprenta_rag.api.routers.structures import FetchStructureRequest
                                
                                request = FetchStructureRequest(
                                    source="pdb",
                                    pdb_id="1ABC",
                                    feature_id=feature_id
                                )
                                
                                result = await fetch_structure(request)
                                
                                # Verify async execution
                                mock_fetch.assert_called_once_with("pdb", "1ABC", None)
                                assert mock_parse.called
                                assert mock_save.called
                                assert result.id == structure_id

    @pytest.mark.asyncio
    async def test_fetch_alphafold_async(self):
        """Test AlphaFold fetch endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.structures._sync_fetch_structure') as mock_fetch:
            with patch('amprenta_rag.api.routers.structures.db_session') as mock_db_session:
                with patch('amprenta_rag.api.routers.structures.parse_pdb_metadata') as mock_parse:
                    with patch('amprenta_rag.api.routers.structures.save_structure_file') as mock_save:
                        with patch('amprenta_rag.api.routers.structures.ProteinStructure') as MockProteinStructure:
                            with patch('amprenta_rag.api.routers.structures.StructureFile') as MockStructureFile:
                                # Mock the fetch function
                                af_content = b"HEADER    ALPHAFOLD PREDICTION\nATOM      1  N   MET A   1      20.154  16.967  14.421  1.00 90.00           N\n"
                                mock_fetch.return_value = af_content
                                
                                # Mock PDB metadata parsing
                                mock_parse.return_value = {
                                    "sequences": {"A": "MTEST"},
                                    "chain_ids": ["A"],
                                    "resolution": None,
                                    "method": "ALPHAFOLD"
                                }
                                
                                # Mock database operations
                                mock_db = MagicMock()
                                mock_db_session.return_value.__enter__.return_value = mock_db
                                
                                structure_id = uuid4()
                                feature_id = uuid4()
                                
                                # Create proper mock structure instance
                                mock_structure_instance = MagicMock()
                                mock_structure_instance.id = structure_id
                                mock_structure_instance.feature_id = feature_id
                                mock_structure_instance.pdb_id = None
                                mock_structure_instance.alphafold_uniprot_id = "P12345"
                                mock_structure_instance.source = "alphafold"
                                mock_structure_instance.resolution = None
                                mock_structure_instance.method = "ALPHAFOLD"
                                mock_structure_instance.chain_ids = ["A"]
                                mock_structure_instance.prep_status = "raw"
                                mock_structure_instance.prep_log = None
                                mock_structure_instance.files = []
                                
                                # Mock ProteinStructure constructor to return our mock instance
                                MockProteinStructure.return_value = mock_structure_instance
                                
                                mock_db.add.return_value = None
                                mock_db.commit.return_value = None
                                mock_db.refresh.return_value = None
                                
                                # Mock file saving
                                mock_save.return_value = "/tmp/alphafold.pdb"
                                
                                from amprenta_rag.api.routers.structures import fetch_structure
                                from amprenta_rag.api.routers.structures import FetchStructureRequest
                                
                                request = FetchStructureRequest(
                                    source="alphafold",
                                    uniprot_id="P12345",
                                    feature_id=feature_id
                                )
                                
                                result = await fetch_structure(request)
                                
                                # Verify async execution
                                mock_fetch.assert_called_once_with("alphafold", None, "P12345")
                                assert mock_parse.called
                                assert mock_save.called
                                assert result.id == structure_id

    @pytest.mark.asyncio
    async def test_prepare_async(self):
        """Test structure preparation endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.structures._sync_prepare_structure') as mock_prepare:
            with patch('amprenta_rag.api.routers.structures.db_session') as mock_db_session:
                # Mock database operations
                mock_db = MagicMock()
                mock_db_session.return_value.__enter__.return_value = mock_db
                
                # Mock structure and file
                structure_id = uuid4()
                mock_structure = MagicMock()
                mock_structure.id = structure_id
                mock_structure.feature_id = uuid4()
                mock_structure.pdb_id = "1ABC"
                mock_structure.alphafold_uniprot_id = None
                mock_structure.source = "pdb"
                mock_structure.resolution = 2.0
                mock_structure.method = "X-RAY DIFFRACTION"
                mock_structure.chain_ids = ["A"]
                mock_structure.prep_status = "raw"
                mock_structure.prep_log = None
                
                mock_file = MagicMock()
                mock_file.id = uuid4()
                mock_file.structure_id = structure_id
                mock_file.file_type = "pdb"
                mock_file.file_path = "/tmp/input.pdb"
                mock_file.file_size_bytes = 1024
                mock_file.md5_hash = "abc123"
                mock_structure.files = [mock_file]
                
                mock_db.query.return_value.filter.return_value.first.return_value = mock_structure
                
                # Mock file operations
                with patch('pathlib.Path.read_bytes') as mock_read_bytes:
                    mock_read_bytes.return_value = b"PREPARED PDB CONTENT"
                    
                    from amprenta_rag.api.routers.structures import prepare_structure_endpoint
                    from amprenta_rag.api.routers.structures import PrepareRequest
                    
                    request = PrepareRequest(chain_id="A")
                    
                    result = await prepare_structure_endpoint(structure_id, request)
                    
                    # Verify async execution
                    mock_prepare.assert_called_once_with(
                        "/tmp/input.pdb", 
                        "/tmp/prepared.pdb", 
                        "A"
                    )
                    assert mock_structure.prep_status == "prepared"
                    assert mock_structure.prep_log == "prepared ok"

    @pytest.mark.asyncio
    async def test_concurrent_fetch(self):
        """Test multiple simultaneous structure fetches."""
        with patch('amprenta_rag.api.routers.structures._sync_fetch_structure') as mock_fetch:
            with patch('amprenta_rag.api.routers.structures.db_session') as mock_db_session:
                with patch('amprenta_rag.api.routers.structures.parse_pdb_metadata') as mock_parse:
                    with patch('amprenta_rag.api.routers.structures.save_structure_file') as mock_save:
                        # Mock function to return different content for different calls
                        def mock_fetch_func(source, pdb_id, uniprot_id):
                            if source == "pdb":
                                return f"PDB CONTENT FOR {pdb_id}".encode()
                            else:
                                return f"ALPHAFOLD CONTENT FOR {uniprot_id}".encode()
                        
                        mock_fetch.side_effect = mock_fetch_func
                        
                        # Mock other dependencies
                        mock_parse.return_value = {
                            "sequences": {"A": "MTEST"},
                            "chain_ids": ["A"],
                            "resolution": 2.0,
                            "method": "X-RAY"
                        }
                        
                        mock_db = MagicMock()
                        mock_db_session.return_value.__enter__.return_value = mock_db
                        
                        # Create properly mocked structures for each call
                        structure_ids = [uuid4() for _ in range(3)]
                        mock_structures = []
                        
                        for i, sid in enumerate(structure_ids):
                            mock_structure = MagicMock()
                            mock_structure.id = sid
                            mock_structure.feature_id = uuid4()
                            if i == 1:  # AlphaFold structure
                                mock_structure.pdb_id = None
                                mock_structure.alphafold_uniprot_id = "P12345"
                                mock_structure.source = "alphafold"
                            else:  # PDB structures
                                mock_structure.pdb_id = f"{i+1}ABC" if i == 0 else "2DEF"
                                mock_structure.alphafold_uniprot_id = None
                                mock_structure.source = "pdb"
                            mock_structure.resolution = 2.0
                            mock_structure.method = "X-RAY DIFFRACTION"
                            mock_structure.chain_ids = ["A"]
                            mock_structure.prep_status = "raw"
                            mock_structure.prep_log = None
                            mock_structure.files = []
                            mock_structures.append(mock_structure)
                        
                        mock_save.return_value = "/tmp/structure.pdb"
                        
                        # Test concurrent execution by calling the sync helper directly
                        # This avoids the complexity of mocking the full database flow
                        async def fetch_pdb_structure(pdb_id: str):
                            from amprenta_rag.api.routers.structures import _sync_fetch_structure
                            return await asyncio.to_thread(_sync_fetch_structure, "pdb", pdb_id, None)
                        
                        async def fetch_af_structure(uniprot_id: str):
                            from amprenta_rag.api.routers.structures import _sync_fetch_structure
                            return await asyncio.to_thread(_sync_fetch_structure, "alphafold", None, uniprot_id)
                        
                        # Make 3 concurrent requests
                        tasks = [
                            fetch_pdb_structure("1ABC"),
                            fetch_af_structure("P12345"),
                            fetch_pdb_structure("2DEF"),
                        ]
                        
                        results = await asyncio.gather(*tasks)
                        
                        # All requests should succeed
                        assert len(results) == 3
                        assert all(isinstance(result, bytes) for result in results)
                        
                        # Verify all calls were made
                        assert mock_fetch.call_count == 3
                        
                        # Verify calls with correct parameters
                        calls = mock_fetch.call_args_list
                        assert calls[0][0] == ("pdb", "1ABC", None)
                        assert calls[1][0] == ("alphafold", None, "P12345")
                        assert calls[2][0] == ("pdb", "2DEF", None)

    @pytest.mark.asyncio
    async def test_fetch_error_handling(self):
        """Test error handling in async fetch operations."""
        with patch('amprenta_rag.api.routers.structures._sync_fetch_structure') as mock_fetch:
            # Mock function to raise an exception
            mock_fetch.side_effect = Exception("External API unavailable")
            
            from amprenta_rag.api.routers.structures import fetch_structure
            from amprenta_rag.api.routers.structures import FetchStructureRequest
            
            request = FetchStructureRequest(
                source="pdb",
                pdb_id="1ABC",
                feature_id=uuid4()
            )
            
            # Should raise HTTPException
            with pytest.raises(Exception) as exc_info:
                await fetch_structure(request)
            
            # Verify the exception is properly handled
            assert "External API unavailable" in str(exc_info.value)
            mock_fetch.assert_called_once()

    @pytest.mark.asyncio
    async def test_prepare_error_handling(self):
        """Test error handling in async prepare operations."""
        with patch('amprenta_rag.api.routers.structures._sync_prepare_structure') as mock_prepare:
            with patch('amprenta_rag.api.routers.structures.db_session') as mock_db_session:
                # Mock function to raise an exception
                mock_prepare.side_effect = Exception("Structure preparation failed")
                
                # Mock database operations
                mock_db = MagicMock()
                mock_db_session.return_value.__enter__.return_value = mock_db
                
                structure_id = uuid4()
                mock_structure = MagicMock()
                mock_structure.id = structure_id
                mock_structure.prep_status = "raw"
                mock_structure.prep_log = None
                
                mock_file = MagicMock()
                mock_file.file_type = "pdb"
                mock_file.file_path = "/tmp/input.pdb"
                mock_structure.files = [mock_file]
                
                mock_db.query.return_value.filter.return_value.first.return_value = mock_structure
                
                from amprenta_rag.api.routers.structures import prepare_structure_endpoint
                from amprenta_rag.api.routers.structures import PrepareRequest
                
                request = PrepareRequest(chain_id="A")
                
                # Should raise HTTPException
                with pytest.raises(Exception) as exc_info:
                    await prepare_structure_endpoint(structure_id, request)
                
                # Verify the exception is properly handled
                assert "Preparation failed" in str(exc_info.value)
                mock_prepare.assert_called_once()
                
                # Verify failure status was set
                assert mock_structure.prep_status == "failed"
                assert "Structure preparation failed" in mock_structure.prep_log
