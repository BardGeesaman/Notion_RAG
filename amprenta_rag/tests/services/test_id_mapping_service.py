"""Tests for ID mapping service."""

from datetime import datetime, timedelta, timezone
from unittest.mock import MagicMock, patch

from amprenta_rag.services.id_mapping_service import (
    get_mapping,
    get_mappings_batch,
    save_mapping,
    get_mapping_stats,
    migrate_gene_protein_map,
    cleanup_expired_mappings,
)


class TestIDMappingService:
    """Test ID mapping service functionality."""
    
    def test_id_mapping_model_create(self):
        """Test IDMapping model creation and basic CRUD."""
        # Test model creation
        mapping_data = {
            "source_type": "gene",
            "source_id": "BRCA1",
            "target_type": "uniprot",
            "target_id": "P38398",
            "organism": "human",
            "confidence": 0.95
        }
        
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock successful save
            mock_db.query.return_value.filter.return_value.first.return_value = None
            mock_db.add.return_value = None
            mock_db.commit.return_value = None
            mock_db.refresh.return_value = None
            
            # Create mapping
            save_mapping(**mapping_data)
            
            # Verify database operations
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
            
    def test_id_mapping_unique_constraint(self):
        """Test that duplicate mappings are handled correctly."""
        mapping_data = {
            "source_type": "gene",
            "source_id": "BRCA1",
            "target_type": "uniprot",
            "target_id": "P38398",
            "organism": "human"
        }
        
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock existing mapping
            existing_mapping = MagicMock()
            existing_mapping.target_id = "OLD_ID"
            mock_db.query.return_value.filter.return_value.first.return_value = existing_mapping
            
            # Save mapping (should update existing)
            save_mapping(**mapping_data)
            
            # Verify update, not insert
            mock_db.add.assert_not_called()
            mock_db.commit.assert_called_once()
            assert existing_mapping.target_id == "P38398"
    
    def test_get_mapping_from_db(self):
        """Test successful mapping lookup from database."""
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock found mapping
            mock_mapping = MagicMock()
            mock_mapping.target_id = "P38398"
            mock_mapping.expires_at = None
            mock_db.query.return_value.filter.return_value.first.return_value = mock_mapping
            
            result = get_mapping("gene", "BRCA1", "uniprot")
            
            assert result == "P38398"
            mock_db.query.assert_called_once()
    
    def test_get_mapping_not_found(self):
        """Test mapping lookup when not found in database."""
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock no mapping found
            mock_db.query.return_value.filter.return_value.first.return_value = None
            
            result = get_mapping("gene", "NONEXISTENT", "uniprot")
            
            assert result is None
    
    def test_get_mapping_expired(self):
        """Test that expired mappings return None."""
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock expired mapping
            mock_mapping = MagicMock()
            mock_mapping.target_id = "P38398"
            mock_mapping.expires_at = datetime.now(timezone.utc) - timedelta(days=1)  # Expired yesterday
            mock_db.query.return_value.filter.return_value.first.return_value = mock_mapping
            
            result = get_mapping("gene", "BRCA1", "uniprot")
            
            assert result is None
    
    def test_save_mapping_new(self):
        """Test saving a new mapping."""
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock no existing mapping
            mock_db.query.return_value.filter.return_value.first.return_value = None
            
            save_mapping("gene", "BRCA2", "uniprot", "P51587")
            
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
    
    def test_save_mapping_update(self):
        """Test updating an existing mapping."""
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock existing mapping
            existing_mapping = MagicMock()
            existing_mapping.target_id = "OLD_ID"
            mock_db.query.return_value.filter.return_value.first.return_value = existing_mapping
            
            save_mapping("gene", "BRCA2", "uniprot", "P51587")
            
            # Should update, not add new
            mock_db.add.assert_not_called()
            mock_db.commit.assert_called_once()
            assert existing_mapping.target_id == "P51587"
    
    def test_get_mappings_batch(self):
        """Test batch mapping lookup."""
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock mappings
            mapping1 = MagicMock()
            mapping1.source_id = "BRCA1"
            mapping1.target_id = "P38398"
            mapping1.expires_at = None
            
            mapping2 = MagicMock()
            mapping2.source_id = "BRCA2"
            mapping2.target_id = "P51587"
            mapping2.expires_at = None
            
            mock_db.query.return_value.filter.return_value.all.return_value = [mapping1, mapping2]
            
            result = get_mappings_batch(["BRCA1", "BRCA2", "NONEXISTENT"], "gene", "uniprot")
            
            expected = {
                "BRCA1": "P38398",
                "BRCA2": "P51587",
                "NONEXISTENT": None
            }
            assert result == expected
    
    def test_get_mapping_stats(self):
        """Test mapping statistics calculation."""
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock query results
            mock_db.query.return_value.count.return_value = 1000
            
            # Mock grouped stats
            source_stats = [
                MagicMock(source_type="gene", count=500),
                MagicMock(source_type="uniprot", count=500)
            ]
            target_stats = [
                MagicMock(target_type="uniprot", count=600),
                MagicMock(target_type="gene", count=400)
            ]
            
            mock_db.query.return_value.group_by.return_value.all.side_effect = [source_stats, target_stats]
            
            # Mock expired and permanent counts
            mock_db.query.return_value.filter.return_value.count.side_effect = [50, 900]
            
            result = get_mapping_stats()
            
            assert result["total_mappings"] == 1000
            assert result["permanent_mappings"] == 900
            assert result["expired_mappings"] == 50
            assert result["by_source_type"]["gene"] == 500
            assert result["by_target_type"]["uniprot"] == 600
    
    def test_uniprot_adapter_parse_line(self):
        """Test UniProt adapter line parsing."""
        from amprenta_rag.sync.adapters.uniprot_mapping import UniProtMappingAdapter
        
        adapter = UniProtMappingAdapter()
        
        # Test valid line parsing
        test_record = {
            "uniprot_ac": "P38398",
            "gene_id": "672",
            "ensembl": "ENSG00000012048",
            "refseq": "NP_009225.1"
        }
        
        # Test transform_record method
        import asyncio
        mappings = asyncio.run(adapter.transform_record(test_record))
        
        # Should create bidirectional mappings
        assert len(mappings) >= 4  # At least gene->uniprot, ensembl->uniprot, refseq->uniprot, and reverse
        
        # Check specific mappings
        gene_to_uniprot = next((m for m in mappings if m["source_type"] == "gene" and m["target_type"] == "uniprot"), None)
        assert gene_to_uniprot is not None
        assert gene_to_uniprot["source_id"] == "672"
        assert gene_to_uniprot["target_id"] == "P38398"
    
    def test_uniprot_adapter_skip_invalid(self):
        """Test that UniProt adapter skips invalid lines."""
        from amprenta_rag.sync.adapters.uniprot_mapping import UniProtMappingAdapter
        
        adapter = UniProtMappingAdapter()
        
        # Test invalid record (no UniProt AC)
        invalid_record = {
            "gene_id": "672",
            "ensembl": "ENSG00000012048"
        }
        
        import asyncio
        mappings = asyncio.run(adapter.transform_record(invalid_record))
        
        # Should return empty list for invalid records
        assert mappings == []
    
    def test_migrate_gene_protein_map(self):
        """Test migration from legacy gene_protein_map table."""
        with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock legacy data
            legacy_data = [
                ("BRCA1", "P38398"),
                ("BRCA2", "P51587"),
                ("TP53", "P04637")
            ]
            mock_db.execute.return_value.fetchall.return_value = legacy_data
            
            # Mock no existing mappings
            mock_db.query.return_value.filter.return_value.first.return_value = None
            
            result = migrate_gene_protein_map()
            
            # Should create 6 mappings (3 forward + 3 reverse)
            assert mock_db.add.call_count == 6
            mock_db.commit.assert_called_once()
            assert result == 6


# Additional test for cleanup functionality
def test_cleanup_expired_mappings():
    """Test cleanup of expired mappings."""
    with patch("amprenta_rag.services.id_mapping_service.db_session") as mock_session:
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db
        
        # Mock 10 expired mappings deleted
        mock_db.query.return_value.filter.return_value.delete.return_value = 10
        
        result = cleanup_expired_mappings()
        
        assert result == 10
        mock_db.commit.assert_called_once()
