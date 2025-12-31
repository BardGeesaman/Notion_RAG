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
    log_refresh_start,
    log_refresh_complete,
    get_last_successful_refresh,
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
    
    @patch("amprenta_rag.services.id_mapping_service.db_session")
    def test_log_refresh_start(self, mock_session):
        """Test logging the start of a refresh operation."""
        from amprenta_rag.database.models import MappingRefreshLog
        
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db
        
        # Mock the created log entry
        mock_log = MappingRefreshLog()
        mock_log.id = "test-log-id"
        mock_log.source = "uniprot"
        mock_log.status = "started"
        mock_log.started_at = datetime.now(timezone.utc)
        
        mock_db.refresh = MagicMock()
        mock_db.expunge = MagicMock()
        
        # Mock add/commit/refresh sequence
        def mock_add(log_entry):
            log_entry.id = "test-log-id"
        
        mock_db.add.side_effect = mock_add
        
        result = log_refresh_start("uniprot")
        
        # Verify database operations
        mock_db.add.assert_called_once()
        mock_db.commit.assert_called_once()
        mock_db.refresh.assert_called_once()
        mock_db.expunge.assert_called_once()
        
        # Verify the log entry was created with correct values
        added_log = mock_db.add.call_args[0][0]
        assert added_log.source == "uniprot"
        assert added_log.status == "started"
        assert added_log.started_at is not None
    
    @patch("amprenta_rag.services.id_mapping_service.db_session")
    def test_log_refresh_complete_success(self, mock_session):
        """Test logging successful completion of a refresh operation."""
        from amprenta_rag.database.models import MappingRefreshLog
        
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db
        
        # Mock existing log entry
        mock_log = MappingRefreshLog()
        mock_log.id = "test-log-id"
        mock_log.source = "uniprot"
        mock_log.status = "started"
        
        mock_db.query.return_value.filter.return_value.first.return_value = mock_log
        
        log_refresh_complete("test-log-id", 1000)
        
        # Verify the log was updated correctly
        assert mock_log.status == "success"
        assert mock_log.records_processed == 1000
        assert mock_log.completed_at is not None
        assert mock_log.error_message is None
        
        mock_db.commit.assert_called_once()
    
    @patch("amprenta_rag.services.id_mapping_service.db_session")
    def test_log_refresh_complete_failure(self, mock_session):
        """Test logging failed completion of a refresh operation."""
        from amprenta_rag.database.models import MappingRefreshLog
        
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db
        
        # Mock existing log entry
        mock_log = MappingRefreshLog()
        mock_log.id = "test-log-id"
        mock_log.source = "uniprot"
        mock_log.status = "started"
        
        mock_db.query.return_value.filter.return_value.first.return_value = mock_log
        
        error_msg = "Connection timeout"
        log_refresh_complete("test-log-id", 500, error_msg)
        
        # Verify the log was updated correctly
        assert mock_log.status == "failed"
        assert mock_log.records_processed == 500
        assert mock_log.completed_at is not None
        assert mock_log.error_message == error_msg
        
        mock_db.commit.assert_called_once()
    
    @patch("amprenta_rag.services.id_mapping_service.db_session")
    def test_get_last_successful_refresh(self, mock_session):
        """Test getting the last successful refresh timestamp."""
        from amprenta_rag.database.models import MappingRefreshLog
        
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db
        
        # Mock successful refresh log
        last_refresh_time = datetime.now(timezone.utc) - timedelta(days=1)
        mock_log = MappingRefreshLog()
        mock_log.completed_at = last_refresh_time
        
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = mock_log
        
        result = get_last_successful_refresh("uniprot")
        
        assert result == last_refresh_time
        
        # Verify query was built correctly
        mock_db.query.assert_called_once()
        filter_call = mock_db.query.return_value.filter
        order_by_call = filter_call.return_value.order_by
        first_call = order_by_call.return_value.first
        
        filter_call.assert_called_once()
        order_by_call.assert_called_once()
        first_call.assert_called_once()
    
    @patch("amprenta_rag.services.id_mapping_service.db_session")
    def test_get_last_successful_refresh_none_found(self, mock_session):
        """Test getting last successful refresh when none exists."""
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db
        
        # Mock no refresh log found
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = None
        
        result = get_last_successful_refresh("uniprot")
        
        assert result is None


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
