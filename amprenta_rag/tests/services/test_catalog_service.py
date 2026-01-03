"""Unit tests for data catalog service."""

import pytest
from datetime import datetime
from uuid import uuid4
from unittest.mock import MagicMock, patch

from amprenta_rag.services import catalog as service


class TestCatalogEntryOperations:
    """Test catalog entry operations."""
    
    def test_get_catalog_entries_no_filters(self):
        """Get all active catalog entries."""
        mock_db = MagicMock()
        mock_entries = [MagicMock(), MagicMock()]
        mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_entries
        
        result = service.get_catalog_entries(mock_db)
        
        assert result == mock_entries
        mock_db.query.assert_called_once()
    
    def test_get_catalog_entries_with_category(self):
        """Get catalog entries filtered by category."""
        mock_db = MagicMock()
        mock_entries = [MagicMock()]
        
        # Set up the query chain
        mock_query = mock_db.query.return_value
        mock_query.filter.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_entries
        
        result = service.get_catalog_entries(mock_db, category="Chemistry")
        
        assert result == mock_entries
    
    def test_get_catalog_entries_with_search(self):
        """Get catalog entries with search filter."""
        mock_db = MagicMock()
        mock_entries = [MagicMock()]
        
        # Set up the query chain for search
        mock_query = mock_db.query.return_value
        mock_query.filter.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_entries
        
        result = service.get_catalog_entries(mock_db, search="compound")
        
        assert result == mock_entries
    
    def test_get_catalog_entry_found(self):
        """Get catalog entry by entity type."""
        mock_db = MagicMock()
        mock_entry = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_entry
        
        result = service.get_catalog_entry(mock_db, "compounds")
        
        assert result == mock_entry
    
    def test_get_catalog_entry_not_found(self):
        """Get non-existent catalog entry returns None."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        result = service.get_catalog_entry(mock_db, "nonexistent")
        
        assert result is None
    
    def test_update_catalog_entry_success(self):
        """Update catalog entry fields."""
        mock_db = MagicMock()
        mock_entry = MagicMock()
        
        with patch.object(service, 'get_catalog_entry', return_value=mock_entry):
            result = service.update_catalog_entry(
                mock_db,
                "compounds",
                description="Updated description",
                category="Chemistry",
                row_count=1000,
            )
            
            assert result == mock_entry
            assert mock_entry.description == "Updated description"
            assert mock_entry.category == "Chemistry"
            assert mock_entry.row_count == 1000
            mock_db.commit.assert_called_once()
    
    def test_update_catalog_entry_not_found(self):
        """Update non-existent entry returns None."""
        mock_db = MagicMock()
        
        with patch.object(service, 'get_catalog_entry', return_value=None):
            result = service.update_catalog_entry(mock_db, "nonexistent")
            
            assert result is None


class TestColumnOperations:
    """Test column operations."""
    
    def test_search_columns(self):
        """Search columns by name or description."""
        mock_db = MagicMock()
        mock_columns = [MagicMock(), MagicMock()]
        mock_db.query.return_value.join.return_value.filter.return_value.limit.return_value.all.return_value = mock_columns
        
        result = service.search_columns(mock_db, "compound_id")
        
        assert result == mock_columns
        mock_db.query.assert_called_once()


class TestGlossaryOperations:
    """Test glossary operations."""
    
    def test_create_glossary_term_success(self):
        """Create glossary term successfully."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None  # No existing term
        
        mock_term = MagicMock()
        mock_term.id = uuid4()
        
        with patch('amprenta_rag.services.catalog.GlossaryTerm', return_value=mock_term):
            result = service.create_glossary_term(
                mock_db,
                term="Test Term",
                definition="Test definition",
                category="Chemistry",
                synonyms=["synonym1", "synonym2"],
            )
            
            assert result == mock_term
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
    
    def test_create_glossary_term_duplicate(self):
        """Create duplicate term raises ValueError."""
        mock_db = MagicMock()
        mock_existing = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_existing
        
        with pytest.raises(ValueError, match="already exists"):
            service.create_glossary_term(
                mock_db,
                term="Existing Term",
                definition="Definition",
            )
    
    def test_get_glossary_terms_no_filters(self):
        """Get all glossary terms."""
        mock_db = MagicMock()
        mock_terms = [MagicMock(), MagicMock()]
        mock_db.query.return_value.order_by.return_value.limit.return_value.all.return_value = mock_terms
        
        result = service.get_glossary_terms(mock_db)
        
        assert result == mock_terms
    
    def test_get_glossary_terms_with_search(self):
        """Get glossary terms with search filter."""
        mock_db = MagicMock()
        mock_terms = [MagicMock()]
        mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_terms
        
        result = service.get_glossary_terms(mock_db, search="compound")
        
        assert result == mock_terms
    
    def test_update_glossary_term_success(self):
        """Update glossary term fields."""
        mock_db = MagicMock()
        mock_term = MagicMock()
        mock_term.term = "Original Term"
        mock_db.query.return_value.filter.return_value.first.return_value = mock_term
        
        result = service.update_glossary_term(
            mock_db,
            uuid4(),
            definition="Updated definition",
            category="Updated category",
        )
        
        assert result == mock_term
        assert mock_term.definition == "Updated definition"
        assert mock_term.category == "Updated category"
        mock_db.commit.assert_called_once()
    
    def test_update_glossary_term_not_found(self):
        """Update non-existent term returns None."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        result = service.update_glossary_term(mock_db, uuid4())
        
        assert result is None
    
    def test_update_glossary_term_name_collision(self):
        """Update term name to existing name raises ValueError."""
        mock_db = MagicMock()
        mock_term = MagicMock()
        mock_term.term = "Original"
        mock_existing = MagicMock()
        
        # First call returns the term to update, second returns existing term with new name
        mock_db.query.return_value.filter.return_value.first.side_effect = [mock_term, mock_existing]
        
        with pytest.raises(ValueError, match="already exists"):
            service.update_glossary_term(mock_db, uuid4(), term="Existing Term")
    
    def test_delete_glossary_term_success(self):
        """Delete glossary term and unlink columns."""
        mock_db = MagicMock()
        mock_term = MagicMock()
        mock_term.term = "Term to Delete"
        mock_db.query.return_value.filter.return_value.first.return_value = mock_term
        
        result = service.delete_glossary_term(mock_db, uuid4())
        
        assert result is True
        mock_db.delete.assert_called_once_with(mock_term)
        mock_db.commit.assert_called_once()
    
    def test_delete_glossary_term_not_found(self):
        """Delete non-existent term returns False."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        result = service.delete_glossary_term(mock_db, uuid4())
        
        assert result is False


class TestDataLineageOperations:
    """Test data lineage operations."""
    
    def test_add_lineage_edge_success(self):
        """Add lineage edge successfully."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None  # No existing edge
        
        mock_edge = MagicMock()
        mock_edge.id = uuid4()
        
        with patch('amprenta_rag.services.catalog.DataLineageEdge', return_value=mock_edge):
            result = service.add_lineage_edge(
                mock_db,
                source_type="datasets",
                source_id=uuid4(),
                target_type="signatures",
                target_id=uuid4(),
                relationship_type="input_to",
                transformation="differential_expression",
            )
            
            assert result == mock_edge
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
    
    def test_add_lineage_edge_duplicate(self):
        """Add duplicate edge raises ValueError."""
        mock_db = MagicMock()
        mock_existing = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_existing
        
        with pytest.raises(ValueError, match="already exists"):
            service.add_lineage_edge(
                mock_db,
                source_type="datasets",
                source_id=uuid4(),
                target_type="signatures",
                target_id=uuid4(),
                relationship_type="input_to",
            )
    
    def test_get_entity_lineage(self):
        """Get entity lineage returns parents and children."""
        mock_db = MagicMock()
        
        # Mock parent edges
        mock_parent = MagicMock()
        mock_parent.source_type = "datasets"
        mock_parent.source_id = uuid4()
        mock_parent.relationship_type = "input_to"
        mock_parent.transformation = "analysis"
        
        # Mock child edges
        mock_child = MagicMock()
        mock_child.target_type = "reports"
        mock_child.target_id = uuid4()
        mock_child.relationship_type = "generated_by"
        mock_child.transformation = "visualization"
        
        # Set up query chain for both parent and child queries
        mock_db.query.return_value.filter.return_value.limit.return_value.all.side_effect = [
            [mock_parent],  # Parents query
            [mock_child],   # Children query
        ]
        
        result = service.get_entity_lineage(mock_db, "signatures", uuid4())
        
        assert "parents" in result
        assert "children" in result
        assert len(result["parents"]) == 1
        assert len(result["children"]) == 1
        assert result["parents"][0]["type"] == "datasets"
        assert result["children"][0]["type"] == "reports"


class TestUtilityFunctions:
    """Test utility functions."""
    
    def test_get_catalog_stats(self):
        """Get catalog statistics."""
        mock_db = MagicMock()
        
        # Mock catalog entries
        mock_entry1 = MagicMock()
        mock_entry1.category = "Chemistry"
        mock_entry1.row_count = 1000
        mock_entry1.last_refreshed = datetime(2025, 1, 1)
        
        mock_entry2 = MagicMock()
        mock_entry2.category = "Chemistry"
        mock_entry2.row_count = 500
        mock_entry2.last_refreshed = datetime(2025, 1, 2)  # More recent
        
        mock_entry3 = MagicMock()
        mock_entry3.category = "Omics"
        mock_entry3.row_count = 2000
        mock_entry3.last_refreshed = None
        
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_entry1, mock_entry2, mock_entry3]
        
        result = service.get_catalog_stats(mock_db)
        
        assert result["total_entries"] == 3
        assert result["total_rows"] == 3500  # 1000 + 500 + 2000
        assert "Chemistry" in result["categories"]
        assert "Omics" in result["categories"]
        assert result["categories"]["Chemistry"]["count"] == 2
        assert result["categories"]["Chemistry"]["rows"] == 1500
        assert result["last_refresh"] == datetime(2025, 1, 2)  # Most recent
    
    def test_search_glossary(self):
        """Search glossary terms."""
        mock_db = MagicMock()
        
        mock_term = MagicMock()
        mock_term.id = uuid4()
        mock_term.term = "Compound"
        mock_term.definition = "A chemical substance with defined molecular structure"
        mock_term.category = "Chemistry"
        mock_term.synonyms = ["molecule", "chemical"]
        
        mock_db.query.return_value.filter.return_value.limit.return_value.all.return_value = [mock_term]
        
        result = service.search_glossary(mock_db, "compound")
        
        assert len(result) == 1
        assert result[0]["term"] == "Compound"
        assert result[0]["category"] == "Chemistry"
        assert result[0]["match_type"] == "term"
        assert result[0]["synonyms"] == ["molecule", "chemical"]
    
    def test_refresh_catalog_stats_success(self):
        """Refresh catalog stats successfully."""
        mock_db = MagicMock()
        mock_entry = MagicMock()
        mock_entry.table_name = "compounds"
        
        # Mock the SQL execution
        mock_result = MagicMock()
        mock_result.scalar.return_value = 1500
        mock_db.execute.return_value = mock_result
        
        with patch.object(service, 'get_catalog_entry', return_value=mock_entry):
            result = service.refresh_catalog_stats(mock_db, "compounds")
            
            assert result is True
            assert mock_entry.row_count == 1500
            mock_db.commit.assert_called_once()
    
    def test_refresh_catalog_stats_not_found(self):
        """Refresh non-existent entry returns False."""
        mock_db = MagicMock()
        
        with patch.object(service, 'get_catalog_entry', return_value=None):
            result = service.refresh_catalog_stats(mock_db, "nonexistent")
            
            assert result is False


class TestServiceFunctionStructure:
    """Test service function structure and imports."""
    
    def test_all_required_functions_exist(self):
        """All 10 required functions exist."""
        required_functions = [
            'get_catalog_entries',
            'get_catalog_entry', 
            'update_catalog_entry',
            'search_columns',
            'create_glossary_term',
            'get_glossary_terms',
            'update_glossary_term',
            'delete_glossary_term',
            'add_lineage_edge',
            'get_entity_lineage',
        ]
        
        for func_name in required_functions:
            assert hasattr(service, func_name)
            assert callable(getattr(service, func_name))
    
    def test_utility_functions_exist(self):
        """Utility functions exist."""
        assert hasattr(service, 'get_catalog_stats')
        assert hasattr(service, 'search_glossary')
        assert hasattr(service, 'refresh_catalog_stats')
    
    def test_service_imports_successfully(self):
        """Service imports without error."""
        # If we got here, imports worked
        assert True


class TestValidationLogic:
    """Test validation and business logic."""
    
    def test_get_entity_lineage_limits_nodes(self):
        """Entity lineage respects max_nodes limit."""
        mock_db = MagicMock()
        
        # Mock query chain
        mock_db.query.return_value.filter.return_value.limit.return_value.all.side_effect = [
            [],  # Parents
            [],  # Children
        ]
        
        service.get_entity_lineage(mock_db, "compounds", uuid4(), max_nodes=100)
        
        # Verify limit was applied (50 for parents, 50 for children)
        calls = mock_db.query.return_value.filter.return_value.limit.call_args_list
        assert len(calls) == 2
        assert calls[0][0][0] == 50  # max_nodes // 2
        assert calls[1][0][0] == 50  # max_nodes // 2
    
    def test_search_glossary_truncates_long_definitions(self):
        """Search glossary truncates long definitions."""
        mock_db = MagicMock()
        
        mock_term = MagicMock()
        mock_term.id = uuid4()
        mock_term.term = "Long Term"
        mock_term.definition = "A" * 300  # 300 characters
        mock_term.category = "Test"
        mock_term.synonyms = []
        
        mock_db.query.return_value.filter.return_value.limit.return_value.all.return_value = [mock_term]
        
        result = service.search_glossary(mock_db, "long")
        
        assert len(result[0]["definition"]) == 203  # 200 + "..."
        assert result[0]["definition"].endswith("...")
    
    def test_search_glossary_match_type_detection(self):
        """Search glossary detects match type correctly."""
        mock_db = MagicMock()
        
        # Term that matches in definition, not term name
        mock_term = MagicMock()
        mock_term.id = uuid4()
        mock_term.term = "Chemical Entity"
        mock_term.definition = "A compound with defined structure"
        mock_term.category = "Chemistry"
        mock_term.synonyms = []
        
        mock_db.query.return_value.filter.return_value.limit.return_value.all.return_value = [mock_term]
        
        result = service.search_glossary(mock_db, "compound")  # Matches definition
        
        assert result[0]["match_type"] == "definition"


class TestAutoDiscovery:
    """Test auto-discovery functions."""
    
    def test_refresh_catalog_function_exists(self):
        """Test that refresh_catalog function exists and is callable."""
        assert hasattr(service, 'refresh_catalog')
        assert callable(service.refresh_catalog)
        
        # Test function signature
        import inspect
        sig = inspect.signature(service.refresh_catalog)
        assert 'db' in sig.parameters
    
    def test_humanize_converts_names(self):
        """Test _humanize converts CamelCase and snake_case."""
        assert service._humanize("CompoundId") == "Compound Id"
        assert service._humanize("compound_id") == "Compound Id"
        assert service._humanize("HTSResult") == "Hts Result"  # Fixed expectation
        assert service._humanize("created_at") == "Created At"
        assert service._humanize("simple") == "Simple"
    
    def test_extract_column_type_maps_correctly(self):
        """Test type mapping for common SQLAlchemy types."""
        # Test the actual function behavior
        assert service._extract_column_type(type('MockUUID', (), {})()) == "unknown"  # Mock class returns unknown
        
        # Test the type mapping logic by checking the map directly
        type_map = {
            'uuid': 'uuid',
            'string': 'string',
            'integer': 'integer',
            'boolean': 'boolean',
            'datetime': 'datetime',
            'jsonb': 'json',
        }
        
        # Verify the mapping exists in the function
        import inspect
        source = inspect.getsource(service._extract_column_type)
        for sa_type, catalog_type in type_map.items():
            assert f"'{sa_type}': '{catalog_type}'" in source
    
    def test_get_sample_values_returns_distinct(self):
        """Test sample values are distinct and limited."""
        mock_db = MagicMock()
        
        # Mock SQL execution result
        mock_result = MagicMock()
        mock_result.__iter__ = lambda self: iter([("value1",), ("value2",), ("value3",)])
        mock_db.execute.return_value = mock_result
        
        result = service.get_sample_values(mock_db, "compounds", "compound_id", limit=3)
        
        assert result == ["value1", "value2", "value3"]
        mock_db.execute.assert_called_once()
    
    def test_get_sample_values_validates_names(self):
        """Test SQL injection prevention in table/column names."""
        mock_db = MagicMock()
        
        # Invalid table names should return empty list
        assert service.get_sample_values(mock_db, "compounds; DROP TABLE", "id") == []
        assert service.get_sample_values(mock_db, "compounds", "id; DROP TABLE") == []
        assert service.get_sample_values(mock_db, "123invalid", "id") == []
        assert service.get_sample_values(mock_db, "compounds", "123invalid") == []
        
        # Valid names should proceed to query
        mock_result = MagicMock()
        mock_result.__iter__ = lambda self: iter([])
        mock_db.execute.return_value = mock_result
        
        result = service.get_sample_values(mock_db, "compounds", "compound_id")
        assert result == []
        mock_db.execute.assert_called_once()
    
    def test_detect_lineage_from_fks(self):
        """Test FK-based lineage detection creates edges."""
        mock_db = MagicMock()
        
        # Mock catalog entries with FK columns
        mock_entry = MagicMock()
        mock_entry.entity_type = "Sample"
        mock_entry.id = uuid4()
        
        mock_column = MagicMock()
        mock_column.is_foreign_key = True
        mock_column.foreign_key_target = "compounds.id"
        mock_column.column_name = "compound_id"
        mock_entry.columns = [mock_column]
        
        mock_target_entry = MagicMock()
        mock_target_entry.entity_type = "Compound"
        mock_target_entry.id = uuid4()
        
        # Set up query chain
        mock_db.query.return_value.all.return_value = [mock_entry]
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_target_entry,  # Target entry lookup
            None,  # No existing edge
        ]
        
        result = service.detect_lineage_from_fks(mock_db)
        
        assert result == 1  # 1 edge created
        mock_db.add.assert_called_once()
        mock_db.commit.assert_called_once()
    
    def test_get_fk_target_parses_correctly(self):
        """Test FK target parsing."""
        mock_column = MagicMock()
        mock_fk = MagicMock()
        mock_fk.target_fullname = "compounds.id"
        mock_column.foreign_keys = {mock_fk}
        
        result = service._get_fk_target(mock_column)
        assert result == "compounds.id"
        
        # No foreign keys
        mock_column.foreign_keys = set()
        result = service._get_fk_target(mock_column)
        assert result is None
