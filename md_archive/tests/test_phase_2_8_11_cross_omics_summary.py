"""Tests for Phase 2.8-2.11: Cross-Omics Summary Modules."""

import inspect
import pytest
from uuid import uuid4, UUID
from unittest.mock import MagicMock, patch


class TestModuleImports:
    """Test all modules import successfully."""

    def test_program_summary_imports(self):
        """Test program_summary module imports."""
        from amprenta_rag.query.cross_omics.program_summary import cross_omics_program_summary
        assert callable(cross_omics_program_summary)

    def test_signature_summary_imports(self):
        """Test signature_summary module imports."""
        from amprenta_rag.query.cross_omics.signature_summary import cross_omics_signature_summary
        assert callable(cross_omics_signature_summary)

    def test_dataset_summary_imports(self):
        """Test dataset_summary module imports."""
        from amprenta_rag.query.cross_omics.dataset_summary import cross_omics_dataset_summary
        assert callable(cross_omics_dataset_summary)

    def test_feature_summary_imports(self):
        """Test feature_summary module imports."""
        from amprenta_rag.query.cross_omics.feature_summary import cross_omics_feature_summary
        assert callable(cross_omics_feature_summary)


class TestNoNotionReferences:
    """Test no Notion references in source code."""

    def test_no_notion_in_program_summary(self):
        """Test no Notion references in program_summary."""
        from amprenta_rag.query.cross_omics import program_summary_postgres as program_summary
        source = inspect.getsource(program_summary)
        assert 'fetch_notion_page' not in source
        assert 'notion_headers' not in source
        assert 'from amprenta_rag.clients.notion_client' not in source

    def test_no_notion_in_signature_summary(self):
        """Test no Notion references in signature_summary."""
        from amprenta_rag.query.cross_omics import signature_summary_postgres as signature_summary
        source = inspect.getsource(signature_summary)
        assert 'fetch_notion_page' not in source
        assert 'notion_headers' not in source
        assert 'from amprenta_rag.clients.notion_client' not in source

    def test_no_notion_in_dataset_summary(self):
        """Test no Notion references in dataset_summary."""
        from amprenta_rag.query.cross_omics import dataset_summary_postgres as dataset_summary
        source = inspect.getsource(dataset_summary)
        assert 'fetch_notion_page' not in source
        assert 'notion_headers' not in source
        assert 'from amprenta_rag.clients.notion_client' not in source

    def test_no_notion_in_feature_summary(self):
        """Test no Notion references in feature_summary."""
        from amprenta_rag.query.cross_omics import feature_summary_postgres as feature_summary
        source = inspect.getsource(feature_summary)
        assert 'fetch_notion_page' not in source
        assert 'notion_headers' not in source
        assert 'from amprenta_rag.clients.notion_client' not in source


class TestFunctionSignatures:
    """Test function signatures accept UUID."""

    def test_program_summary_accepts_uuid(self):
        """Test program_summary accepts UUID parameter."""
        from amprenta_rag.query.cross_omics.program_summary import cross_omics_program_summary
        sig = inspect.signature(cross_omics_program_summary)
        params = list(sig.parameters.keys())
        assert 'program_id' in params
        assert 'db' in params
        # Verify program_id type hint is UUID
        param = sig.parameters['program_id']
        assert param.annotation == UUID or 'UUID' in str(param.annotation)

    def test_signature_summary_accepts_uuid(self):
        """Test signature_summary accepts UUID parameter."""
        from amprenta_rag.query.cross_omics.signature_summary import cross_omics_signature_summary
        sig = inspect.signature(cross_omics_signature_summary)
        params = list(sig.parameters.keys())
        assert 'signature_id' in params
        assert 'db' in params
        # Verify signature_id type hint is UUID
        param = sig.parameters['signature_id']
        assert param.annotation == UUID or 'UUID' in str(param.annotation)

    def test_dataset_summary_accepts_uuid(self):
        """Test dataset_summary accepts UUID parameter."""
        from amprenta_rag.query.cross_omics.dataset_summary import cross_omics_dataset_summary
        sig = inspect.signature(cross_omics_dataset_summary)
        params = list(sig.parameters.keys())
        assert 'dataset_id' in params
        assert 'db' in params
        # Verify dataset_id type hint is UUID
        param = sig.parameters['dataset_id']
        assert param.annotation == UUID or 'UUID' in str(param.annotation)

    def test_feature_summary_accepts_name(self):
        """Test feature_summary accepts feature_name parameter."""
        from amprenta_rag.query.cross_omics.feature_summary import cross_omics_feature_summary
        sig = inspect.signature(cross_omics_feature_summary)
        params = list(sig.parameters.keys())
        assert 'feature_name' in params
        assert 'feature_type' in params
        assert 'db' in params


class TestNonExistentEntity:
    """Test handling of non-existent entities."""

    @patch('amprenta_rag.query.cross_omics.program_summary.get_db')
    def test_program_not_found(self, mock_get_db):
        """Test program_summary handles non-existent program."""
        from amprenta_rag.query.cross_omics.program_summary import cross_omics_program_summary

        # Mock database session
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        mock_get_db.return_value = iter([mock_db])

        result = cross_omics_program_summary(program_id=uuid4())

        assert isinstance(result, str)
        assert 'not found' in result.lower() or 'error' in result.lower()
        mock_db.close.assert_called_once()

    @patch('amprenta_rag.query.cross_omics.signature_summary.get_db')
    def test_signature_not_found(self, mock_get_db):
        """Test signature_summary handles non-existent signature."""
        from amprenta_rag.query.cross_omics.signature_summary import cross_omics_signature_summary

        # Mock database session
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        mock_get_db.return_value = iter([mock_db])

        result = cross_omics_signature_summary(signature_id=uuid4())

        assert isinstance(result, str)
        assert 'not found' in result.lower() or 'error' in result.lower()
        mock_db.close.assert_called_once()

    @patch('amprenta_rag.query.cross_omics.dataset_summary.get_db')
    def test_dataset_not_found(self, mock_get_db):
        """Test dataset_summary handles non-existent dataset."""
        from amprenta_rag.query.cross_omics.dataset_summary import cross_omics_dataset_summary

        # Mock database session
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        mock_get_db.return_value = iter([mock_db])

        result = cross_omics_dataset_summary(dataset_id=uuid4())

        assert isinstance(result, str)
        assert 'not found' in result.lower() or 'error' in result.lower()
        mock_db.close.assert_called_once()

    @patch('amprenta_rag.query.cross_omics.feature_summary.get_db')
    def test_feature_not_found(self, mock_get_db):
        """Test feature_summary handles non-existent feature."""
        from amprenta_rag.query.cross_omics.feature_summary import cross_omics_feature_summary

        # Mock database session
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        mock_get_db.return_value = iter([mock_db])

        # Mock Pinecone query to return empty results
        with patch('amprenta_rag.query.cross_omics.feature_summary.query_pinecone', return_value=[]):
            result = cross_omics_feature_summary(feature_name="NonExistentFeature", feature_type="gene")

        assert isinstance(result, str)
        assert 'not found' in result.lower() or 'no sufficient' in result.lower() or 'no relevant' in result.lower()
        mock_db.close.assert_called_once()


class TestScopeFixes:
    """Test scope fixes prevent UnboundLocalError."""

    @patch('amprenta_rag.query.cross_omics.program_summary.get_db')
    @patch('amprenta_rag.query.cross_omics.program_summary.retrieve_chunks_for_objects')
    def test_program_summary_no_unbound_local_error(self, mock_retrieve, mock_get_db):
        """Test program_summary doesn't raise UnboundLocalError on exception."""
        from amprenta_rag.query.cross_omics.program_summary import cross_omics_program_summary

        # Mock database session that raises exception
        mock_db = MagicMock()
        mock_db.query.side_effect = Exception("Database error")
        mock_get_db.return_value = iter([mock_db])
        mock_retrieve.return_value = []

        # Should not raise UnboundLocalError - variables should be initialized
        # Other exceptions are OK, but UnboundLocalError indicates scope issue
        try:
            result = cross_omics_program_summary(program_id=uuid4())
            # If we get here, no UnboundLocalError occurred
            assert isinstance(result, str)
        except UnboundLocalError as e:
            pytest.fail(f"UnboundLocalError raised - scope fix not working: {e}")
        except Exception:
            # Other exceptions are OK - the important thing is no UnboundLocalError
            pass
        finally:
            # Verify session cleanup happens even on exception
            mock_db.close.assert_called_once()

    @patch('amprenta_rag.query.cross_omics.signature_summary.get_db')
    @patch('amprenta_rag.query.cross_omics.signature_summary.query_pinecone')
    def test_signature_summary_no_unbound_local_error(self, mock_query, mock_get_db):
        """Test signature_summary doesn't raise UnboundLocalError on exception."""
        from amprenta_rag.query.cross_omics.signature_summary import cross_omics_signature_summary

        # Mock database session that raises exception
        mock_db = MagicMock()
        mock_db.query.side_effect = Exception("Database error")
        mock_get_db.return_value = iter([mock_db])
        mock_query.return_value = []

        # Should not raise UnboundLocalError
        # Other exceptions are OK, but UnboundLocalError indicates scope issue
        try:
            result = cross_omics_signature_summary(signature_id=uuid4())
            assert isinstance(result, str)
        except UnboundLocalError as e:
            pytest.fail(f"UnboundLocalError raised - scope fix not working: {e}")
        except Exception:
            # Other exceptions are OK - the important thing is no UnboundLocalError
            pass
        finally:
            # Verify session cleanup happens even on exception
            mock_db.close.assert_called_once()


class TestPostgresIntegration:
    """Test Postgres integration works correctly."""

    @patch('amprenta_rag.query.cross_omics.program_summary.get_db')
    @patch('amprenta_rag.query.cross_omics.program_summary.retrieve_chunks_for_objects')
    @patch('amprenta_rag.query.cross_omics.program_summary.synthesize_cross_omics_summary')
    def test_program_summary_with_postgres_data(self, mock_synthesize, mock_retrieve, mock_get_db):
        """Test program_summary works with Postgres data."""
        from amprenta_rag.query.cross_omics.program_summary import cross_omics_program_summary

        # Mock program with relationships
        mock_program = MagicMock()
        mock_program.name = "Test Program"
        mock_program.experiments = []
        mock_program.datasets = []

        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_program
        mock_get_db.return_value = iter([mock_db])
        mock_retrieve.return_value = []
        mock_synthesize.return_value = "Test summary"

        result = cross_omics_program_summary(program_id=uuid4())

        assert isinstance(result, str)
        # Should return early with "No sufficient multi-omics context" message
        assert 'no sufficient' in result.lower() or result == "Test summary"
        mock_db.close.assert_called_once()

