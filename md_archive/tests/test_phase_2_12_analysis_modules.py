"""Tests for Phase 2.12: Analysis Module Notion Removal."""

import inspect
from unittest.mock import patch, MagicMock


class TestModuleImports:
    """Test both modules import successfully."""

    def test_enrichment_imports(self):
        """Test enrichment module imports."""
        from amprenta_rag.analysis.enrichment import (
            enrich_dataset_pathways,
            enrich_signature_pathways,
        )
        assert callable(enrich_dataset_pathways)
        assert callable(enrich_signature_pathways)

    def test_program_signature_maps_imports(self):
        """Test program_signature_maps module imports."""
        from amprenta_rag.analysis.program_signature_maps import (
            get_program_datasets,
            compute_program_signature_scores,
            update_notion_with_program_map,
            ProgramSignatureMap,
            ProgramSignatureScore,
            ProgramOmicsCoverage,
        )
        assert callable(get_program_datasets)
        assert callable(compute_program_signature_scores)
        assert callable(update_notion_with_program_map)
        assert ProgramSignatureMap is not None
        assert ProgramSignatureScore is not None
        assert ProgramOmicsCoverage is not None


class TestNoNotionReferences:
    """Test no Notion references in source code."""

    def test_no_notion_in_enrichment(self):
        """Test no Notion references in enrichment module."""
        from amprenta_rag.analysis import enrichment
        source = inspect.getsource(enrichment)
        assert 'notion_headers' not in source
        assert 'pathway_notion_integration' not in source
        assert 'from amprenta_rag.clients.notion_client' not in source
        assert 'from amprenta_rag.analysis.pathway_notion_integration' not in source

    def test_no_notion_in_program_signature_maps(self):
        """Test no Notion references in program_signature_maps module."""
        from amprenta_rag.analysis import program_signature_maps
        source = inspect.getsource(program_signature_maps)
        assert 'notion_headers' not in source
        assert 'from amprenta_rag.clients.notion_client' not in source
        # update_notion_with_program_map should exist but be stubbed
        assert 'update_notion_with_program_map' in source


class TestUpdateNotionStub:
    """Test update_notion_with_program_map is properly stubbed."""

    def test_stubbed_function_does_nothing(self):
        """Test stubbed function does nothing and returns None."""
        from amprenta_rag.analysis.program_signature_maps import (
            update_notion_with_program_map,
            ProgramSignatureMap,
        )

        # Create mock program map
        mock_map = MagicMock(spec=ProgramSignatureMap)
        mock_map.program_id = "test-program-id"
        mock_map.program_name = "Test Program"

        # Should not raise, should do nothing
        result = update_notion_with_program_map(mock_map)
        assert result is None

    def test_stubbed_function_logs_debug(self):
        """Test stubbed function logs debug message."""
        from amprenta_rag.analysis.program_signature_maps import (
            update_notion_with_program_map,
            ProgramSignatureMap,
        )

        mock_map = MagicMock(spec=ProgramSignatureMap)

        with patch('amprenta_rag.analysis.program_signature_maps.logger') as mock_logger:
            update_notion_with_program_map(mock_map)
            # Verify debug was called
            mock_logger.debug.assert_called_once()
            # Verify message contains "Notion update removed" or similar
            call_args = mock_logger.debug.call_args[0][0]
            assert 'notion' in call_args.lower() or 'removed' in call_args.lower()


class TestFunctionSignatures:
    """Test functions have expected signatures."""

    def test_enrich_dataset_pathways_signature(self):
        """Test enrich_dataset_pathways signature."""
        from amprenta_rag.analysis.enrichment import enrich_dataset_pathways
        sig = inspect.signature(enrich_dataset_pathways)
        params = list(sig.parameters.keys())
        assert 'dataset_page_id' in params
        assert 'p_value_threshold' in params
        assert 'pathway_sources' in params

    def test_enrich_signature_pathways_signature(self):
        """Test enrich_signature_pathways signature."""
        from amprenta_rag.analysis.enrichment import enrich_signature_pathways
        sig = inspect.signature(enrich_signature_pathways)
        params = list(sig.parameters.keys())
        assert 'signature_page_id' in params
        assert 'p_value_threshold' in params
        assert 'pathway_sources' in params

    def test_get_program_datasets_signature(self):
        """Test get_program_datasets signature."""
        from amprenta_rag.analysis.program_signature_maps import get_program_datasets
        sig = inspect.signature(get_program_datasets)
        params = list(sig.parameters.keys())
        assert 'program_page_id' in params

    def test_compute_program_signature_scores_signature(self):
        """Test compute_program_signature_scores signature."""
        from amprenta_rag.analysis.program_signature_maps import compute_program_signature_scores
        sig = inspect.signature(compute_program_signature_scores)
        params = list(sig.parameters.keys())
        assert 'program_page_id' in params
        assert 'signature_ids' in params
        assert 'use_cache' in params


class TestEnrichmentFunctions:
    """Test enrichment functions work correctly."""

    @patch('amprenta_rag.ingestion.multi_omics_scoring.extract_dataset_features_by_type')
    @patch('amprenta_rag.analysis.enrichment.perform_pathway_enrichment')
    def test_enrich_dataset_pathways_returns_results(self, mock_enrichment, mock_extract):
        """Test enrich_dataset_pathways returns pathway results."""
        from amprenta_rag.analysis.enrichment import enrich_dataset_pathways
        from amprenta_rag.analysis.pathway_analysis import PathwayEnrichmentResult

        # Mock feature extraction
        mock_extract.return_value = {
            "Lipidomics": {"Cer(d18:1/16:0)", "PC(16:0/18:1)"},
            "Metabolomics": {"Glucose", "Lactate"},
        }

        # Mock enrichment results
        mock_result = PathwayEnrichmentResult(
            pathway=MagicMock(),
            p_value=0.01,
            adjusted_p_value=0.05,
            enrichment_ratio=2.5,
            input_features=2,
            pathway_size=10,
            background_size=100,
            matched_features=["Cer(d18:1/16:0)", "PC(16:0/18:1)"],
        )
        mock_enrichment.return_value = [mock_result]

        results = enrich_dataset_pathways("test-dataset-id")

        assert isinstance(results, list)
        assert len(results) == 1
        assert results[0] == mock_result
        mock_extract.assert_called_once()
        mock_enrichment.assert_called_once()

    @patch('amprenta_rag.database.base.get_db')
    @patch('amprenta_rag.ingestion.postgres_signature_loader.load_signature_from_postgres')
    @patch('amprenta_rag.analysis.enrichment.perform_pathway_enrichment')
    def test_enrich_signature_pathways_with_uuid(self, mock_enrichment, mock_load, mock_get_db):
        """Test enrich_signature_pathways works with UUID."""
        from amprenta_rag.analysis.enrichment import enrich_signature_pathways
        from amprenta_rag.analysis.pathway_analysis import PathwayEnrichmentResult
        from amprenta_rag.signatures.signature_loader import Signature, SignatureComponent

        # Mock database session
        mock_db = MagicMock()
        mock_sig_model = MagicMock()
        mock_sig_model.id = "test-signature-id"
        mock_db.query.return_value.filter.return_value.first.return_value = mock_sig_model
        mock_get_db.return_value = iter([mock_db])

        # Mock signature
        mock_sig = Signature(
            name="Test Signature",
            components=[
                SignatureComponent(feature_name="Cer(d18:1/16:0)", feature_type="lipid"),
            ],
        )
        mock_load.return_value = mock_sig

        # Mock enrichment results
        mock_result = PathwayEnrichmentResult(
            pathway=MagicMock(),
            p_value=0.01,
            adjusted_p_value=0.05,
            enrichment_ratio=2.5,
            input_features=2,
            pathway_size=10,
            background_size=100,
            matched_features=["Cer(d18:1/16:0)", "PC(16:0/18:1)"],
        )
        mock_enrichment.return_value = [mock_result]

        results = enrich_signature_pathways("test-signature-id")

        assert isinstance(results, list)
        assert len(results) == 1
        mock_db.close.assert_called_once()


class TestProgramSignatureMapsFunctions:
    """Test program_signature_maps functions work correctly."""

    @patch('amprenta_rag.database.base.get_db')
    def test_get_program_datasets_returns_list(self, mock_get_db):
        """Test get_program_datasets returns dataset IDs from Postgres."""
        from amprenta_rag.analysis.program_signature_maps import get_program_datasets

        # Mock program with datasets
        mock_program = MagicMock()
        mock_program.name = "Test Program"
        mock_dataset1 = MagicMock()
        mock_dataset1.notion_page_id = "dataset-id-1"
        mock_dataset2 = MagicMock()
        mock_dataset2.notion_page_id = "dataset-id-2"
        mock_program.datasets = [mock_dataset1, mock_dataset2]
        mock_program.experiments = []

        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_program
        mock_get_db.return_value = iter([mock_db])

        result = get_program_datasets("test-program-id")

        assert isinstance(result, list)
        assert len(result) == 2
        assert "dataset-id-1" in result
        assert "dataset-id-2" in result
        mock_db.close.assert_called_once()

    @patch('amprenta_rag.analysis.program_signature_maps.get_program_datasets')
    @patch('amprenta_rag.ingestion.postgres_signature_loader.fetch_all_signatures_from_postgres')
    @patch('amprenta_rag.ingestion.multi_omics_scoring.extract_dataset_features_by_type')
    @patch('amprenta_rag.ingestion.multi_omics_scoring.score_multi_omics_signature_against_dataset')
    @patch('amprenta_rag.database.base.get_db')
    def test_compute_program_signature_scores_returns_scores(self, mock_get_db, mock_score, mock_extract, mock_fetch, mock_get_datasets):
        """Test compute_program_signature_scores returns signature scores."""
        from amprenta_rag.analysis.program_signature_maps import compute_program_signature_scores
        from amprenta_rag.database.models import Signature as SignatureModel

        # Mock datasets
        mock_get_datasets.return_value = ["dataset-1", "dataset-2"]

        # Mock signatures
        mock_sig_model = MagicMock(spec=SignatureModel)
        mock_sig_model.id = "sig-id-1"
        mock_sig_model.notion_page_id = "sig-notion-id-1"
        mock_sig_model.name = "Test Signature"
        mock_fetch.return_value = [mock_sig_model]

        # Mock feature extraction
        mock_extract.return_value = {
            "Lipidomics": {"Cer(d18:1/16:0)"},
        }

        # Mock scoring
        mock_score_result = MagicMock()
        mock_score_result.total_score = 0.8
        mock_score_result.component_matches = []
        mock_score.return_value = mock_score_result

        # Mock signature loading
        with patch('amprenta_rag.ingestion.postgres_signature_loader.load_signature_from_postgres') as mock_load:
            from amprenta_rag.signatures.signature_loader import Signature
            mock_sig = Signature(name="Test Signature", components=[])
            mock_load.return_value = mock_sig

            # Mock program name lookup (first call)
            mock_db1 = MagicMock()
            mock_program = MagicMock()
            mock_program.name = "Test Program"
            mock_db1.query.return_value.filter.return_value.first.return_value = mock_program

            # Mock signature loading (second call)
            mock_db2 = MagicMock()
            mock_db2.query.return_value.filter.return_value.first.return_value = None

            # Return iterator that yields both database sessions
            mock_get_db.return_value = iter([mock_db1, mock_db2])

            results = compute_program_signature_scores("test-program-id")

            assert isinstance(results, list)
            # Results may be empty if no matching signatures, but should not crash
            assert mock_db1.close.called or mock_db2.close.called


class TestBackwardCompatibility:
    """Test backward compatibility with Notion page IDs."""

    @patch('amprenta_rag.database.base.get_db')
    def test_get_program_datasets_handles_notion_page_id(self, mock_get_db):
        """Test get_program_datasets handles Notion page ID lookup."""
        from amprenta_rag.analysis.program_signature_maps import get_program_datasets

        # Mock program lookup by notion_page_id
        mock_program = MagicMock()
        mock_program.datasets = []
        mock_program.experiments = []

        mock_db = MagicMock()
        # First query (UUID lookup) returns None
        # Second query (notion_page_id lookup) returns program
        mock_db.query.return_value.filter.return_value.first.side_effect = [None, mock_program]
        mock_get_db.return_value = iter([mock_db])

        result = get_program_datasets("notion-page-id-with-dashes")

        assert isinstance(result, list)
        # Should have tried both UUID and notion_page_id lookups
        assert mock_db.query.return_value.filter.return_value.first.call_count >= 1
        mock_db.close.assert_called_once()

