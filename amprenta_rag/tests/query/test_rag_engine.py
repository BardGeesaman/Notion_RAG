"""
Tests for the RAG query engine module.

Coverage:
- RAGQueryResult and MatchSummary dataclasses
- Match processing (summarize_match, filter_matches)
- Chunk collection
- Answer synthesis
- Main query_rag function
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from amprenta_rag.query.rag.models import MatchSummary, RAGQueryResult
from amprenta_rag.query.rag.match_processing import (
    summarize_match,
    filter_matches,
    _get_source,
    _format_snippet,
    _get_tags,
)


class TestMatchSummary:
    """Tests for MatchSummary dataclass."""

    def test_match_summary_creation(self):
        """MatchSummary should be created with all required fields."""
        summary = MatchSummary(
            id="test-id",
            score=0.95,
            source="literature",
            title="Test Paper",
            snippet="This is a test snippet.",
            tags=["ALS", "biomarker"],
            metadata={"key": "value"},
        )
        assert summary.id == "test-id"
        assert summary.score == 0.95
        assert summary.source == "literature"
        assert summary.title == "Test Paper"
        assert summary.snippet == "This is a test snippet."
        assert summary.tags == ["ALS", "biomarker"]
        assert summary.metadata == {"key": "value"}


class TestRAGQueryResult:
    """Tests for RAGQueryResult dataclass."""

    def test_rag_query_result_creation(self):
        """RAGQueryResult should be created with all required fields."""
        match = MatchSummary(
            id="m1",
            score=0.9,
            source="lit",
            title="Title",
            snippet="Snippet",
            tags=[],
            metadata={},
        )
        result = RAGQueryResult(
            query="test query",
            matches=[match],
            filtered_matches=[match],
            context_chunks=["chunk1", "chunk2"],
            answer="Test answer",
        )
        assert result.query == "test query"
        assert len(result.matches) == 1
        assert len(result.filtered_matches) == 1
        assert len(result.context_chunks) == 2
        assert result.answer == "Test answer"


class TestMatchProcessingHelpers:
    """Tests for match processing helper functions."""

    def test_get_source_with_source_key(self):
        """_get_source should extract 'source' key."""
        meta = {"source": "literature"}
        assert _get_source(meta) == "literature"

    def test_get_source_with_source_type_key(self):
        """_get_source should fallback to 'source_type' key."""
        meta = {"source_type": "omics"}
        assert _get_source(meta) == "omics"

    def test_get_source_unknown_fallback(self):
        """_get_source should return 'Unknown' when no source found."""
        meta = {}
        assert _get_source(meta) == "Unknown"

    def test_format_snippet_truncates(self):
        """_format_snippet should truncate to 200 chars."""
        long_text = "x" * 300
        meta = {"snippet": long_text}
        result = _format_snippet(meta)
        assert len(result) == 200

    def test_format_snippet_removes_newlines(self):
        """_format_snippet should replace newlines with spaces."""
        meta = {"snippet": "line1\nline2\nline3"}
        result = _format_snippet(meta)
        assert "\n" not in result
        assert "line1 line2 line3" == result

    def test_format_snippet_empty(self):
        """_format_snippet should handle empty/missing snippet."""
        meta = {}
        result = _format_snippet(meta)
        assert result == ""

    def test_get_tags_with_tags_key(self):
        """_get_tags should extract 'tags' key."""
        meta = {"tags": ["tag1", "tag2"]}
        assert _get_tags(meta) == ["tag1", "tag2"]

    def test_get_tags_with_zotero_tags_key(self):
        """_get_tags should fallback to 'zotero_tags' key."""
        meta = {"zotero_tags": ["zotero1", "zotero2"]}
        assert _get_tags(meta) == ["zotero1", "zotero2"]

    def test_get_tags_string_to_list(self):
        """_get_tags should convert string to list."""
        meta = {"tags": "single-tag"}
        assert _get_tags(meta) == ["single-tag"]

    def test_get_tags_empty(self):
        """_get_tags should return empty list when no tags."""
        meta = {}
        assert _get_tags(meta) == []


class TestSummarizeMatch:
    """Tests for summarize_match function."""

    def test_summarize_match_with_object(self):
        """summarize_match should handle object-style matches."""
        mock_match = MagicMock()
        mock_match.id = "match-123"
        mock_match.score = 0.85
        mock_match.metadata = {
            "source": "literature",
            "title": "Test Title",
            "snippet": "Test snippet text",
            "tags": ["ALS"],
        }

        result = summarize_match(mock_match)

        assert result.id == "match-123"
        assert result.score == 0.85
        assert result.source == "literature"
        assert result.title == "Test Title"
        assert "Test snippet" in result.snippet
        assert result.tags == ["ALS"]

    def test_summarize_match_with_empty_metadata(self):
        """summarize_match should handle matches with empty metadata."""
        # Use a simple class that mimics Pinecone match objects
        class SimpleMatch:
            def __init__(self):
                self.id = "empty-meta"
                self.score = 0.75
                self.metadata = {}
            
            def get(self, key, default=None):
                # Fallback for dict-style access
                return getattr(self, key, default)

        result = summarize_match(SimpleMatch())

        assert result.id == "empty-meta"
        assert result.score == 0.75
        assert result.source == "Unknown"
        assert result.title == "(untitled)"

    def test_summarize_match_missing_title(self):
        """summarize_match should use '(untitled)' for missing title."""
        mock_match = MagicMock()
        mock_match.id = "no-title"
        mock_match.score = 0.5
        mock_match.metadata = {"source": "literature"}  # has source but no title

        result = summarize_match(mock_match)

        assert result.title == "(untitled)"
        assert result.source == "literature"


class TestFilterMatches:
    """Tests for filter_matches function."""

    @pytest.fixture
    def sample_matches(self):
        """Create sample MatchSummary objects for testing."""
        return [
            MatchSummary(
                id="m1",
                score=0.9,
                source="lit",
                title="Paper 1",
                snippet="About ALS research",
                tags=["ALS", "biomarker"],
                metadata={},
            ),
            MatchSummary(
                id="m2",
                score=0.8,
                source="omics",
                title="Dataset 2",
                snippet="Proteomics data",
                tags=["proteomics", "ALS"],
                metadata={},
            ),
            MatchSummary(
                id="m3",
                score=0.7,
                source="lit",
                title="Paper 3",
                snippet="Cancer research",
                tags=["cancer", "oncology"],
                metadata={},
            ),
        ]

    def test_filter_matches_no_filter(self, sample_matches):
        """filter_matches with no tag should return all matches."""
        result = filter_matches(sample_matches)
        assert len(result) == 3

    def test_filter_matches_by_tag(self, sample_matches):
        """filter_matches should filter by tag substring."""
        result = filter_matches(sample_matches, tag="ALS")
        assert len(result) == 2
        assert all("ALS" in m.tags for m in result)

    def test_filter_matches_case_insensitive(self, sample_matches):
        """filter_matches should be case-insensitive."""
        result = filter_matches(sample_matches, tag="als")
        assert len(result) == 2

    def test_filter_matches_no_match(self, sample_matches):
        """filter_matches should return empty list when no matches."""
        result = filter_matches(sample_matches, tag="nonexistent")
        assert len(result) == 0


class TestQueryRAG:
    """Tests for the main query_rag function."""

    @patch("amprenta_rag.query.rag.query.collect_hybrid_chunks")
    @patch("amprenta_rag.query.rag.query.synthesize_answer")
    @patch("amprenta_rag.query.rag.query.query_pinecone")
    def test_query_rag_basic(self, mock_pinecone, mock_synthesize, mock_chunks):
        """query_rag should orchestrate pinecone query and synthesis."""
        from amprenta_rag.query.rag.query import query_rag

        # Setup mocks - query_pinecone returns list of match objects
        mock_match = MagicMock()
        mock_match.id = "p1"
        mock_match.score = 0.9
        mock_match.metadata = {
            "source": "lit",
            "title": "Test",
            "snippet": "Test content",
        }
        mock_pinecone.return_value = [mock_match]
        mock_chunks.return_value = ["chunk content"]
        mock_synthesize.return_value = "Synthesized answer"

        result = query_rag("test query", top_k=5)

        assert isinstance(result, RAGQueryResult)
        assert result.query == "test query"
        assert len(result.matches) == 1
        mock_pinecone.assert_called_once()

    @patch("amprenta_rag.query.rag.query.query_pinecone")
    def test_query_rag_empty_results(self, mock_pinecone):
        """query_rag should handle empty results gracefully."""
        from amprenta_rag.query.rag.query import query_rag

        mock_pinecone.return_value = []  # Empty list, not dict

        result = query_rag("no results query", top_k=5)

        assert isinstance(result, RAGQueryResult)
        assert len(result.matches) == 0
        assert "No matches" in result.answer


class TestChunkCollection:
    """Tests for chunk collection functionality."""

    def test_collect_chunks_basic(self):
        """collect_chunks should extract snippets from matches."""
        from amprenta_rag.query.rag.chunk_collection import collect_chunks

        matches = [
            MatchSummary(
                id="c1",
                score=0.9,
                source="lit",
                title="T1",
                snippet="Chunk 1 content",
                tags=[],
                metadata={},
            ),
            MatchSummary(
                id="c2",
                score=0.8,
                source="lit",
                title="T2",
                snippet="Chunk 2 content",
                tags=[],
                metadata={},
            ),
        ]

        chunks = collect_chunks(matches)

        assert len(chunks) == 2
        assert "Chunk 1 content" in chunks[0]
        assert "Chunk 2 content" in chunks[1]

