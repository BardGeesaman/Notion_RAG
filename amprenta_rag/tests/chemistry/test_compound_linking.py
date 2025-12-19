"""
Tests for compound_linking module (PostgreSQL version).
"""
from unittest.mock import MagicMock, patch
from uuid import uuid4


from amprenta_rag.chemistry.compound_linking import (
    link_compound_to_signature,
    link_compound_to_program,
    get_compounds_for_signature,
    get_compounds_for_program,
)


class TestLinkCompoundToSignature:
    @patch("amprenta_rag.chemistry.compound_linking.get_db")
    def test_link_success(self, mock_get_db):
        mock_compound = MagicMock()
        mock_compound.external_ids = {}

        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_compound
        mock_gen = MagicMock()
        mock_gen.__next__ = MagicMock(return_value=mock_db)
        mock_gen.close = MagicMock()
        mock_get_db.return_value = mock_gen

        result = link_compound_to_signature("C1", "S1")
        assert result is True
        assert "S1" in mock_compound.external_ids.get("signatures", [])
        mock_db.commit.assert_called_once()

    @patch("amprenta_rag.chemistry.compound_linking.get_db")
    def test_compound_not_found(self, mock_get_db):
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        mock_gen = MagicMock()
        mock_gen.__next__ = MagicMock(return_value=mock_db)
        mock_gen.close = MagicMock()
        mock_get_db.return_value = mock_gen

        result = link_compound_to_signature("INVALID", "S1")
        assert result is False


class TestLinkCompoundToProgram:
    @patch("amprenta_rag.chemistry.compound_linking.get_db")
    def test_link_success(self, mock_get_db):
        mock_compound = MagicMock()
        mock_compound.programs = []
        mock_program = MagicMock()

        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.side_effect = [mock_compound, mock_program]
        mock_gen = MagicMock()
        mock_gen.__next__ = MagicMock(return_value=mock_db)
        mock_gen.close = MagicMock()
        mock_get_db.return_value = mock_gen

        result = link_compound_to_program("C1", str(uuid4()))
        assert result is True
        mock_db.commit.assert_called_once()

    @patch("amprenta_rag.chemistry.compound_linking.get_db")
    def test_compound_not_found(self, mock_get_db):
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        mock_gen = MagicMock()
        mock_gen.__next__ = MagicMock(return_value=mock_db)
        mock_gen.close = MagicMock()
        mock_get_db.return_value = mock_gen

        result = link_compound_to_program("INVALID", str(uuid4()))
        assert result is False


class TestGetCompoundsForSignature:
    @patch("amprenta_rag.chemistry.compound_linking.get_db")
    def test_returns_matching_compounds(self, mock_get_db):
        mock_c1 = MagicMock()
        mock_c1.compound_id = "C1"
        mock_c1.external_ids = {"signatures": ["S1", "S2"]}

        mock_c2 = MagicMock()
        mock_c2.compound_id = "C2"
        mock_c2.external_ids = {"signatures": ["S3"]}

        mock_db = MagicMock()
        mock_db.query.return_value.all.return_value = [mock_c1, mock_c2]
        mock_gen = MagicMock()
        mock_gen.__next__ = MagicMock(return_value=mock_db)
        mock_gen.close = MagicMock()
        mock_get_db.return_value = mock_gen

        result = get_compounds_for_signature("S1")
        assert result == ["C1"]

    @patch("amprenta_rag.chemistry.compound_linking.get_db")
    def test_returns_empty_when_no_match(self, mock_get_db):
        mock_db = MagicMock()
        mock_db.query.return_value.all.return_value = []
        mock_gen = MagicMock()
        mock_gen.__next__ = MagicMock(return_value=mock_db)
        mock_gen.close = MagicMock()
        mock_get_db.return_value = mock_gen

        result = get_compounds_for_signature("NONEXISTENT")
        assert result == []


class TestGetCompoundsForProgram:
    @patch("amprenta_rag.chemistry.compound_linking.get_db")
    def test_returns_program_compounds(self, mock_get_db):
        mock_c1 = MagicMock()
        mock_c1.compound_id = "C1"
        mock_program = MagicMock()
        mock_program.compounds = [mock_c1]

        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_program
        mock_gen = MagicMock()
        mock_gen.__next__ = MagicMock(return_value=mock_db)
        mock_gen.close = MagicMock()
        mock_get_db.return_value = mock_gen

        result = get_compounds_for_program(str(uuid4()))
        assert result == ["C1"]

    @patch("amprenta_rag.chemistry.compound_linking.get_db")
    def test_program_not_found(self, mock_get_db):
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        mock_gen = MagicMock()
        mock_gen.__next__ = MagicMock(return_value=mock_db)
        mock_gen.close = MagicMock()
        mock_get_db.return_value = mock_gen

        result = get_compounds_for_program(str(uuid4()))
        assert result == []
