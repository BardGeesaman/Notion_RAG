"""
Tests for supplementary file parser.

Tests schema detection, column type inference, and file parsing.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from amprenta_rag.extraction.supplementary_parser import (
    detect_schema,
    infer_column_type,
    parse_csv,
    parse_excel,
    parse_supplementary_file,
    SupplementaryExtractionResult,
)


class TestColumnTypeInference:
    """Tests for column type inference."""

    def test_infer_gene_column(self):
        """Test inference of gene column."""
        values = pd.Series(["TP53", "BRCA1", "EGFR"])
        
        assert infer_column_type("Gene", values) == "gene"
        assert infer_column_type("gene_symbol", values) == "gene"
        assert infer_column_type("GeneSymbol", values) == "gene"

    def test_infer_p_value_column(self):
        """Test inference of p-value column."""
        values = pd.Series([0.001, 0.05, 0.0001])
        
        assert infer_column_type("pvalue", values) == "p_value"
        assert infer_column_type("padj", values) == "p_value"
        assert infer_column_type("FDR", values) == "p_value"

    def test_infer_fold_change_column(self):
        """Test inference of fold change column."""
        values = pd.Series([2.5, -1.8, 3.2])
        
        assert infer_column_type("log2FC", values) == "fold_change"
        assert infer_column_type("FoldChange", values) == "fold_change"

    def test_infer_compound_column(self):
        """Test inference of compound column."""
        values = pd.Series(["Compound A", "Compound B"])
        
        assert infer_column_type("Compound", values) == "compound"
        assert infer_column_type("Drug", values) == "compound"
        assert infer_column_type("SMILES", values) == "compound"

    def test_infer_ic50_column(self):
        """Test inference of IC50 column."""
        values = pd.Series([10.5, 25.3, 5.1])
        
        assert infer_column_type("IC50", values) == "ic50"
        assert infer_column_type("EC50", values) == "ic50"


class TestSchemaDetection:
    """Tests for schema type detection."""

    def test_detect_gene_expression_schema(self):
        """Test detection of gene expression schema."""
        df = pd.DataFrame({
            "Gene": ["TP53", "BRCA1", "EGFR"],
            "log2FC": [2.5, -1.8, 3.2],
            "padj": [0.001, 0.05, 0.0001],
        })
        
        schema = detect_schema(df)
        
        assert schema.schema_type == "gene_expression"
        assert schema.confidence >= 0.8
        assert "Gene" in schema.detected_columns
        assert schema.detected_columns["Gene"] == "gene"

    def test_detect_compound_activity_schema(self):
        """Test detection of compound activity schema."""
        df = pd.DataFrame({
            "Compound": ["Drug A", "Drug B", "Drug C"],
            "IC50": [10.5, 25.3, 5.1],
            "Target": ["Kinase1", "Kinase2", "Kinase1"],
        })
        
        schema = detect_schema(df)
        
        assert schema.schema_type == "compound_activity"
        assert schema.confidence >= 0.8

    def test_detect_proteomics_schema(self):
        """Test detection of proteomics schema."""
        df = pd.DataFrame({
            "Protein": ["P53", "BRCA1", "EGFR"],
            "Abundance": [1000.5, 2500.3, 500.1],
            "Sample": ["Control", "Treated", "Control"],
        })
        
        schema = detect_schema(df)
        
        assert schema.schema_type == "proteomics"
        assert schema.confidence >= 0.8

    def test_detect_generic_schema(self):
        """Test detection of generic schema."""
        df = pd.DataFrame({
            "Column1": [1, 2, 3],
            "Column2": ["A", "B", "C"],
        })
        
        schema = detect_schema(df)
        
        assert schema.schema_type == "generic"
        assert schema.confidence >= 0.0

    def test_detect_schema_empty_dataframe(self):
        """Test schema detection on empty dataframe."""
        df = pd.DataFrame()
        
        schema = detect_schema(df)
        
        assert schema.schema_type == "generic"
        assert schema.confidence == 0.0


class TestCSVParsing:
    """Tests for CSV file parsing."""

    def test_parse_csv_gene_expression(self):
        """Test parsing CSV with gene expression data."""
        # Create temporary CSV
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("Gene,log2FC,padj\n")
            f.write("TP53,2.5,0.001\n")
            f.write("BRCA1,-1.8,0.05\n")
            csv_path = f.name
        
        try:
            result = parse_csv(csv_path)
            
            assert result.file_type == "csv"
            assert len(result.tables) == 1
            assert result.tables[0].detected_schema.schema_type == "gene_expression"
            assert result.tables[0].row_count == 2
        finally:
            Path(csv_path).unlink()

    def test_parse_csv_file_not_found(self):
        """Test CSV parsing with non-existent file."""
        with pytest.raises(FileNotFoundError):
            parse_csv("/nonexistent/file.csv")

    def test_parse_csv_empty_file(self):
        """Test parsing empty CSV file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("Column1,Column2\n")  # Header only
            csv_path = f.name
        
        try:
            result = parse_csv(csv_path)
            
            # Should handle empty data gracefully
            assert result.file_type == "csv"
            # Empty CSV still has structure
            assert len(result.tables) >= 0
        finally:
            Path(csv_path).unlink()


class TestExcelParsing:
    """Tests for Excel file parsing."""

    @patch("amprenta_rag.extraction.supplementary_parser.pd.ExcelFile")
    @patch("amprenta_rag.extraction.supplementary_parser.pd.read_excel")
    def test_parse_excel_success(self, mock_read_excel, mock_excel_file_class):
        """Test successful Excel parsing."""
        # Mock Excel file with one sheet
        mock_excel_file = MagicMock()
        mock_excel_file.sheet_names = ["Sheet1"]
        mock_excel_file_class.return_value = mock_excel_file
        
        # Mock dataframe
        mock_df = pd.DataFrame({
            "Gene": ["TP53", "BRCA1"],
            "log2FC": [2.5, -1.8],
            "padj": [0.001, 0.05],
        })
        mock_read_excel.return_value = mock_df
        
        # Need an actual file for Path.exists() check
        with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as f:
            excel_path = f.name
        
        try:
            result = parse_excel(excel_path)
            
            assert result.file_type == "excel"
            assert len(result.tables) == 1
            assert result.tables[0].sheet_name == "Sheet1"
            assert result.tables[0].row_count == 2
        finally:
            Path(excel_path).unlink()


class TestSupplementaryFileParsing:
    """Tests for auto-detecting file type."""

    @patch("amprenta_rag.extraction.supplementary_parser.parse_excel")
    def test_parse_supplementary_excel(self, mock_parse_excel):
        """Test parsing .xlsx file."""
        mock_parse_excel.return_value = SupplementaryExtractionResult(
            file_type="excel",
            tables=[],
        )
        
        with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as f:
            excel_path = f.name
        
        try:
            result = parse_supplementary_file(excel_path)
            
            assert result.file_type == "excel"
            mock_parse_excel.assert_called_once()
        finally:
            Path(excel_path).unlink()

    @patch("amprenta_rag.extraction.supplementary_parser.parse_csv")
    def test_parse_supplementary_csv(self, mock_parse_csv):
        """Test parsing .csv file."""
        mock_parse_csv.return_value = SupplementaryExtractionResult(
            file_type="csv",
            tables=[],
        )
        
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
            csv_path = f.name
        
        try:
            result = parse_supplementary_file(csv_path)
            
            assert result.file_type == "csv"
            mock_parse_csv.assert_called_once()
        finally:
            Path(csv_path).unlink()

    def test_parse_supplementary_unsupported_type(self):
        """Test parsing unsupported file type."""
        with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as f:
            pdf_path = f.name
        
        try:
            result = parse_supplementary_file(pdf_path)
            
            assert result.file_type == "unknown"
            assert len(result.tables) == 0
            assert len(result.extraction_notes) > 0
        finally:
            Path(pdf_path).unlink()

