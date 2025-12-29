"""
Supplementary file parser for scientific publications.

Auto-detects table schemas (gene expression, compound activity, proteomics)
and extracts structured data from Excel and CSV files.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from pydantic import BaseModel, Field

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class DetectedSchema(BaseModel):
    """Detected schema type and confidence."""

    schema_type: str = Field(..., description="gene_expression, compound_activity, proteomics, or generic")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Confidence in schema detection")
    detected_columns: Dict[str, str] = Field(default_factory=dict, description="Mapping of detected column roles")


class ExtractedTable(BaseModel):
    """Extracted table data with schema information."""

    sheet_name: Optional[str] = Field(None, description="Sheet name (for Excel files)")
    detected_schema: DetectedSchema
    row_count: int
    columns: List[str]
    data: List[Dict[str, Any]] = Field(default_factory=list, description="First 100 rows as dict")
    full_data_truncated: bool = Field(default=False, description="Whether data was truncated")


class SupplementaryExtractionResult(BaseModel):
    """Result of supplementary file extraction."""

    file_type: str = Field(..., description="csv or excel")
    tables: List[ExtractedTable] = Field(default_factory=list)
    extraction_notes: List[str] = Field(default_factory=list, description="Warnings and notes")


def infer_column_type(column_name: str, values: pd.Series) -> str:
    """
    Infer the type of a column based on name and values.

    Args:
        column_name: Name of the column
        values: Pandas Series of column values

    Returns:
        Column type: gene, protein, p_value, fold_change, ic50, sample, treatment, or generic
    """
    col_lower = column_name.lower()
    
    # Gene/protein identifiers
    if any(kw in col_lower for kw in ["gene", "symbol", "ensembl", "entrez"]):
        return "gene"
    if any(kw in col_lower for kw in ["protein", "uniprot"]):
        return "protein"
    
    # Statistical values
    if any(kw in col_lower for kw in ["pval", "p.val", "p_val", "padj", "fdr", "qval"]):
        return "p_value"
    if any(kw in col_lower for kw in ["log2fc", "logfc", "fold", "fc"]):
        return "fold_change"
    
    # Compound activity
    if any(kw in col_lower for kw in ["ic50", "ec50", "ki", "kd"]):
        return "ic50"
    if "compound" in col_lower or "drug" in col_lower or "smiles" in col_lower:
        return "compound"
    if "target" in col_lower:
        return "target"
    
    # Sample/condition identifiers
    if any(kw in col_lower for kw in ["sample", "replicate", "condition", "treatment"]):
        return "sample"
    
    # Check value patterns if name is ambiguous
    if pd.api.types.is_numeric_dtype(values):
        # Check if values look like p-values (0-1 range)
        try:
            values_numeric = pd.to_numeric(values, errors='coerce').dropna()
            if len(values_numeric) > 0:
                if values_numeric.min() >= 0 and values_numeric.max() <= 1.0:
                    return "p_value"
        except Exception:
            pass
    
    return "generic"


def detect_schema(df: pd.DataFrame) -> DetectedSchema:
    """
    Auto-detect table schema type based on columns.

    Args:
        df: Pandas DataFrame

    Returns:
        DetectedSchema with type, confidence, and column mappings
    """
    if df.empty:
        return DetectedSchema(schema_type="generic", confidence=0.0, detected_columns={})
    
    # Infer column types
    column_types = {}
    for col in df.columns:
        col_type = infer_column_type(col, df[col])
        column_types[col] = col_type
    
    # Count column type occurrences
    type_counts = {}
    for col_type in column_types.values():
        type_counts[col_type] = type_counts.get(col_type, 0) + 1
    
    # Detect schema based on column patterns
    schema_type = "generic"
    confidence = 0.5
    
    # Gene expression: has gene + p_value + fold_change
    if "gene" in type_counts and "p_value" in type_counts and "fold_change" in type_counts:
        schema_type = "gene_expression"
        confidence = 0.9
    
    # Compound activity: has compound + ic50/target
    elif "compound" in type_counts and ("ic50" in type_counts or "target" in type_counts):
        schema_type = "compound_activity"
        confidence = 0.85
    
    # Proteomics: has protein + (p_value or fold_change or sample)
    elif "protein" in type_counts and (
        "p_value" in type_counts or "fold_change" in type_counts or "sample" in type_counts
    ):
        schema_type = "proteomics"
        confidence = 0.85
    
    # Generic with some structure
    elif len([t for t in type_counts if t != "generic"]) >= 2:
        confidence = 0.6
    
    logger.info(
        "[SUPPLEMENTARY_PARSER] Detected schema: %s (confidence=%.2f)",
        schema_type,
        confidence,
    )
    
    return DetectedSchema(
        schema_type=schema_type,
        confidence=confidence,
        detected_columns=column_types,
    )


def parse_excel(file_path: str) -> SupplementaryExtractionResult:
    """
    Parse Excel file and extract structured data.

    Args:
        file_path: Path to Excel file

    Returns:
        SupplementaryExtractionResult with extracted tables
    """
    file_path_obj = Path(file_path)
    
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Excel file not found: {file_path}")
    
    tables = []
    notes = []
    
    try:
        # Read all sheets
        excel_file = pd.ExcelFile(file_path)
        
        for sheet_name in excel_file.sheet_names:
            try:
                df = pd.read_excel(file_path, sheet_name=sheet_name)
                
                if df.empty:
                    notes.append(f"Sheet '{sheet_name}' is empty")
                    continue
                
                # Detect schema
                schema = detect_schema(df)
                
                # Extract data (limit to first 100 rows)
                data_rows = df.head(100).to_dict(orient="records")
                
                table = ExtractedTable(
                    sheet_name=sheet_name,
                    detected_schema=schema,
                    row_count=len(df),
                    columns=list(df.columns),
                    data=data_rows,
                    full_data_truncated=(len(df) > 100),
                )
                
                tables.append(table)
                
                logger.info(
                    "[SUPPLEMENTARY_PARSER] Extracted sheet '%s': %d rows, schema=%s",
                    sheet_name,
                    len(df),
                    schema.schema_type,
                )
            
            except Exception as e:
                notes.append(f"Failed to parse sheet '{sheet_name}': {str(e)}")
                logger.warning("[SUPPLEMENTARY_PARSER] Sheet parsing failed: %r", e)
        
        return SupplementaryExtractionResult(
            file_type="excel",
            tables=tables,
            extraction_notes=notes,
        )
    
    except Exception as e:
        logger.error("[SUPPLEMENTARY_PARSER] Excel parsing failed: %r", e)
        return SupplementaryExtractionResult(
            file_type="excel",
            tables=[],
            extraction_notes=[f"Excel parsing failed: {str(e)}"],
        )


def parse_csv(file_path: str) -> SupplementaryExtractionResult:
    """
    Parse CSV file and extract structured data.

    Args:
        file_path: Path to CSV file

    Returns:
        SupplementaryExtractionResult with extracted table
    """
    file_path_obj = Path(file_path)
    
    if not file_path_obj.exists():
        raise FileNotFoundError(f"CSV file not found: {file_path}")
    
    notes = []
    
    try:
        # Auto-detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
        
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
        
        if df.empty:
            notes.append("CSV file is empty")
            return SupplementaryExtractionResult(
                file_type="csv",
                tables=[],
                extraction_notes=notes,
            )
        
        # Detect schema
        schema = detect_schema(df)
        
        # Extract data (limit to first 100 rows)
        data_rows = df.head(100).to_dict(orient="records")
        
        table = ExtractedTable(
            sheet_name=None,
            detected_schema=schema,
            row_count=len(df),
            columns=list(df.columns),
            data=data_rows,
            full_data_truncated=(len(df) > 100),
        )
        
        logger.info(
            "[SUPPLEMENTARY_PARSER] Extracted CSV: %d rows, schema=%s",
            len(df),
            schema.schema_type,
        )
        
        return SupplementaryExtractionResult(
            file_type="csv",
            tables=[table],
            extraction_notes=notes,
        )
    
    except Exception as e:
        logger.error("[SUPPLEMENTARY_PARSER] CSV parsing failed: %r", e)
        return SupplementaryExtractionResult(
            file_type="csv",
            tables=[],
            extraction_notes=[f"CSV parsing failed: {str(e)}"],
        )


def parse_supplementary_file(file_path: str) -> SupplementaryExtractionResult:
    """
    Parse supplementary file (auto-detect Excel or CSV).

    Args:
        file_path: Path to supplementary file

    Returns:
        SupplementaryExtractionResult with extracted data
    """
    file_path_obj = Path(file_path)
    suffix = file_path_obj.suffix.lower()
    
    if suffix in [".xlsx", ".xls"]:
        return parse_excel(file_path)
    elif suffix in [".csv", ".tsv", ".txt"]:
        return parse_csv(file_path)
    else:
        logger.error("[SUPPLEMENTARY_PARSER] Unsupported file type: %s", suffix)
        return SupplementaryExtractionResult(
            file_type="unknown",
            tables=[],
            extraction_notes=[f"Unsupported file type: {suffix}"],
        )

