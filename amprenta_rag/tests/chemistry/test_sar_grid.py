"""Tests for SAR grid functionality."""

from __future__ import annotations

import pytest
import pandas as pd

from amprenta_rag.chemistry.sar_analysis import build_sar_matrix


def test_build_sar_matrix_basic() -> None:
    """Test basic SAR matrix construction."""
    compounds = [
        {"R1": "H", "R2": "Me", "ic50": 100.0},
        {"R1": "H", "R2": "Et", "ic50": 50.0},
        {"R1": "Me", "R2": "Me", "ic50": 200.0},
        {"R1": "Me", "R2": "Et", "ic50": 75.0},
    ]
    
    matrix = build_sar_matrix(compounds, x_axis="R1", y_axis="R2", activity_col="ic50")
    
    # Check shape: 2 R2 values (Et, Me) x 2 R1 values (H, Me)
    assert matrix.shape == (2, 2)
    
    # Check values are correct
    assert matrix.loc["Me", "H"] == 100.0  # R1=H, R2=Me
    assert matrix.loc["Et", "H"] == 50.0   # R1=H, R2=Et
    assert matrix.loc["Me", "Me"] == 200.0 # R1=Me, R2=Me
    assert matrix.loc["Et", "Me"] == 75.0  # R1=Me, R2=Et
    
    # Check index and column names
    assert matrix.index.name == "R2"
    assert matrix.columns.name == "R1"


def test_build_sar_matrix_aggregation() -> None:
    """Test that duplicate R-group combinations are averaged."""
    compounds = [
        {"R1": "H", "R2": "Me", "ic50": 100.0},
        {"R1": "H", "R2": "Me", "ic50": 200.0},  # Duplicate combination
        {"R1": "Me", "R2": "Et", "ic50": 50.0},
        {"R1": "Me", "R2": "Et", "ic50": 150.0}, # Duplicate combination
    ]
    
    matrix = build_sar_matrix(compounds, x_axis="R1", y_axis="R2", activity_col="ic50")
    
    # Check aggregated values (mean)
    assert matrix.loc["Me", "H"] == 150.0  # (100 + 200) / 2
    assert matrix.loc["Et", "Me"] == 100.0 # (50 + 150) / 2


def test_build_sar_matrix_missing_rgroup() -> None:
    """Test handling of compounds with missing R-groups."""
    compounds = [
        {"R1": "H", "R2": "Me", "ic50": 100.0},
        {"R1": "H", "R2": None, "ic50": 50.0},    # Missing R2
        {"R1": None, "R2": "Et", "ic50": 75.0},   # Missing R1
        {"R1": "Me", "R2": "Et", "ic50": 200.0},
    ]
    
    matrix = build_sar_matrix(compounds, x_axis="R1", y_axis="R2", activity_col="ic50")
    
    # Should only include compounds with both R1 and R2 defined
    assert matrix.shape == (2, 2)  # Only 2 complete compounds
    assert matrix.loc["Me", "H"] == 100.0
    assert matrix.loc["Et", "Me"] == 200.0


def test_build_sar_matrix_missing_activity() -> None:
    """Test handling of compounds with missing activity data."""
    compounds = [
        {"R1": "H", "R2": "Me", "ic50": 100.0},
        {"R1": "H", "R2": "Et", "ic50": None},    # Missing activity
        {"R1": "Me", "R2": "Me", "ic50": 200.0},
    ]
    
    matrix = build_sar_matrix(compounds, x_axis="R1", y_axis="R2", activity_col="ic50")
    
    # Should only include compounds with activity data
    assert matrix.shape == (1, 2)  # Only 2 compounds with activity
    assert matrix.loc["Me", "H"] == 100.0
    assert matrix.loc["Me", "Me"] == 200.0


def test_build_sar_matrix_empty_input() -> None:
    """Test handling of empty compound list."""
    matrix = build_sar_matrix([], x_axis="R1", y_axis="R2", activity_col="ic50")
    
    assert matrix.empty
    assert isinstance(matrix, pd.DataFrame)


def test_build_sar_matrix_missing_columns() -> None:
    """Test handling of missing required columns."""
    compounds = [
        {"R1": "H", "ic50": 100.0},  # Missing R2
        {"R2": "Me", "ic50": 50.0},  # Missing R1
    ]
    
    matrix = build_sar_matrix(compounds, x_axis="R1", y_axis="R2", activity_col="ic50")
    
    # Should return empty DataFrame when required columns are missing
    assert matrix.empty
    assert isinstance(matrix, pd.DataFrame)


def test_build_sar_matrix_custom_axes() -> None:
    """Test SAR matrix with different axis configurations."""
    compounds = [
        {"R1": "H", "R3": "Ph", "pic50": 7.0},
        {"R1": "Me", "R3": "Ph", "pic50": 8.0},
        {"R1": "H", "R3": "Py", "pic50": 6.5},
        {"R1": "Me", "R3": "Py", "pic50": 7.5},
    ]
    
    matrix = build_sar_matrix(compounds, x_axis="R3", y_axis="R1", activity_col="pic50")
    
    # Check shape and axis names
    assert matrix.shape == (2, 2)  # 2 R1 values x 2 R3 values
    assert matrix.index.name == "R1"
    assert matrix.columns.name == "R3"
    
    # Check values
    assert matrix.loc["H", "Ph"] == 7.0
    assert matrix.loc["Me", "Ph"] == 8.0
    assert matrix.loc["H", "Py"] == 6.5
    assert matrix.loc["Me", "Py"] == 7.5


def test_build_sar_matrix_single_rgroup() -> None:
    """Test handling when only one R-group position varies."""
    compounds = [
        {"R1": "H", "R2": "Me", "ic50": 100.0},
        {"R1": "Et", "R2": "Me", "ic50": 50.0},
        {"R1": "Ph", "R2": "Me", "ic50": 25.0},
    ]
    
    matrix = build_sar_matrix(compounds, x_axis="R1", y_axis="R2", activity_col="ic50")
    
    # Should have 1 row (R2=Me) and 3 columns (different R1 values)
    assert matrix.shape == (1, 3)
    assert list(matrix.columns) == ["Et", "H", "Ph"]  # Sorted alphabetically
    assert matrix.loc["Me", "H"] == 100.0
    assert matrix.loc["Me", "Et"] == 50.0
    assert matrix.loc["Me", "Ph"] == 25.0


# Import and test the ic50_to_pic50 function from the dashboard module
def test_ic50_to_pic50() -> None:
    """Test IC50 to pIC50 conversion function."""
    from scripts.dashboard.pages.chemistry.sar_analysis import ic50_to_pic50
    
    # Test normal conversions
    assert ic50_to_pic50(100.0) == pytest.approx(7.0, rel=1e-9)  # -log10(100e-9) = 7
    assert ic50_to_pic50(10.0) == pytest.approx(8.0, rel=1e-9)   # -log10(10e-9) = 8
    assert ic50_to_pic50(1.0) == pytest.approx(9.0, rel=1e-9)    # -log10(1e-9) = 9
    
    # Test edge cases
    assert ic50_to_pic50(None) is None
    assert ic50_to_pic50(0.0) is None
    assert ic50_to_pic50(-1.0) is None
    
    # Test very small values
    assert ic50_to_pic50(0.001) == pytest.approx(12.0, rel=1e-9)  # -log10(0.001e-9) = 12


def test_ic50_to_pic50_precision() -> None:
    """Test pIC50 calculation precision with various IC50 values."""
    from scripts.dashboard.pages.chemistry.sar_analysis import ic50_to_pic50
    
    # Test cases with known pIC50 values
    test_cases = [
        (1000.0, 6.0),    # 1 μM
        (100.0, 7.0),     # 100 nM
        (10.0, 8.0),      # 10 nM
        (1.0, 9.0),       # 1 nM
        (0.1, 10.0),      # 100 pM
    ]
    
    for ic50_nm, expected_pic50 in test_cases:
        result = ic50_to_pic50(ic50_nm)
        assert result == pytest.approx(expected_pic50, rel=1e-9), f"IC50 {ic50_nm} nM should give pIC50 {expected_pic50}"
