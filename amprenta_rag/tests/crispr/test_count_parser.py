from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from amprenta_rag.crispr.count_parser import parse_count_matrix, validate_count_matrix


def test_validate_count_matrix_requires_sgrna_and_gene():
    df = pd.DataFrame({"sgRNA": ["g1"], "sampleA": [10]})
    with pytest.raises(ValueError, match="Gene"):
        validate_count_matrix(df)

    df2 = pd.DataFrame({"Gene": ["TP53"], "sampleA": [10]})
    with pytest.raises(ValueError, match="sgRNA"):
        validate_count_matrix(df2)


def test_validate_count_matrix_requires_sample_columns():
    df = pd.DataFrame({"sgRNA": ["g1"], "Gene": ["TP53"]})
    with pytest.raises(ValueError, match="sample"):
        validate_count_matrix(df)


def test_parse_count_matrix_reads_and_validates(tmp_path: Path):
    p = tmp_path / "counts.tsv"
    p.write_text("sgRNA\tGene\tcontrol\ttreatment\nsg1\tTP53\t10\t20\n", encoding="utf-8")
    df = parse_count_matrix(str(p))
    assert list(df.columns) == ["sgRNA", "Gene", "control", "treatment"]
    assert df.shape == (1, 4)


