from __future__ import annotations

import pandas as pd
import pytest

from amprenta_rag.multi_omics.data_extractor import align_samples


def test_align_samples_row_oriented_reorders_and_renames():
    mat_rna = pd.DataFrame(
        {"RNA_B": [1.0, 2.0], "RNA_A": [3.0, 4.0]},
        index=["G1", "G2"],
    )
    mat_prot = pd.DataFrame(
        {"PROT_A": [10.0, 20.0], "PROT_B": [30.0, 40.0]},
        index=["P1", "P2"],
    )
    mapping = [
        {"sample_id": "S1", "transcriptomics": "RNA_A", "proteomics": "PROT_A"},
        {"sample_id": "S2", "transcriptomics": "RNA_B", "proteomics": "PROT_B"},
    ]
    out = align_samples({"transcriptomics": mat_rna, "proteomics": mat_prot}, mapping)
    assert list(out["transcriptomics"].columns) == ["S1", "S2"]
    assert list(out["proteomics"].columns) == ["S1", "S2"]
    # Values should follow mapping order
    assert out["transcriptomics"].loc["G1", "S1"] == pytest.approx(3.0)
    assert out["transcriptomics"].loc["G1", "S2"] == pytest.approx(1.0)
    assert out["proteomics"].loc["P1", "S1"] == pytest.approx(10.0)
    assert out["proteomics"].loc["P1", "S2"] == pytest.approx(30.0)


def test_align_samples_raises_on_missing_mapped_column():
    mat = pd.DataFrame({"A": [1.0]}, index=["G1"])
    mapping = [{"sample_id": "S1", "transcriptomics": "MISSING"}]
    with pytest.raises(ValueError, match="missing mapped sample columns"):
        align_samples({"transcriptomics": mat}, mapping)


