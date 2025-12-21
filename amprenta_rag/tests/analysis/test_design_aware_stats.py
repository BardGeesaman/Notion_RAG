from __future__ import annotations

import numpy as np
import pandas as pd

from amprenta_rag.analysis import design_aware_stats as das


def test_run_design_aware_analysis_case_control(monkeypatch) -> None:
    df = pd.DataFrame([[1, 2], [3, 4]], index=["f1", "f2"], columns=["c1", "c2"])

    monkeypatch.setattr(das, "ttest_independent", lambda a, b: {"pvalue": 0.01})
    monkeypatch.setattr(das, "adjust_pvalues", lambda pvals, method="fdr_bh": ([0.02] * len(pvals), [True] * len(pvals)))

    result = das.run_design_aware_analysis(
        features_df=df,
        sample_groups={"control": ["c1"], "case": ["c2"]},
        design_type="case_control",
    )

    assert result["test_type"] == "t-test"
    assert result["p_values"]["f1"] == 0.01
    assert result["significant_features"] == ["f1", "f2"]


def test_run_design_aware_analysis_dose_response(monkeypatch) -> None:
    df = pd.DataFrame([[1, 2, 3]], index=["f1"], columns=["d1", "d2", "d3"])

    monkeypatch.setattr(das, "pearson_corr", lambda doses, means: {"pvalue": 0.05})
    monkeypatch.setattr(das, "adjust_pvalues", lambda pvals, method="fdr_bh": ([0.05], [False]))

    result = das.run_design_aware_analysis(
        features_df=df,
        sample_groups={"d1": ["d1"], "d2": ["d2"], "d3": ["d3"]},
        design_type="dose_response",
    )

    assert result["test_type"] == "correlation"
    assert result["p_values"]["f1"] == 0.05
    assert result["significant_features"] == []


def test_run_design_aware_analysis_observational(monkeypatch) -> None:
    df = pd.DataFrame([[np.nan, np.nan]], index=["f1"], columns=["s1", "s2"])
    monkeypatch.setattr(das, "adjust_pvalues", lambda pvals, method="fdr_bh": ([np.nan], [False]))

    result = das.run_design_aware_analysis(
        features_df=df,
        sample_groups={},
        design_type="observational",
    )

    assert result["test_type"] == "none"
    assert np.isnan(result["p_values"]["f1"])

