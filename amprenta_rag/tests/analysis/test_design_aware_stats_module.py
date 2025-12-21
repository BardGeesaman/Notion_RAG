from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd

from amprenta_rag.analysis import design_aware_stats as das


def _make_df():
    # Two features x three samples
    data = [[1, 2, 3], [2, 2, 2]]
    return pd.DataFrame(data, index=["f1", "f2"], columns=["s1", "s2", "s3"])


def test_extract_group_values_filters_missing():
    df = pd.DataFrame({"s1": [1, 2]})
    arr = das._extract_group_values(df, ["s2"])
    assert arr.size == 0


def test_ttest_by_feature_returns_nan_when_missing():
    df = _make_df()
    pvals = das._ttest_by_feature(df, ["s1"], ["s2"])
    assert all(np.isnan(v) for v in pvals.values())


def test_anova_by_feature_basic(monkeypatch):
    df = _make_df()
    # Provide a fake anova returning p=0.05
    monkeypatch.setattr(das, "anova_oneway", lambda *a, **k: {"pvalue": 0.05})
    pvals = das._anova_by_feature(df, {"g1": ["s1", "s2"], "g2": ["s3", "s4"]})
    assert "f1" in pvals and isinstance(pvals["f1"], float)


def test_dose_response_parses_and_falls_back(monkeypatch):
    df = _make_df()
    groups: Dict[str, List[str]] = {"1uM": ["s1"], "NaN": ["s2"], "3": ["s3"]}
    monkeypatch.setattr(das, "pearson_corr", lambda doses, means: {"pvalue": 0.01})
    pvals = das._dose_response_by_feature(df, groups)
    assert all(isinstance(v, float) for v in pvals.values())


def test_run_design_aware_analysis_paths(monkeypatch):
    df = _make_df()
    monkeypatch.setattr(das, "_ttest_by_feature", lambda df, a, b: {"f1": 0.04})
    monkeypatch.setattr(das, "_anova_by_feature", lambda df, g: {"f1": 0.02})
    monkeypatch.setattr(das, "_dose_response_by_feature", lambda df, g: {"f1": 0.03})
    monkeypatch.setattr(das, "adjust_pvalues", lambda pv, method="fdr_bh": ([0.05], [True]))

    out_case = das.run_design_aware_analysis(df, {"control": ["s1"], "case": ["s2"]}, "case_control")
    assert out_case["test_type"] == "t-test"

    out_time = das.run_design_aware_analysis(df, {"timepoints": {"t0": ["s1"], "t1": ["s2"]}}, "time_course")
    assert out_time["test_type"] == "anova"

    out_intervention = das.run_design_aware_analysis(df, {"treated": ["s1"], "untreated": ["s2"]}, "intervention")
    assert out_intervention["test_type"] == "t-test"

    out_dose = das.run_design_aware_analysis(df, {"1": ["s1"], "2": ["s2"]}, "dose_response")
    assert out_dose["test_type"] == "correlation"

    out_multi = das.run_design_aware_analysis(df, {"g1": ["s1"], "g2": ["s2"]}, "multi_factorial")
    assert out_multi["test_type"].startswith("two-way")

    out_obs = das.run_design_aware_analysis(df, {"g1": ["s1"]}, "observational")
    assert out_obs["test_type"] == "none"

