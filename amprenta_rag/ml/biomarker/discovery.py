"""Orchestration for biomarker discovery methods and consensus ranking."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.database.models import Experiment as ExperimentModel
from amprenta_rag.database.session import db_session
from amprenta_rag.ml.biomarker.datasets import load_biomarker_dataset
from amprenta_rag.ml.biomarker.importance import CVFeatureImportance
from amprenta_rag.ml.biomarker.stability import StabilitySelector
from amprenta_rag.ml.biomarker.statistical import fdr_correction, t_test_selection


def _as_sample_list(exp: ExperimentModel, group: Any) -> List[str]:
    if isinstance(group, (list, tuple)):
        return [str(s) for s in group]
    if isinstance(group, dict):
        for k in ("samples", "sample_ids", "ids"):
            if k in group and isinstance(group[k], (list, tuple)):
                return [str(s) for s in group[k]]
        return []
    if isinstance(group, str):
        sg = getattr(exp, "sample_groups", None) or {}
        if isinstance(sg, dict):
            val = sg.get(group)
            if isinstance(val, list):
                return [str(s) for s in val]
    return []


def _pick_omics_type(exp: ExperimentModel) -> str:
    dsets = getattr(exp, "datasets", []) or []
    for ds in dsets:
        ot = getattr(ds, "omics_type", None)
        if ot:
            return str(ot)
    return "omics"


@dataclass
class BiomarkerDiscoveryResult:
    consensus_ranking: List[Dict[str, Any]]
    method_results: Dict[str, Any]


class BiomarkerDiscoveryService:
    def discover(
        self,
        experiment_id: UUID | str,
        group1: Any,
        group2: Any,
        methods: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        exp_id = UUID(str(experiment_id)) if not isinstance(experiment_id, UUID) else experiment_id
        methods = methods or ["statistical", "stability", "importance"]
        methods = [str(m).lower() for m in methods]

        with db_session() as db:
            exp = db.query(ExperimentModel).filter(ExperimentModel.id == exp_id).first()
            if not exp:
                raise ValueError("Experiment not found")
            g1 = _as_sample_list(exp, group1)
            g2 = _as_sample_list(exp, group2)
            omics_type = _pick_omics_type(exp)

        X, y, feature_names = load_biomarker_dataset(
            experiment_id=exp_id,
            group1_samples=g1,
            group2_samples=g2,
            omics_type=omics_type,
            min_coverage=0.5,
        )

        method_rankings: Dict[str, List[str]] = {}
        method_results: Dict[str, Any] = {}

        if "statistical" in methods:
            stats_rows = t_test_selection(X, y, feature_names=feature_names)
            pvals = [r[2] for r in stats_rows]
            try:
                _, p_adj = fdr_correction(pvals)
            except Exception:
                p_adj = pvals
            out = []
            for (fname, tstat, pval), padj in zip(stats_rows, p_adj):
                out.append({"feature": fname, "t_stat": float(tstat), "p_value": float(pval), "p_adj": float(padj)})
            method_results["statistical"] = out
            method_rankings["statistical"] = [r["feature"] for r in out]

        if "stability" in methods:
            sel = StabilitySelector(n_bootstrap=50, threshold=0.0)
            sel.fit(X, y, feature_names=feature_names)
            ranked = sel.get_ranked_features()
            method_results["stability"] = [
                {"feature": f, "frequency": float(freq), "mean_coef": float(coef)} for (f, freq, coef) in ranked
            ]
            method_rankings["stability"] = [f for (f, _, _) in ranked]

        if "importance" in methods:
            imp = CVFeatureImportance(n_folds=5, n_estimators=100)
            imp.fit(X, y, feature_names=feature_names)
            ranked2 = imp.get_ranked_features()
            method_results["importance"] = [
                {"feature": f, "mean_importance": float(m), "std_importance": float(s)} for (f, m, s) in ranked2
            ]
            method_rankings["importance"] = [f for (f, _, _) in ranked2]

        consensus = self._consensus_ranking(method_rankings)
        return {"consensus_ranking": consensus, "method_results": method_results}

    def _consensus_ranking(self, results_dict: Dict[str, List[str]]) -> List[Dict[str, Any]]:
        # Average rank across provided methods (lower is better)
        rank_sums: Dict[str, float] = {}
        rank_counts: Dict[str, int] = {}

        for _, ranked in results_dict.items():
            for i, feat in enumerate(ranked):
                if not feat:
                    continue
                rank_sums[feat] = rank_sums.get(feat, 0.0) + float(i + 1)
                rank_counts[feat] = rank_counts.get(feat, 0) + 1

        rows: List[Dict[str, Any]] = []
        for feat, total in rank_sums.items():
            c = rank_counts.get(feat, 1)
            rows.append({"feature": feat, "avg_rank": float(total / float(c)), "methods": int(c)})

        rows.sort(key=lambda r: r["avg_rank"])
        return rows


__all__ = ["BiomarkerDiscoveryService", "BiomarkerDiscoveryResult"]


