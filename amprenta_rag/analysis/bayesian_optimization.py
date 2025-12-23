"""Bayesian optimization helpers (BoTorch).

This module provides a small, backend-agnostic recommendation function that
uses a GP surrogate + qExpectedImprovement to rank candidate compounds.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence


def _require_botorch():
    try:
        import torch  # type: ignore
        from botorch.acquisition.monte_carlo import qExpectedImprovement  # type: ignore
        from botorch.fit import fit_gpytorch_mll  # type: ignore
        from botorch.models import SingleTaskGP  # type: ignore
        from botorch.sampling.normal import SobolQMCNormalSampler  # type: ignore
        from gpytorch.mlls import ExactMarginalLogLikelihood  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError(
            "BoTorch is required for Bayesian optimization. "
            "Install with: botorch>=0.9 ax-platform>=0.3 (and torch)."
        ) from e
    return torch, qExpectedImprovement, fit_gpytorch_mll, SingleTaskGP, SobolQMCNormalSampler, ExactMarginalLogLikelihood


def _as_float_matrix(rows: Sequence[Dict[str, Any]], key: str = "features") -> List[List[float]]:
    out: List[List[float]] = []
    for r in rows:
        feats = r.get(key)
        if not isinstance(feats, list) or not feats:
            raise ValueError(f"Each entry must have non-empty list '{key}'")
        out.append([float(x) for x in feats])
    return out


def recommend_next_compounds(
    tested_compounds: List[Dict],  # {"features": [...], "activity": float}
    candidate_pool: List[Dict],  # {"id": str, "features": [...]}
    batch_size: int = 10,
    objectives: List[str] = ["potency"],
) -> List[Dict]:
    """
    BoTorch GP surrogate + qEI acquisition.

    Assumptions:
    - Single objective: maximize activity (potency proxy).
    - Features are numeric vectors with consistent dimensionality.

    Returns:
      ranked list of dicts:
        {"id": str, "acquisition": float, "pred_mean": float, "pred_std": float}
    """
    if not objectives:
        raise ValueError("objectives must be non-empty")
    if objectives != ["potency"]:
        # MVP: only potency supported; expand later to multi-objective.
        raise ValueError("Only objectives=['potency'] is supported in Phase 2 MVP")
    if batch_size <= 0:
        raise ValueError("batch_size must be > 0")
    if not tested_compounds or len(tested_compounds) < 5:
        raise ValueError("Need at least 5 tested_compounds to fit a GP")
    if not candidate_pool:
        return []

    torch, qEI, fit_gpytorch_mll, SingleTaskGP, Sampler, MLL = _require_botorch()

    X_train = _as_float_matrix(tested_compounds, key="features")
    y_train = []
    for r in tested_compounds:
        if "activity" not in r:
            raise ValueError("Each tested_compounds entry must include 'activity'")
        y_train.append(float(r["activity"]))

    X_cand = _as_float_matrix(candidate_pool, key="features")
    ids = [str(r.get("id", "")) for r in candidate_pool]

    d = len(X_train[0])
    if any(len(x) != d for x in X_train) or any(len(x) != d for x in X_cand):
        raise ValueError("All feature vectors must have the same dimensionality")

    device = torch.device("cpu")
    dtype = torch.double
    X = torch.tensor(X_train, device=device, dtype=dtype)
    Y = torch.tensor(y_train, device=device, dtype=dtype).unsqueeze(-1)

    # Standardize features and outputs for stability.
    X_mean = X.mean(dim=0, keepdim=True)
    X_std = X.std(dim=0, keepdim=True).clamp_min(1e-12)
    Xs = (X - X_mean) / X_std

    Y_mean = Y.mean()
    Y_std = Y.std().clamp_min(1e-12)
    Ys = (Y - Y_mean) / Y_std

    model = SingleTaskGP(Xs, Ys)
    mll = MLL(model.likelihood, model)
    fit_gpytorch_mll(mll)

    Xc = torch.tensor(X_cand, device=device, dtype=dtype)
    Xcs = (Xc - X_mean) / X_std

    best_f = Ys.max()
    sampler = Sampler(sample_shape=torch.Size([256]))
    acq = qEI(model=model, best_f=best_f, sampler=sampler)

    # qEI expects q-batch; evaluate candidates individually.
    with torch.no_grad():
        acq_vals = acq(Xcs.unsqueeze(1)).squeeze(-1)  # (n,)
        post = model.posterior(Xcs)
        mu = post.mean.squeeze(-1)  # standardized
        var = post.variance.squeeze(-1).clamp_min(0.0)
        sigma = torch.sqrt(var)

    # Unstandardize predictive mean/std back to activity scale.
    mu_y = (mu * Y_std + Y_mean).cpu().numpy().tolist()
    sigma_y = (sigma * Y_std).cpu().numpy().tolist()
    acq_list = acq_vals.cpu().numpy().tolist()

    ranked = []
    for cid, a, m, s in zip(ids, acq_list, mu_y, sigma_y):
        if not cid:
            continue
        ranked.append({"id": cid, "acquisition": float(a), "pred_mean": float(m), "pred_std": float(s)})

    ranked.sort(key=lambda r: r["acquisition"], reverse=True)
    return ranked[:batch_size]


__all__ = ["recommend_next_compounds"]


