"""Bayesian optimization helpers (BoTorch).

This module provides a small, backend-agnostic recommendation function that
uses a GP surrogate + qExpectedImprovement to rank candidate compounds.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence


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


def _require_botorch_multi():
    """Import multi-objective BoTorch components."""
    try:
        import torch  # type: ignore
        from botorch.models import SingleTaskGP  # type: ignore
        from botorch.models.model_list_gp_regression import ModelListGP  # type: ignore
        from botorch.acquisition.multi_objective import qNoisyExpectedHypervolumeImprovement  # type: ignore
        from botorch.utils.multi_objective import is_non_dominated  # type: ignore
        from botorch.fit import fit_gpytorch_mll  # type: ignore
        from gpytorch.mlls import ExactMarginalLogLikelihood  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("BoTorch multi-objective requires: botorch>=0.9") from e
    return (torch, SingleTaskGP, ModelListGP, qNoisyExpectedHypervolumeImprovement, 
            is_non_dominated, fit_gpytorch_mll, ExactMarginalLogLikelihood)


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


def compute_pareto_front(objectives_matrix):
    """Compute Pareto front using BoTorch is_non_dominated.
    
    Args:
        objectives_matrix: List of lists or 2D array, shape (n_compounds, n_objectives)
        
    Returns:
        Dict with "indices" (list of int) and "values" (list of lists)
    """
    torch, _, _, _, is_non_dominated, _, _ = _require_botorch_multi()
    obj_tensor = torch.tensor(objectives_matrix, dtype=torch.double)
    mask = is_non_dominated(obj_tensor)
    indices = mask.nonzero(as_tuple=True)[0].tolist()
    values = obj_tensor[mask].tolist()
    return {"indices": indices, "values": values}


def recommend_multi_objective(
    tested_compounds: List[Dict],
    candidate_pool: List[Dict],
    objectives: List[str],
    objective_directions: List[str],
    batch_size: int = 10,
    ref_point: Optional[List[float]] = None,
) -> Dict:
    """Multi-objective Bayesian optimization using qNEHVI.
    
    Args:
        tested_compounds: List of dicts with "features" (list) and objective values
        candidate_pool: List of dicts with "id" and "features"
        objectives: List of objective names (e.g., ["potency", "herg"])
        objective_directions: List of "maximize" or "minimize" for each objective
        batch_size: Number of recommendations to return
        ref_point: Optional reference point for hypervolume (auto-computed if None)
        
    Returns:
        Dict with:
            "recommendations": List of dicts with id, acquisition, pred_means, pred_stds
            "pareto_front": Dict with indices and values
    """
    # Validation
    if len(tested_compounds) < 5:
        raise ValueError("Need at least 5 tested_compounds to fit GPs")
    if len(objectives) != len(objective_directions):
        raise ValueError("objectives and objective_directions must have same length")
    if len(objectives) < 2:
        raise ValueError("Multi-objective requires at least 2 objectives")
    for direction in objective_directions:
        if direction not in ["maximize", "minimize"]:
            raise ValueError(f"Invalid direction: {direction}. Must be 'maximize' or 'minimize'")
    
    # Check all objectives present
    for compound in tested_compounds:
        for obj in objectives:
            if obj not in compound:
                raise ValueError(f"Objective '{obj}' missing from tested compound")
    
    if not candidate_pool:
        return {"recommendations": [], "pareto_front": {"indices": [], "values": []}}
    
    torch, SingleTaskGP, ModelListGP, qNEHVI, is_non_dominated, fit_gpytorch_mll, MLL = _require_botorch_multi()
    
    # Extract features and objectives
    X_train = _as_float_matrix(tested_compounds, key="features")
    Y_train = []
    for compound in tested_compounds:
        y_vals = []
        for obj, direction in zip(objectives, objective_directions):
            val = float(compound[obj])
            # Flip sign for minimize objectives (turn into maximization)
            if direction == "minimize":
                val = -val
            y_vals.append(val)
        Y_train.append(y_vals)
    
    X_cand = _as_float_matrix(candidate_pool, key="features")
    ids = [str(r.get("id", "")) for r in candidate_pool]
    
    d = len(X_train[0])
    n_obj = len(objectives)
    
    if any(len(x) != d for x in X_train) or any(len(x) != d for x in X_cand):
        raise ValueError("All feature vectors must have the same dimensionality")
    
    device = torch.device("cpu")
    dtype = torch.double
    X = torch.tensor(X_train, device=device, dtype=dtype)
    Y = torch.tensor(Y_train, device=device, dtype=dtype)  # (n, n_obj)
    
    # Standardize features
    X_mean = X.mean(dim=0, keepdim=True)
    X_std = X.std(dim=0, keepdim=True).clamp_min(1e-12)
    Xs = (X - X_mean) / X_std
    
    # Build one GP per objective
    y_means_list = []
    y_stds_list = []
    models = []
    for i in range(n_obj):
        y_i = Y[:, i].unsqueeze(-1)
        y_mean = y_i.mean()
        y_std = y_i.std().clamp_min(1e-12)
        y_means_list.append(y_mean)
        y_stds_list.append(y_std)
        ys_i = (y_i - y_mean) / y_std
        
        gp = SingleTaskGP(Xs, ys_i)
        mll = MLL(gp.likelihood, gp)
        fit_gpytorch_mll(mll)
        models.append(gp)
    
    model = ModelListGP(*models)
    
    # Auto-compute reference point if not provided
    if ref_point is None:
        ref_point = []
        for i in range(n_obj):
            y_vals = Y[:, i].cpu().numpy()
            ref_val = float(y_vals.min() - 0.1 * y_vals.std())
            ref_point.append(ref_val)
    
    ref_point_tensor = torch.tensor(ref_point, device=device, dtype=dtype)
    
    # Prepare candidates
    Xc = torch.tensor(X_cand, device=device, dtype=dtype)
    Xcs = (Xc - X_mean) / X_std
    
    # qNEHVI acquisition
    acq = qNEHVI(
        model=model,
        ref_point=ref_point_tensor,
        X_baseline=Xs,
        prune_baseline=True,
    )
    
    # Evaluate acquisition for each candidate
    with torch.no_grad():
        acq_vals = []
        pred_means = []
        pred_stds = []
        
        for i in range(len(Xcs)):
            x_i = Xcs[i:i+1].unsqueeze(1)  # (1, 1, d)
            acq_val = acq(x_i).item()
            acq_vals.append(acq_val)
            
            # Get predictions for each objective
            means = []
            stds = []
            for j, gp in enumerate(models):
                post = gp.posterior(Xcs[i:i+1])
                mean_std = post.mean.item()  # standardized
                std_std = torch.sqrt(post.variance.clamp_min(0.0)).item()
                # Unstandardize to original scale
                mean = mean_std * y_stds_list[j].item() + y_means_list[j].item()
                std = std_std * y_stds_list[j].item()
                means.append(mean)
                stds.append(std)
            pred_means.append(means)
            pred_stds.append(stds)
    
    # Build recommendations
    ranked = []
    for cid, acq_val, means, stds in zip(ids, acq_vals, pred_means, pred_stds):
        if not cid:
            continue
        
        pred_means_dict = {}
        pred_stds_dict = {}
        for j, obj in enumerate(objectives):
            pred_means_dict[obj] = float(means[j])
            pred_stds_dict[obj] = float(stds[j])
        
        ranked.append({
            "id": cid,
            "acquisition": float(acq_val),
            "pred_means": pred_means_dict,
            "pred_stds": pred_stds_dict,
        })
    
    ranked.sort(key=lambda r: r["acquisition"], reverse=True)
    recommendations = ranked[:batch_size]
    
    # Compute Pareto front from tested compounds
    pareto_front = compute_pareto_front(Y_train)
    
    return {
        "recommendations": recommendations,
        "pareto_front": pareto_front,
    }


__all__ = ["recommend_next_compounds", "compute_pareto_front", "recommend_multi_objective"]


