"""Train baseline ADMET models using TDC datasets.

Note: This script requires optional heavy dependencies (tdc, xgboost, rdkit).
Run it in an environment where those are installed.
"""

from __future__ import annotations

from typing import List


def _require_training_deps():
    try:
        import numpy as np  # type: ignore
        from sklearn.metrics import mean_squared_error, roc_auc_score  # type: ignore
        from tdc.single_pred import ADME, Tox  # type: ignore
        from xgboost import XGBClassifier, XGBRegressor  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError(
            "Training requires: tdc, xgboost, scikit-learn, numpy (and rdkit for feature extraction). "
            f"Underlying error: {type(e).__name__}: {e}"
        ) from e
    return np, roc_auc_score, mean_squared_error, ADME, Tox, XGBClassifier, XGBRegressor


def _featurize_smiles(predictor, smiles_list: List[str]):
    import numpy as np

    feats = [predictor._get_features(s) for s in smiles_list]  # noqa: SLF001
    mask = [x is not None for x in feats]
    X = np.array([x for x, ok in zip(feats, mask) if ok])
    return X, mask


def train_all() -> None:
    from amprenta_rag.ml.registry import get_registry
    from amprenta_rag.ml.admet.predictor import ADMETPredictor

    np, roc_auc_score, mean_squared_error, ADME, Tox, XGBClassifier, XGBRegressor = _require_training_deps()

    registry = get_registry()
    predictor = ADMETPredictor()

    # 1. hERG (classification)
    print("Training hERG model...")
    herg_data = Tox(name="hERG")
    train, _, test = herg_data.get_split()

    X_train, mask_train = _featurize_smiles(predictor, train["Drug"].tolist())
    y_train = train["Y"].values[mask_train]

    X_test, mask_test = _featurize_smiles(predictor, test["Drug"].tolist())
    y_test = test["Y"].values[mask_test]

    model_herg = XGBClassifier(n_estimators=200, max_depth=6, learning_rate=0.05)
    model_herg.fit(X_train, y_train)

    auc = roc_auc_score(y_test, model_herg.predict_proba(X_test)[:, 1])
    print(f"hERG AUC: {auc:.3f}")

    registry.register_model(
        name="admet_herg_xgb",
        version="v1.0.0",
        model_type="admet_classification",
        framework="xgboost",
        model_object=model_herg,
        metrics={"auc": float(auc)},
        description="hERG inhibition classifier (TDC dataset)",
    )

    # 2. LogS (regression)
    print("Training LogS model...")
    logs_data = ADME(name="Solubility_AqSolDB")
    train, _, test = logs_data.get_split()

    X_train, mask_train = _featurize_smiles(predictor, train["Drug"].tolist())
    y_train = train["Y"].values[mask_train].astype(float)

    X_test, mask_test = _featurize_smiles(predictor, test["Drug"].tolist())
    y_test = test["Y"].values[mask_test].astype(float)

    model_logs = XGBRegressor(n_estimators=400, max_depth=6, learning_rate=0.05)
    model_logs.fit(X_train, y_train)

    mse = mean_squared_error(y_test, model_logs.predict(X_test))
    rmse = float(np.sqrt(mse))
    print(f"LogS RMSE: {rmse:.3f}")

    registry.register_model(
        name="admet_logs_xgb",
        version="v1.0.0",
        model_type="admet_regression",
        framework="xgboost",
        model_object=model_logs,
        metrics={"rmse": float(rmse)},
        description="Aqueous solubility regressor (TDC dataset)",
    )

    # 3. LogP (regression)
    print("Training LogP model...")
    logp_data = ADME(name="Lipophilicity_AstraZeneca")
    train, _, test = logp_data.get_split()

    X_train, mask_train = _featurize_smiles(predictor, train["Drug"].tolist())
    y_train = train["Y"].values[mask_train].astype(float)

    X_test, mask_test = _featurize_smiles(predictor, test["Drug"].tolist())
    y_test = test["Y"].values[mask_test].astype(float)

    model_logp = XGBRegressor(n_estimators=400, max_depth=6, learning_rate=0.05)
    model_logp.fit(X_train, y_train)

    mse = mean_squared_error(y_test, model_logp.predict(X_test))
    rmse = float(np.sqrt(mse))
    print(f"LogP RMSE: {rmse:.3f}")

    registry.register_model(
        name="admet_logp_xgb",
        version="v1.0.0",
        model_type="admet_regression",
        framework="xgboost",
        model_object=model_logp,
        metrics={"rmse": float(rmse)},
        description="Lipophilicity regressor (TDC dataset)",
    )

    print("All models trained and registered!")


if __name__ == "__main__":
    train_all()


