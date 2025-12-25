from __future__ import annotations

import argparse
import json
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("pandas is required (pip install pandas)") from e
    return pd


def _require_sklearn():
    try:
        from sklearn.model_selection import train_test_split  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("scikit-learn is required (pip install scikit-learn)") from e
    return train_test_split


@dataclass(frozen=True)
class EndpointSpec:
    name: str
    sql: str
    kind: str  # "classification" | "regression"


BASE_SELECT = """
SELECT
  act.activity_id AS activity_id,
  act.molregno AS molregno,
  md.chembl_id AS molecule_chembl_id,
  ass.assay_id AS assay_id,
  ass.chembl_id AS assay_chembl_id,
  ass.assay_type AS assay_type,
  td.pref_name AS target_pref_name,
  td.chembl_id AS target_chembl_id,
  act.standard_type AS standard_type,
  act.standard_value AS standard_value,
  act.standard_units AS standard_units,
  act.standard_relation AS standard_relation,
  cs.canonical_smiles AS canonical_smiles,
  cp.full_mwt AS full_mwt,
  cp.alogp AS alogp,
  cp.psa AS psa
FROM activities act
JOIN assays ass ON ass.assay_id = act.assay_id
JOIN target_dictionary td ON td.tid = ass.tid
JOIN molecule_dictionary md ON md.molregno = act.molregno
LEFT JOIN compound_structures cs ON cs.molregno = act.molregno
LEFT JOIN compound_properties cp ON cp.molregno = act.molregno
WHERE
  act.standard_relation = '='
  AND act.standard_value IS NOT NULL
  AND act.standard_type IS NOT NULL
"""


ENDPOINTS: Dict[str, EndpointSpec] = {
    "herg": EndpointSpec(
        name="herg",
        kind="classification",
        sql=(
            BASE_SELECT
            + """
  AND td.pref_name LIKE '%KCNH2%'
  AND act.standard_type = 'IC50'
  AND ass.assay_type = 'B'
"""
        ),
    ),
    "logs": EndpointSpec(
        name="logs",
        kind="regression",
        sql=(
            BASE_SELECT
            + """
  AND act.standard_type = 'Solubility'
  AND act.standard_units IN ('uM', 'ug/mL')
"""
        ),
    ),
    "logp": EndpointSpec(
        name="logp",
        kind="regression",
        sql=(
            BASE_SELECT
            + """
  AND act.standard_type IN ('LogP', 'LogD')
"""
        ),
    ),
}


def _connect(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    return conn


def _run_query(conn: sqlite3.Connection, sql: str, *, limit: Optional[int] = None):
    pd = _require_pandas()
    q = sql.strip()
    if limit is not None:
        q = f"{q}\nLIMIT {int(limit)}"
    return pd.read_sql_query(q, conn)


def _filter_assays_min_size(df, *, min_compounds: int = 50):
    if df.empty:
        return df
    counts = df.groupby("assay_id")["molecule_chembl_id"].nunique().rename("n_compounds")
    keep = counts[counts >= int(min_compounds)].index
    out = df[df["assay_id"].isin(keep)].copy()
    out = out.merge(counts.reset_index(), on="assay_id", how="left")
    return out


def _dedup_median_per_compound(df, *, value_col: str = "standard_value"):
    pd = _require_pandas()
    if df.empty:
        return df
    df = df.copy()
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    df = df[df[value_col].notna()]

    # Keep one row per molecule_chembl_id using median for value + first non-null metadata columns.
    num = df.groupby("molecule_chembl_id")[value_col].median().rename("median_value")
    meta_cols = [
        "molregno",
        "canonical_smiles",
        "full_mwt",
        "alogp",
        "psa",
    ]
    meta = (
        df.sort_values(["molecule_chembl_id", "assay_id", "activity_id"])
        .groupby("molecule_chembl_id", as_index=True)[meta_cols]
        .first()
    )
    out = pd.concat([meta, num], axis=1).reset_index()
    return out


def _herg_to_uM(value: float, units: str | None, full_mwt: float | None) -> Optional[float]:  # noqa: ARG001
    if value is None:
        return None
    u = (units or "").strip()
    try:
        v = float(value)
    except Exception:  # noqa: BLE001
        return None

    if u == "nM":
        return v / 1000.0
    if u == "uM":
        return v
    # Unknown: skip
    return None


def _solubility_to_uM(value: float, units: str | None, full_mwt: float | None) -> Optional[float]:
    if value is None:
        return None
    u = (units or "").strip()
    try:
        v = float(value)
    except Exception:  # noqa: BLE001
        return None

    if u == "uM":
        return v
    if u == "ug/mL":
        # ug/mL == mg/L. Convert mg/L -> uM: uM = (mg/L) * 1000 / MW
        if full_mwt is None:
            return None
        try:
            mw = float(full_mwt)
        except Exception:  # noqa: BLE001
            return None
        if mw <= 0:
            return None
        return (v * 1000.0) / mw
    return None


def _downsample_majority_to_ratio(df, *, label_col: str, max_ratio: float = 2.0, seed: int = 42):
    pd = _require_pandas()
    if df.empty:
        return df
    vc = df[label_col].value_counts(dropna=False)
    if len(vc) < 2:
        return df
    maj_label = int(vc.idxmax())
    min_label = int(vc.idxmin())
    n_min = int(vc[min_label])
    n_maj_keep = int(min(float(max_ratio) * n_min, float(vc[maj_label])))
    df_min = df[df[label_col] == min_label]
    df_maj = df[df[label_col] == maj_label].sample(n=n_maj_keep, random_state=seed)
    out = pd.concat([df_min, df_maj], ignore_index=True).sample(frac=1.0, random_state=seed)
    return out


def _split(df, *, kind: str, label_col: str = "y", seed: int = 42):
    train_test_split = _require_sklearn()
    if df.empty:
        return df, df, df

    stratify = df[label_col] if kind == "classification" else None
    train, tmp = train_test_split(df, test_size=0.30, random_state=seed, stratify=stratify)
    stratify_tmp = tmp[label_col] if kind == "classification" else None
    cal, test = train_test_split(tmp, test_size=0.50, random_state=seed, stratify=stratify_tmp)
    return train, cal, test


def _ensure_parquet_support() -> None:
    try:
        import pyarrow  # noqa: F401
    except Exception as e:  # noqa: BLE001
        raise ImportError("Parquet output requires pyarrow (pip install pyarrow)") from e


def _write_parquet(df, path: Path) -> None:
    _ensure_parquet_support()
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)


def curate_endpoint(conn: sqlite3.Connection, spec: EndpointSpec, *, min_assay_compounds: int = 50) -> Tuple[Any, Dict[str, Any]]:
    pd = _require_pandas()

    raw = _run_query(conn, spec.sql)
    raw["standard_units"] = raw["standard_units"].astype(str)

    filtered = _filter_assays_min_size(raw, min_compounds=min_assay_compounds)

    # Endpoint-specific normalization
    if spec.name == "herg":
        filtered["value_uM"] = filtered.apply(
            lambda r: _herg_to_uM(r["standard_value"], r["standard_units"], r.get("full_mwt")),
            axis=1,
        )
        filtered = filtered[filtered["value_uM"].notna()].copy()
        # label: IC50 < 10uM => 1 else 0
        filtered["y"] = (filtered["value_uM"].astype(float) < 10.0).astype(int)
        # downsample majority to 2:1
        filtered = _downsample_majority_to_ratio(filtered, label_col="y", max_ratio=2.0, seed=42)

        # dedup per compound using median IC50_uM
        med = filtered.groupby("molecule_chembl_id")["value_uM"].median().rename("y_value_uM")
        meta = (
            filtered.sort_values(["molecule_chembl_id", "assay_id", "activity_id"])
            .groupby("molecule_chembl_id", as_index=True)[["molregno", "canonical_smiles", "full_mwt", "alogp", "psa"]]
            .first()
        )
        out = pd.concat([meta, med], axis=1).reset_index()
        out["y"] = (out["y_value_uM"].astype(float) < 10.0).astype(int)
    elif spec.name == "logs":
        filtered["value_uM"] = filtered.apply(
            lambda r: _solubility_to_uM(r["standard_value"], r["standard_units"], r.get("full_mwt")),
            axis=1,
        )
        filtered = filtered[filtered["value_uM"].notna()].copy()

        out = _dedup_median_per_compound(filtered.rename(columns={"value_uM": "standard_value_uM"}), value_col="standard_value_uM")
        out = out.rename(columns={"median_value": "y"})
        out["endpoint_units"] = "uM"
    else:
        # logp/logd: use standard_value as y
        filtered["standard_value"] = pd.to_numeric(filtered["standard_value"], errors="coerce")
        filtered = filtered[filtered["standard_value"].notna()].copy()
        out = _dedup_median_per_compound(filtered, value_col="standard_value")
        out = out.rename(columns={"median_value": "y"})

    summary: Dict[str, Any] = {
        "endpoint": spec.name,
        "kind": spec.kind,
        "raw_rows": int(len(raw)),
        "post_assay_filter_rows": int(len(filtered)),
        "dedup_rows": int(len(out)),
    }
    if spec.kind == "classification" and not out.empty:
        vc = out["y"].value_counts().to_dict()
        summary["label_counts"] = {str(k): int(v) for k, v in vc.items()}
    return out, summary


def main() -> None:
    ap = argparse.ArgumentParser(description="Curate ChEMBL ADMET datasets (hERG/LogS/LogP) from ChEMBL SQLite.")
    ap.add_argument("--db", default="data/chembl/chembl_33.db", help="Path to chembl_33.db")
    ap.add_argument("--out-dir", default="data/admet", help="Output directory for parquet + summary.json")
    ap.add_argument("--min-assay-compounds", type=int, default=50, help="Minimum distinct compounds per assay")
    ap.add_argument("--dry-run", action="store_true", help="Print SQL queries and run small sanity LIMIT checks if DB exists")
    ap.add_argument("--limit", type=int, default=0, help="Optional limit (per endpoint) for faster local iteration")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    db_path = Path(args.db)

    manifest: Dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "db_path": str(db_path),
        "dry_run": bool(args.dry_run),
        "endpoints": {},
    }

    if args.dry_run:
        print("DRY RUN: curate chembl admet")
        for name, spec in ENDPOINTS.items():
            print(f"\n--- {name} SQL ---\n{spec.sql.strip()}\n")

        if not db_path.exists():
            manifest["db_found"] = False
            (out_dir / "summary.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
            print(f"DB not found at {db_path}; wrote dry-run summary to {out_dir / 'summary.json'}")
            return

        manifest["db_found"] = True
        conn = _connect(db_path)
        try:
            for name, spec in ENDPOINTS.items():
                df = _run_query(conn, spec.sql, limit=5)
                manifest["endpoints"][name] = {"sanity_rows": int(len(df))}
        finally:
            conn.close()
        (out_dir / "summary.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
        print(f"Wrote dry-run summary to {out_dir / 'summary.json'}")
        return

    if not db_path.exists():
        raise FileNotFoundError(f"ChEMBL SQLite not found at: {db_path} (run scripts/download_chembl_sqlite.py)")

    conn = _connect(db_path)
    try:
        for name, spec in ENDPOINTS.items():
            df, summary = curate_endpoint(conn, spec, min_assay_compounds=int(args.min_assay_compounds))
            if args.limit and int(args.limit) > 0 and len(df) > int(args.limit):
                df = df.sample(n=int(args.limit), random_state=42).reset_index(drop=True)
                summary["limited_to"] = int(args.limit)

            # standard columns across endpoints
            df["endpoint"] = spec.name

            train, cal, test = _split(df, kind=spec.kind, label_col="y", seed=42)

            _write_parquet(train, out_dir / f"{spec.name}_train.parquet")
            _write_parquet(cal, out_dir / f"{spec.name}_cal.parquet")
            _write_parquet(test, out_dir / f"{spec.name}_test.parquet")

            manifest["endpoints"][name] = {
                **summary,
                "split_sizes": {"train": int(len(train)), "cal": int(len(cal)), "test": int(len(test))},
            }

            if spec.kind == "classification" and not df.empty:
                manifest["endpoints"][name]["split_label_counts"] = {
                    "train": {str(k): int(v) for k, v in train["y"].value_counts().to_dict().items()},
                    "cal": {str(k): int(v) for k, v in cal["y"].value_counts().to_dict().items()},
                    "test": {str(k): int(v) for k, v in test["y"].value_counts().to_dict().items()},
                }

        (out_dir / "summary.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
        print(f"Wrote: {out_dir / 'summary.json'}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()


