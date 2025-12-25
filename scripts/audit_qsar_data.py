"""Audit biochemical assay data availability for per-target QSAR training.

This script groups BiochemicalResult by target and reports:
- distinct compounds per target
- active compounds per target (IC50 < threshold, normalized to nM)
- active ratio and whether target meets training criteria
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

def _ensure_repo_on_path() -> None:
    # When executed as a script ("python scripts/xyz.py"), Python's sys.path[0]
    # is the scripts/ directory, not the repo root. Ensure repo root is importable.
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))


_ensure_repo_on_path()

from amprenta_rag.database.session import db_session  # noqa: E402
from amprenta_rag.models.chemistry import BiochemicalResult  # noqa: E402


def normalize_to_nm(value: float, unit: str) -> float:
    """Normalize concentration units to nM.

    Supported units: nM, uM/μM, mM, M (case-insensitive).
    """
    if value is None:
        raise ValueError("value is None")
    if unit is None:
        raise ValueError("unit is None")

    u = str(unit).strip().lower()
    u = u.replace("µ", "μ")  # normalize micro sign

    factors = {
        "nm": 1.0,
        "um": 1_000.0,
        "μm": 1_000.0,
        "mm": 1_000_000.0,
        "m": 1_000_000_000.0,
    }
    if u not in factors:
        raise ValueError(f"Unsupported unit: {unit}")
    return float(value) * factors[u]


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Audit biochemical IC50 data for QSAR feasibility.")
    p.add_argument("--threshold", type=float, default=1000.0, help="IC50 threshold in nM (default: 1000).")
    p.add_argument("--min-compounds", type=int, default=100, help="Minimum distinct compounds per target (default: 100).")
    p.add_argument("--min-active-ratio", type=float, default=0.2, help="Minimum active ratio per target (default: 0.2).")
    p.add_argument("--output", type=str, default="data/qsar_audit.csv", help="Output CSV path.")
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with db_session() as db:
            rows = (
                db.query(
                    BiochemicalResult.target,
                    BiochemicalResult.compound_id,
                    BiochemicalResult.ic50,
                    BiochemicalResult.units,
                )
                .filter(BiochemicalResult.ic50.isnot(None))
                .all()
            )
    except Exception as e:  # noqa: BLE001
        print(f"ERROR: could not connect/query database: {e}", file=sys.stderr)
        return 1

    if not rows:
        print("WARNING: No BiochemicalResult IC50 data found; writing empty audit CSV.")
        df_out = pd.DataFrame(
            columns=["target", "compound_count", "active_count", "active_ratio", "meets_criteria"]
        )
        df_out.to_csv(out_path, index=False)
        print(f"Wrote: {out_path}")
        return 0

    df = pd.DataFrame(rows, columns=["target", "compound_id", "ic50", "units"])
    df["target"] = df["target"].fillna("UNKNOWN").astype(str)

    def _to_nm(row) -> float:
        try:
            return normalize_to_nm(float(row["ic50"]), str(row["units"]))
        except Exception:
            return float("nan")

    df["ic50_nm"] = df.apply(_to_nm, axis=1)
    df = df.dropna(subset=["ic50_nm"])

    if df.empty:
        print("WARNING: No IC50 rows could be normalized to nM; writing empty audit CSV.")
        df_out = pd.DataFrame(
            columns=["target", "compound_count", "active_count", "active_ratio", "meets_criteria"]
        )
        df_out.to_csv(out_path, index=False)
        print(f"Wrote: {out_path}")
        return 0

    # For each (target, compound), take the minimum IC50 (best potency) to define activity.
    per_cmp = (
        df.groupby(["target", "compound_id"], as_index=False)["ic50_nm"].min()
    )
    per_cmp["is_active"] = per_cmp["ic50_nm"] < float(args.threshold)

    grouped = per_cmp.groupby("target").agg(
        compound_count=("compound_id", "nunique"),
        active_count=("is_active", "sum"),
    )
    grouped["active_ratio"] = grouped["active_count"] / grouped["compound_count"]
    grouped["meets_criteria"] = (grouped["compound_count"] >= int(args.min_compounds)) & (
        grouped["active_ratio"] >= float(args.min_active_ratio)
    )

    out = grouped.reset_index().sort_values(
        by=["meets_criteria", "compound_count", "active_ratio"],
        ascending=[False, False, False],
    )

    out = out[["target", "compound_count", "active_count", "active_ratio", "meets_criteria"]]
    out.to_csv(out_path, index=False)

    # Print summary table to stdout (top 25).
    pd.set_option("display.max_columns", 10)
    pd.set_option("display.width", 140)
    print(out.head(25).to_string(index=False))
    meets = int(out["meets_criteria"].sum())
    print(f"\nTargets meeting criteria: {meets} / {len(out)}")
    print(f"Wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


