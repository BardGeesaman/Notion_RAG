from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any, Dict, List


def _ensure_repo_on_path() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))


def _parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Train per-target QSAR ensemble models.")
    p.add_argument(
        "--targets",
        type=str,
        default="common",
        help='Comma-separated target names or ChEMBL IDs, or "common" (default: common).',
    )
    p.add_argument("--source", choices=["chembl", "local"], default="chembl")
    p.add_argument("--threshold", type=float, default=1000.0, help="IC50 threshold (nM) for active label.")
    p.add_argument("--n-models", type=int, default=5, help="Bootstrap ensemble size.")
    p.add_argument("--output-dir", type=str, default="models/qsar", help="Artifact output directory.")
    p.add_argument("--dry-run", action="store_true", help="Print what would be trained without training.")
    return p.parse_args(argv)


def _parse_targets(arg: str, common_targets: Dict[str, str]) -> List[str]:
    s = (arg or "").strip()
    if not s or s.lower() == "common":
        return list(common_targets.keys())
    parts = [p.strip() for p in s.split(",")]
    return [p for p in parts if p]


def _print_summary(rows: List[Dict[str, Any]]) -> None:
    if not rows:
        print("No targets processed.")
        return

    cols = ["target", "status", "auc", "accuracy", "ece", "train_size", "source", "artifact"]
    widths = {c: max(len(c), max(len(str(r.get(c, ""))) for r in rows)) for c in cols}

    def fmt_row(r: Dict[str, Any]) -> str:
        return " | ".join(str(r.get(c, "")).ljust(widths[c]) for c in cols)

    header = " | ".join(c.ljust(widths[c]) for c in cols)
    sep = "-+-".join("-" * widths[c] for c in cols)
    print("\nSummary")
    print(header)
    print(sep)
    for r in rows:
        print(fmt_row(r))


def main(argv: List[str] | None = None) -> int:
    _ensure_repo_on_path()
    args = _parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Ensure ML registry artifacts (joblib) are stored under --output-dir for this run.
    os.environ["ML_ARTIFACT_PATH"] = str(out_dir)

    from amprenta_rag.ml.qsar.datasets import COMMON_TARGETS

    targets = _parse_targets(args.targets, COMMON_TARGETS)

    if args.dry_run:
        print("DRY RUN: would train QSAR models with:")
        print(f"- targets: {targets}")
        print(f"- source: {args.source}")
        print(f"- threshold_nm: {float(args.threshold)}")
        print(f"- n_models: {int(args.n_models)}")
        print(f"- output_dir: {out_dir}")
        return 0

    # Import after env var set so registry uses correct artifact path.
    from amprenta_rag.ml.qsar.trainer import train_target_model

    rows: List[Dict[str, Any]] = []

    for t in targets:
        try:
            out = train_target_model(
                target=str(t),
                source=str(args.source),
                threshold_nm=float(args.threshold),
                n_models=int(args.n_models),
                calibration_method="isotonic",
                register=True,
            )
            metrics = out.get("metrics") if isinstance(out, dict) else {}
            meta = out.get("metadata") if isinstance(out, dict) else {}
            rows.append(
                {
                    "target": str(t),
                    "status": "ok",
                    "auc": round(float(metrics.get("auc", float("nan"))), 3) if isinstance(metrics, dict) else "",
                    "accuracy": round(float(metrics.get("accuracy", float("nan"))), 3) if isinstance(metrics, dict) else "",
                    "ece": round(float(metrics.get("ece", float("nan"))), 3) if isinstance(metrics, dict) else "",
                    "train_size": int(meta.get("train_size", 0)) if isinstance(meta, dict) else "",
                    "source": str(args.source),
                    "artifact": f"{out_dir}/{out.get('model_name')}_{out.get('version')}.joblib",
                }
            )
        except Exception as e:  # noqa: BLE001
            rows.append(
                {
                    "target": str(t),
                    "status": f"error:{type(e).__name__}",
                    "auc": "",
                    "accuracy": "",
                    "ece": "",
                    "train_size": "",
                    "source": str(args.source),
                    "artifact": "",
                }
            )
            print(f"[WARN] Training failed for {t}: {e}", file=sys.stderr)
            continue

    _print_summary(rows)
    ok = sum(1 for r in rows if str(r.get("status")) == "ok")
    print(f"\nTrained successfully: {ok}/{len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


