"""Run all seed scripts in dependency order."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# Ensure repo root is on sys.path when running as `python scripts/seed_all.py`
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from amprenta_rag.database.models import (  # noqa: E402
    Dataset,
    DockingRun,
    Feature,
    Program,
    ProteinStructure,
    Signature,
    SingleCellDataset,
    VariantSet,
)
from amprenta_rag.database.session import db_session  # noqa: E402
from amprenta_rag.models.chemistry import Compound, HTSCampaign  # noqa: E402

from scripts.seed_utils import get_common_parser, safe_count, with_progress  # noqa: E402

SeederSpec = Tuple[str, str, str, Sequence[object]]


SEEDERS: List[SeederSpec] = [
    ("core", "scripts.seed_core_scaffolding", "main", (Program, Dataset, Compound)),
    ("transcriptomics", "scripts.seed_transcriptomics_data", "main", (Feature,)),
    ("proteomics", "scripts.seed_proteomics_data", "main", (Dataset,)),
    ("metabolomics", "scripts.seed_metabolomics_data", "main", (Dataset,)),
    ("lipidomics", "scripts.seed_lipidomics_data", "main", (Dataset,)),
    ("hts", "scripts.seed_hts_data", "main", (HTSCampaign,)),
    ("sar", "scripts.seed_sar_data", "main", (Compound,)),
    ("signatures", "scripts.seed_signatures", "main", (Signature,)),
    ("genomics", "scripts.seed_genomics_data", "main", (VariantSet,)),
    ("single_cell", "scripts.seed_single_cell_data", "main", (SingleCellDataset,)),
    ("structures", "scripts.seed_structures_data", "main", (ProteinStructure, DockingRun)),
]


def _build_args(size: str, reset: bool, seed: int, dry_run: bool) -> List[str]:
    args = ["--size", size, "--seed", str(seed)]
    if reset:
        args.append("--reset")
    if dry_run:
        args.append("--dry-run")
    return args


def _run_seeder(module_name: str, main_name: str, args: List[str]) -> None:
    mod = importlib.import_module(module_name)
    if not hasattr(mod, main_name):
        raise RuntimeError(f"{module_name} missing {main_name}")
    # simulate CLI by patching sys.argv
    old_argv = sys.argv
    sys.argv = [module_name] + args
    try:
        getattr(mod, main_name)()
    finally:
        sys.argv = old_argv


def _schema_ready(required_models: Sequence[object]) -> bool:
    # Check if all required tables exist; if any are missing, skip the seeder.
    with db_session() as db:
        for m in required_models:
            cnt = safe_count(db, m)  # None indicates missing table / inaccessible
            if cnt is None:
                return False
    return True


def _counts_snapshot(required_models: Sequence[object]) -> Dict[str, Optional[int]]:
    out: Dict[str, Optional[int]] = {}
    with db_session() as db:
        for m in required_models:
            tablename = getattr(m, "__tablename__", getattr(m, "__name__", m.__class__.__name__))
            out[str(tablename)] = safe_count(db, m)
    return out


def main() -> None:
    parser = get_common_parser("Run all seeders in dependency order.")
    parser.add_argument("--only", choices=[s[0] for s in SEEDERS], help="Run only the specified seeder.")
    args = parser.parse_args()

    results: Dict[str, str] = {}
    selected: List[SeederSpec] = []
    for spec in SEEDERS:
        if args.only and args.only != spec[0]:
            continue
        selected.append(spec)

    for name, module_name, main_name, required_models in with_progress(selected, desc="Seeding"):
        if not _schema_ready(required_models):
            results[name] = "skipped (schema missing)"
            print(f"=== Seeding {name} (skipped: schema missing) ===")
            continue
        if args.only and args.only != name:
            continue
        try:
            print(f"=== Seeding {name} ===")
            before = _counts_snapshot(required_models)
            _run_seeder(module_name, main_name, _build_args(args.size, args.reset, args.seed, args.dry_run))
            after = _counts_snapshot(required_models)
            results[name] = "ok"
            # Lightweight post-seed counts per domain
            if not args.dry_run:
                print("[counts]")
                for k in sorted(after.keys()):
                    b = before.get(k)
                    a = after.get(k)
                    if a is None:
                        continue
                    delta = (a - b) if (a is not None and b is not None) else None
                    if delta is None:
                        print(f"  {k}: {a}")
                    else:
                        print(f"  {k}: {a} (Δ {delta:+d})")
        except Exception as exc:  # pragma: no cover - runtime guard
            results[name] = f"error: {exc}"
            print(f"[ERROR] Seeder {name} failed: {exc}")

    print("=== Seed summary ===")
    for name in (args.only and [args.only] or [s[0] for s in SEEDERS]):
        status = results.get(name, "skipped")
        print(f"{name}: {status}")


if __name__ == "__main__":
    main()

