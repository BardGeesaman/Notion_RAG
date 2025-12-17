"""Run all seed scripts in dependency order."""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from typing import Dict, List

# Ensure repo root is on sys.path when running as `python scripts/seed_all.py`
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SEEDERS = [
    ("core", "scripts.seed_core_scaffolding", "main"),
    ("transcriptomics", "scripts.seed_transcriptomics_data", "main"),
    ("proteomics", "scripts.seed_proteomics_data", "main"),
    ("metabolomics", "scripts.seed_metabolomics_data", "main"),
    ("lipidomics", "scripts.seed_lipidomics_data", "main"),
    ("hts", "scripts.seed_hts_data", "main"),
    ("sar", "scripts.seed_sar_data", "main"),
    ("signatures", "scripts.seed_signatures", "main"),
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Run all seeders in dependency order.")
    parser.add_argument("--size", choices=["small", "medium", "large"], default="small")
    parser.add_argument("--reset", action="store_true")
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--only", choices=[s[0] for s in SEEDERS], help="Run only the specified seeder.")
    args = parser.parse_args()

    results: Dict[str, str] = {}
    for name, module_name, main_name in SEEDERS:
        if args.only and args.only != name:
            continue
        try:
            print(f"=== Seeding {name} ===")
            _run_seeder(module_name, main_name, _build_args(args.size, args.reset, args.seed, args.dry_run))
            results[name] = "ok"
        except Exception as exc:  # pragma: no cover - runtime guard
            results[name] = f"error: {exc}"
            print(f"[ERROR] Seeder {name} failed: {exc}")

    print("=== Seed summary ===")
    for name in (args.only and [args.only] or [s[0] for s in SEEDERS]):
        status = results.get(name, "skipped")
        print(f"{name}: {status}")


if __name__ == "__main__":
    main()

