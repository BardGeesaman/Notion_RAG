"""Seed core scaffolding data (Programs, Experiments, Datasets, Users).

Idempotent, deterministic seeding for demo/testing environments.
"""

from __future__ import annotations

import argparse
import random
import re
from typing import List, Tuple

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Program, Experiment, Dataset, User
from amprenta_rag.auth.password import hash_password


SIZE_PRESETS = {
    "small": (1, 2, 4),   # programs, experiments, datasets
    "medium": (3, 6, 18),
    "large": (10, 30, 100),
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed core scaffolding demo data.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small", help="Seed size preset.")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo entities before seeding.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for deterministic generation.")
    parser.add_argument("--dry-run", action="store_true", help="Show actions without committing.")
    return parser.parse_args()


def _reset_demo_data() -> Tuple[int, int, int, int]:
    prog_deleted = exp_deleted = ds_deleted = user_deleted = 0
    with db_session() as db:
        ds_deleted = db.query(Dataset).filter(Dataset.name.like("DEMO_DS_%")).delete(synchronize_session=False)
        exp_deleted = db.query(Experiment).filter(Experiment.name.like("DEMO_EXP_%")).delete(synchronize_session=False)
        prog_deleted = db.query(Program).filter(Program.name.like("DEMO_PROG_%")).delete(synchronize_session=False)
        user_deleted = db.query(User).filter(User.email.like("%@demo.local")).delete(synchronize_session=False)
        db.commit()
    return prog_deleted, exp_deleted, ds_deleted, user_deleted


def _get_or_create_user(email: str, name: str, dry_run: bool) -> int | None:
    """Return user ID (or None if dry_run)."""
    with db_session() as db:
        user = db.query(User).filter(User.email == email).first()
        if user:
            return user.id
        # Current User model requires username + password_hash.
        base_username = re.sub(r"[^a-zA-Z0-9_]+", "_", email.split("@", 1)[0]).strip("_") or "demo"
        username = base_username
        # Ensure username uniqueness if the base already exists
        suffix = 1
        while db.query(User).filter(User.username == username).first() is not None:
            suffix += 1
            username = f"{base_username}_{suffix}"

        user = User(
            username=username,
            email=email,
            password_hash=hash_password("demo"),
            role="admin",
            is_active=True,
        )
        db.add(user)
        # Flush so user.id is available to reference as created_by_id even in dry-run.
        db.flush()
        user_id = user.id  # capture ID before session closes
        if dry_run:
            db.rollback()
            return None
        else:
            db.commit()
            return user_id


def _seed_entities(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int, int]:
    random.seed(seed_value)
    prog_target, exp_target, ds_target = SIZE_PRESETS[size]
    prog_created = exp_created = ds_created = 0

    creator_id = _get_or_create_user("demo-seed@demo.local", "Demo Seed User", dry_run)

    with db_session() as db:
        # Programs
        programs: List[Program] = []
        for idx in range(1, prog_target + 1):
            name = f"DEMO_PROG_{idx:03d}"
            existing = db.query(Program).filter(Program.name == name).first()
            if existing:
                programs.append(existing)
                continue
            p = Program(name=name, description=f"Demo program {idx}", created_by_id=creator_id)
            db.add(p)
            programs.append(p)
            prog_created += 1

        if not dry_run:
            db.commit()
            for p in programs:
                db.refresh(p)

        # Experiments
        experiments: List[Experiment] = []
        exps_per_prog = max(1, exp_target // len(programs))
        exp_counter = 1
        for prog in programs:
            for _ in range(exps_per_prog):
                if exp_counter > exp_target:
                    break
                name = f"DEMO_EXP_{exp_counter:03d}"
                existing = db.query(Experiment).filter(Experiment.name == name).first()
                if existing:
                    experiments.append(existing)
                else:
                    e = Experiment(name=name, description=f"Demo experiment {exp_counter}", created_by_id=creator_id)
                    e.programs.append(prog)
                    db.add(e)
                    experiments.append(e)
                    exp_created += 1
                exp_counter += 1

        if not dry_run:
            db.commit()
            for e in experiments:
                db.refresh(e)

        # Datasets
        datasets: List[Dataset] = []
        ds_per_exp = max(1, ds_target // max(1, len(experiments)))
        ds_counter = 1
        omics_types = ["transcriptomics", "proteomics", "metabolomics", "lipidomics"]

        for exp in experiments:
            for _ in range(ds_per_exp):
                if ds_counter > ds_target:
                    break
                name = f"DEMO_DS_{ds_counter:03d}"
                existing = db.query(Dataset).filter(Dataset.name == name).first()
                if existing:
                    datasets.append(existing)
                else:
                    ds = Dataset(
                        name=name,
                        description=f"Demo dataset {ds_counter}",
                        omics_type=random.choice(omics_types),
                        created_by_id=creator_id,
                    )
                    ds.experiments.append(exp)
                    # link to the first program of this experiment, if present
                    if exp.programs:
                        ds.programs.append(exp.programs[0])
                    db.add(ds)
                    datasets.append(ds)
                    ds_created += 1
                ds_counter += 1

        if not dry_run:
            db.commit()
        else:
            # Ensure db_session() doesn't auto-commit pending inserts in dry-run.
            db.rollback()

    return prog_created, exp_created, ds_created


def main() -> None:
    args = _parse_args()
    if args.reset:
        p, e, d, u = _reset_demo_data()
        print(f"Reset: programs={p}, experiments={e}, datasets={d}, users={u}")

    prog_created, exp_created, ds_created = _seed_entities(args.size, args.seed, args.dry_run)

    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"Programs created: {prog_created}")
    print(f"Experiments created: {exp_created}")
    print(f"Datasets created: {ds_created}")


if __name__ == "__main__":
    main()

