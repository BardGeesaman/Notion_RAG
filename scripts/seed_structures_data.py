"""Seed synthetic structural biology data (structures, pockets, docking runs/poses)."""

from __future__ import annotations

import argparse
import random
from datetime import datetime, timezone
from typing import Dict, List, Tuple

from sqlalchemy.orm import Session

from amprenta_rag.database.models import BindingSite, DockingPose, DockingRun, ProteinStructure, StructureFile
from amprenta_rag.database.session import db_session
from amprenta_rag.models.chemistry import Compound


SIZE_PRESETS: Dict[str, Tuple[int, int, int, int]] = {
    # structures, sites, runs, poses
    "small": (5, 10, 2, 20),
    "medium": (20, 40, 10, 100),
    "large": (100, 200, 50, 500),
}

METHODS = ["X-ray diffraction", "Cryo-EM", "NMR"]
CHAINS = [["A"], ["A", "B"], ["A", "B", "C"]]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed structural biology demo data.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo structures first.")
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _distribute(total: int, parts: int) -> List[int]:
    if parts <= 0:
        return []
    base = total // parts
    rem = total % parts
    return [base + (1 if i < rem else 0) for i in range(parts)]


def _reset_demo() -> Tuple[int, int, int, int]:
    """
    Delete demo structural objects.

    Order matters (FK constraints):
      DockingPose -> DockingRun -> BindingSite/StructureFile -> ProteinStructure
    """
    with db_session() as db:
        demo_structs = db.query(ProteinStructure).filter(ProteinStructure.source == "demo").all()
        struct_ids = [s.id for s in demo_structs]

        demo_runs = []
        if struct_ids:
            demo_runs = db.query(DockingRun).filter(DockingRun.structure_id.in_(struct_ids)).all()
        run_ids = [r.id for r in demo_runs]

        poses_deleted = 0
        if run_ids:
            poses_deleted = (
                db.query(DockingPose)
                .filter(DockingPose.docking_run_id.in_(run_ids))
                .delete(synchronize_session=False)
            )

        runs_deleted = 0
        if struct_ids:
            runs_deleted = (
                db.query(DockingRun)
                .filter(DockingRun.structure_id.in_(struct_ids))
                .delete(synchronize_session=False)
            )

        sites_deleted = 0
        if struct_ids:
            sites_deleted = (
                db.query(BindingSite)
                .filter(BindingSite.structure_id.in_(struct_ids))
                .delete(synchronize_session=False)
            )

        if struct_ids:
            (
                db.query(StructureFile)
                .filter(StructureFile.structure_id.in_(struct_ids))
                .delete(synchronize_session=False)
            )

        structs_deleted = 0
        if struct_ids:
            structs_deleted = (
                db.query(ProteinStructure)
                .filter(ProteinStructure.id.in_(struct_ids))
                .delete(synchronize_session=False)
            )

        # Delete demo docking compounds last (after poses)
        db.query(Compound).filter(Compound.compound_id.like("DEMO_DOCK_CMP_%")).delete(
            synchronize_session=False
        )

        db.commit()

    return structs_deleted, sites_deleted, runs_deleted, poses_deleted


def _ensure_compounds(db: Session, n: int) -> List[Compound]:
    existing = (
        db.query(Compound)
        .filter(Compound.compound_id.like("DEMO_DOCK_CMP_%"))
        .order_by(Compound.compound_id.asc())
        .all()
    )
    out: List[Compound] = list(existing)
    idx = len(existing) + 1
    while len(out) < n:
        cid = f"DEMO_DOCK_CMP_{idx:05d}"
        c = Compound(compound_id=cid, smiles="C")  # valid minimal placeholder
        db.add(c)
        out.append(c)
        idx += 1
    db.flush()
    return out[:n]


def _seed(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int, int, int]:
    rng = random.Random(seed_value)
    n_structs, total_sites, total_runs, total_poses = SIZE_PRESETS[size]

    sites_per_struct = _distribute(total_sites, n_structs)
    runs_per_struct = _distribute(total_runs, n_structs)
    poses_per_run = _distribute(total_poses, total_runs) if total_runs > 0 else []

    created_structs = created_sites = created_runs = created_poses = 0

    with db_session() as db:
        now = datetime.now(timezone.utc)

        # Create enough compounds to support docking poses.
        compounds = _ensure_compounds(db, max(1, total_poses))

        # Create structures + binding sites
        structures: List[ProteinStructure] = []
        binding_sites_by_struct: Dict[str, List[BindingSite]] = {}
        for i in range(n_structs):
            pdb_id = f"DEMO{i+1:04d}"
            existing = (
                db.query(ProteinStructure)
                .filter(ProteinStructure.source == "demo", ProteinStructure.pdb_id == pdb_id)
                .first()
            )
            if existing:
                # Idempotency: do not duplicate structures or their children on re-run without --reset.
                continue

            s = ProteinStructure(
                feature_id=None,
                pdb_id=pdb_id,
                alphafold_uniprot_id=None,
                source="demo",
                resolution=round(rng.uniform(1.2, 3.5), 2),
                method=rng.choice(METHODS),
                chain_ids=rng.choice(CHAINS),
                prep_status="prepared",
                prep_log="Seeded demo structure (no real file).",
                created_at=now,
                updated_at=now,
            )
            db.add(s)
            db.flush()
            created_structs += 1
            structures.append(s)

            # Add a placeholder structure file record for UX completeness.
            sf = StructureFile(
                structure_id=s.id,
                file_type="pdb",
                file_path=f"data/structures/{s.id}/raw.pdb",
                file_size_bytes=rng.randint(50_000, 3_000_000),
                md5_hash=None,
                created_at=now,
            )
            db.add(sf)

            bs_list: List[BindingSite] = []
            for rank in range(1, sites_per_struct[i] + 1):
                cx = rng.uniform(-10, 10)
                cy = rng.uniform(-10, 10)
                cz = rng.uniform(-10, 10)
                bs = BindingSite(
                    structure_id=s.id,
                    pocket_rank=rank,
                    score=round(rng.uniform(5, 50), 3),
                    volume=round(rng.uniform(200, 2000), 2),
                    center_x=round(cx, 3),
                    center_y=round(cy, 3),
                    center_z=round(cz, 3),
                    residues=None,
                    pocket_pdb_path=None,
                    detection_method="fpocket",
                    detected_at=now,
                )
                db.add(bs)
                bs_list.append(bs)
                created_sites += 1
            db.flush()
            binding_sites_by_struct[str(s.id)] = bs_list

        # Create runs + poses
        pose_idx = 0
        run_global_idx = 0
        for s_idx, s in enumerate(structures):
            for _ in range(runs_per_struct[s_idx]):
                nposes = poses_per_run[run_global_idx] if run_global_idx < len(poses_per_run) else 0
                run_global_idx += 1

                sites = binding_sites_by_struct.get(str(s.id)) or []
                bs = rng.choice(sites) if sites else None
                cx = bs.center_x if bs is not None else rng.uniform(-10, 10)
                cy = bs.center_y if bs is not None else rng.uniform(-10, 10)
                cz = bs.center_z if bs is not None else rng.uniform(-10, 10)

                status = rng.choices(
                    ["completed", "running", "failed"],
                    weights=[0.75, 0.15, 0.10],
                    k=1,
                )[0]
                total_compounds = nposes
                completed = nposes if status == "completed" else max(0, int(nposes * rng.uniform(0.0, 0.9)))

                selected_compounds = compounds[pose_idx : pose_idx + nposes]
                compound_ids_json = [str(c.id) for c in selected_compounds]

                dr = DockingRun(
                    structure_id=s.id,
                    binding_site_id=bs.id if bs is not None else None,
                    center_x=float(cx) if cx is not None else None,
                    center_y=float(cy) if cy is not None else None,
                    center_z=float(cz) if cz is not None else None,
                    size_x=20.0,
                    size_y=20.0,
                    size_z=20.0,
                    status=status,
                    compound_ids=compound_ids_json,
                    exhaustiveness=8,
                    total_compounds=total_compounds,
                    completed_compounds=completed,
                    started_at=now if status in ("running", "completed", "failed") else None,
                    completed_at=now if status == "completed" else None,
                    error_log="Seeded failure." if status == "failed" else None,
                )
                db.add(dr)
                db.flush()
                created_runs += 1

                for rank in range(1, nposes + 1):
                    c = compounds[pose_idx]
                    pose_idx += 1
                    dp = DockingPose(
                        docking_run_id=dr.id,
                        compound_id=c.id,
                        binding_affinity=round(rng.uniform(-12.0, -4.0), 3),
                        pose_rank=rank,
                        rmsd_lb=round(rng.uniform(0.0, 2.0), 3),
                        rmsd_ub=round(rng.uniform(0.0, 3.0), 3),
                        pose_pdbqt_path=f"data/docking/{dr.id}/{c.id}.pdbqt",
                        docked_at=now,
                    )
                    db.add(dp)
                    created_poses += 1

        if not dry_run:
            db.commit()
        else:
            db.rollback()

    return created_structs, created_sites, created_runs, created_poses


def main() -> None:
    args = _parse_args()
    if args.reset:
        structs_deleted, sites_deleted, runs_deleted, poses_deleted = _reset_demo()
        print(f"Reset: deleted ProteinStructure rows = {structs_deleted}")
        print(f"Reset: deleted BindingSite rows = {sites_deleted}")
        print(f"Reset: deleted DockingRun rows = {runs_deleted}")
        print(f"Reset: deleted DockingPose rows = {poses_deleted}")

    created_structs, created_sites, created_runs, created_poses = _seed(
        args.size, args.seed, args.dry_run
    )
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"ProteinStructure created: {created_structs}")
    print(f"BindingSite created: {created_sites}")
    print(f"DockingRun created: {created_runs}")
    print(f"DockingPose created: {created_poses}")


if __name__ == "__main__":
    main()


