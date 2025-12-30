"""Docking service (Vina) orchestration."""

from __future__ import annotations

import logging
import os
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from uuid import UUID

from sqlalchemy import update

from amprenta_rag.database.models import Compound, DockingPose, DockingRun
from amprenta_rag.database.session import db_session
from amprenta_rag.structural.ligand_prep import prepare_ligand
from amprenta_rag.structural.receptor_prep import prepare_receptor
from amprenta_rag.structural.vina_parser import parse_vina_output
from amprenta_rag.structural.vina_runner import run_vina


logger = logging.getLogger(__name__)


class DockingService:
    """Background docking runner for DockingRun jobs."""

    def __init__(self, max_workers: int = 4) -> None:
        self.max_workers = max_workers

    def start_run(self, run_id: UUID) -> None:
        # Submit job via Celery or fallback to threading
        import os
        if os.environ.get("USE_CELERY", "true").lower() == "true":
            from amprenta_rag.jobs.tasks.docking import run_docking
            run_docking.delay(str(run_id))
        else:
            # Fallback to threading for gradual rollout
            thread = threading.Thread(target=self._run_async, args=(run_id,), daemon=True)
            thread.start()

    def _run_async(self, run_id: UUID) -> None:
        try:
            with db_session() as db:
                run = db.query(DockingRun).filter_by(id=run_id).first()
                if not run:
                    return
                run.status = "running"
                run.started_at = datetime.now(timezone.utc)
                run.error_log = None

            base = Path(os.environ.get("DOCKING_STORE_BASE", "data/docking")) / str(run_id)
            base.mkdir(parents=True, exist_ok=True)

            # Load run + structure for receptor prep
            with db_session() as db:
                run = db.query(DockingRun).filter_by(id=run_id).first()
                if not run or not run.structure:
                    return
                receptor_pdbqt = base / "receptor.pdbqt"
                if not receptor_pdbqt.exists():
                    ok = prepare_receptor(run.structure, str(receptor_pdbqt))
                    if not ok:
                        run.status = "failed"
                        run.error_log = "Failed to prepare receptor (need structure file + obabel)."
                        run.completed_at = datetime.now(timezone.utc)
                        db.add(run)
                        return

                # Select compounds to dock:
                # - prefer explicit run.compound_ids (list of UUID strings)
                # - otherwise cap at 100 for safety
                compounds = []
                if isinstance(getattr(run, "compound_ids", None), list) and run.compound_ids:
                    try:
                        ids = [UUID(str(x)) for x in run.compound_ids]
                        compounds = db.query(Compound).filter(Compound.id.in_(ids)).all()
                    except Exception:
                        compounds = []
                if not compounds:
                    limit = int(run.total_compounds) if run.total_compounds else 100
                    compounds = db.query(Compound).limit(limit).all()
                    # Persist the chosen set for reproducibility
                    run.compound_ids = [str(c.id) for c in compounds]
                run.total_compounds = len(compounds)
                run.completed_compounds = 0
                db.add(run)

            center = (run.center_x or 0.0, run.center_y or 0.0, run.center_z or 0.0)
            size = (run.size_x or 20.0, run.size_y or 20.0, run.size_z or 20.0)
            exhaustiveness = int(getattr(run, "exhaustiveness", 8) or 8)

            # Process in parallel
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {
                    executor.submit(
                        self._process_compound,
                        run_id,
                        c.id,
                        str(receptor_pdbqt),
                        center,
                        size,
                        exhaustiveness,
                        str(base),
                    ): c.id
                    for c in compounds
                }
                for fut in as_completed(futures):
                    cid = futures[fut]
                    try:
                        pose = fut.result()
                        if pose is not None:
                            with db_session() as db:
                                db.add(pose)
                    except Exception as e:  # noqa: BLE001
                        logger.exception("[DOCKING] compound %s failed", cid)
                        with db_session() as db:
                            run2 = db.query(DockingRun).filter_by(id=run_id).first()
                            if run2:
                                run2.error_log = (run2.error_log or "") + f"\nCompound {cid}: {e}"
                                db.add(run2)
                    finally:
                        with db_session() as db:
                            db.execute(
                                update(DockingRun)
                                .where(DockingRun.id == run_id)
                                .values(completed_compounds=DockingRun.completed_compounds + 1)
                            )

            with db_session() as db:
                run = db.query(DockingRun).filter_by(id=run_id).first()
                if not run:
                    return
                run.status = "completed" if not run.error_log else "failed"
                run.completed_at = datetime.now(timezone.utc)
                db.add(run)

        except Exception as e:  # noqa: BLE001
            logger.exception("[DOCKING] run %s failed", run_id)
            with db_session() as db:
                run = db.query(DockingRun).filter_by(id=run_id).first()
                if run:
                    run.status = "failed"
                    run.error_log = str(e)[:2000]
                    run.completed_at = datetime.now(timezone.utc)
                    db.add(run)

    def _process_compound(
        self,
        run_id: UUID,
        compound_id: UUID,
        receptor_pdbqt: str,
        center: tuple[float, float, float],
        size: tuple[float, float, float],
        exhaustiveness: int,
        run_dir: str,
    ) -> Optional[DockingPose]:
        # Load compound in this thread
        with db_session() as db:
            run = db.query(DockingRun).filter_by(id=run_id).first()
            compound = db.query(Compound).filter_by(id=compound_id).first()
            if not run or not compound:
                return None
            smiles = compound.smiles

        cdir = Path(run_dir) / str(compound_id)
        cdir.mkdir(parents=True, exist_ok=True)
        ligand_pdbqt = cdir / "ligand.pdbqt"
        if not ligand_pdbqt.exists():
            ok = prepare_ligand(smiles, str(ligand_pdbqt))
            if not ok:
                raise RuntimeError("Ligand preparation failed")

        ok = run_vina(
            receptor_pdbqt=str(receptor_pdbqt),
            ligand_pdbqt=str(ligand_pdbqt),
            center=center,
            size=size,
            exhaustiveness=int(exhaustiveness or 8),
            output_dir=str(cdir),
        )
        if not ok:
            raise RuntimeError("Vina execution failed")

        res = parse_vina_output(str(cdir / "out.pdbqt"), str(cdir / "vina.log"))
        pose = DockingPose(
            docking_run_id=run_id,
            compound_id=compound_id,
            binding_affinity=res.binding_affinity,
            rmsd_lb=res.rmsd_lb,
            rmsd_ub=res.rmsd_ub,
            pose_rank=1,
            pose_pdbqt_path=str(cdir / "out.pdbqt"),
            docked_at=datetime.now(timezone.utc),
        )
        return pose


__all__ = ["DockingService"]


