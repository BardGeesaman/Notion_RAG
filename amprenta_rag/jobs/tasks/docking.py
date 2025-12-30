"""Celery tasks for molecular docking."""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from uuid import UUID

from amprenta_rag.jobs.celery_app import celery_app


@celery_app.task(bind=True, max_retries=2, default_retry_delay=120, queue='high')
def run_docking(self, run_id: str) -> dict:
    """Execute AutoDock Vina docking run."""
    from sqlalchemy import update
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import DockingRun, Compound
    from amprenta_rag.structural.receptor_prep import prepare_receptor
    
    run_uuid = UUID(run_id)
    
    try:
        # Mark running
        with db_session() as db:
            run = db.query(DockingRun).filter_by(id=run_uuid).first()
            if not run:
                return {"status": "failed", "error": "Run not found", "run_id": run_id}
            run.status = "running"
            run.started_at = datetime.now(timezone.utc)
            run.error_log = None
            db.add(run)
            db.commit()

        base = Path(os.environ.get("DOCKING_STORE_BASE", "data/docking")) / str(run_uuid)
        base.mkdir(parents=True, exist_ok=True)

        # Load run + structure for receptor prep
        with db_session() as db:
            run = db.query(DockingRun).filter_by(id=run_uuid).first()
            if not run or not run.structure:
                return {"status": "failed", "error": "Run or structure not found", "run_id": run_id}
            
            receptor_pdbqt = base / "receptor.pdbqt"
            if not receptor_pdbqt.exists():
                ok = prepare_receptor(run.structure, str(receptor_pdbqt))
                if not ok:
                    run.status = "failed"
                    run.error_log = "Failed to prepare receptor (need structure file + obabel)."
                    run.completed_at = datetime.now(timezone.utc)
                    db.add(run)
                    db.commit()
                    return {"status": "failed", "error": "Receptor preparation failed", "run_id": run_id}

            # Select compounds to dock
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
            db.commit()

        center = (run.center_x or 0.0, run.center_y or 0.0, run.center_z or 0.0)
        size = (run.size_x or 20.0, run.size_y or 20.0, run.size_z or 20.0)
        exhaustiveness = int(getattr(run, "exhaustiveness", 8) or 8)

        # Import the _process_compound method from DockingService
        from amprenta_rag.structural.docking_service import DockingService
        service = DockingService()

        # Process compounds in parallel
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {
                executor.submit(
                    service._process_compound,
                    run_uuid,
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
                            db.commit()
                except Exception as e:
                    with db_session() as db:
                        run2 = db.query(DockingRun).filter_by(id=run_uuid).first()
                        if run2:
                            run2.error_log = (run2.error_log or "") + f"\nCompound {cid}: {e}"
                            db.add(run2)
                            db.commit()
                finally:
                    with db_session() as db:
                        db.execute(
                            update(DockingRun)
                            .where(DockingRun.id == run_uuid)
                            .values(completed_compounds=DockingRun.completed_compounds + 1)
                        )
                        db.commit()

        # Mark complete
        with db_session() as db:
            run = db.query(DockingRun).filter_by(id=run_uuid).first()
            if not run:
                return {"status": "failed", "error": "Run not found", "run_id": run_id}
            run.status = "completed" if not run.error_log else "failed"
            run.completed_at = datetime.now(timezone.utc)
            db.add(run)
            db.commit()

        return {"status": "completed", "run_id": run_id, "compounds_processed": len(compounds)}

    except Exception as exc:
        # Update run status to failed
        try:
            with db_session() as db:
                run = db.query(DockingRun).filter_by(id=run_uuid).first()
                if run is not None:
                    run.status = "failed"
                    run.completed_at = datetime.now(timezone.utc)
                    run.error_log = (run.error_log or "") + f"\nTask error: {exc}"
                    db.add(run)
                    db.commit()
        except Exception:
            pass  # Don't let DB errors prevent retry logic
        
        # Retry with exponential backoff
        if self.request.retries >= self.max_retries:
            return {"status": "failed", "error": str(exc), "run_id": run_id}
        
        self.retry(exc=exc, countdown=120 * (2 ** self.request.retries))
