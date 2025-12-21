"""Automated discovery service for repository scanning."""
import logging
from datetime import datetime
from typing import Optional, List, Any, cast
from uuid import UUID
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import DiscoveryJob, DiscoveredStudy, Experiment
from amprenta_rag.ingestion.repositories.discovery import fetch_study_metadata
from amprenta_rag.ingestion.repositories import GEORepository, MetabolightsRepository

logger = logging.getLogger(__name__)

REPOSITORY_CLASSES = {
    "GEO": GEORepository,
    "MetaboLights": MetabolightsRepository,
}


def run_discovery_job(
    repository: str,
    query: str,
    max_results: int = 50,
    user_id: Optional[str] = None
) -> str:
    """
    Run a discovery job to scan a repository for studies.

    Returns: job_id
    """
    db_gen = get_db()
    db = next(db_gen)
    try:
        # Create job record
        job = DiscoveryJob(
            repository=repository,
            query=query,
            status="running",
            started_at=datetime.utcnow(),
            created_by_id=cast(Any, UUID(user_id)) if user_id and user_id != "test" else None,
        )
        db.add(job)
        db.commit()
        job_id = str(job.id)

        logger.info(f"[DISCOVERY] Started job {job_id}: {repository} query='{query}'")

        try:
            # Get repository class
            repo_class = REPOSITORY_CLASSES.get(repository)
            if not repo_class:
                raise ValueError(f"Unknown repository: {repository}")

            repo = repo_class()
            results = repo.search_studies(query, max_results=max_results)

            studies_found = 0
            for result in results:
                # Check if already discovered or imported
                existing = db.query(DiscoveredStudy).filter(
                    DiscoveredStudy.study_id == result.study_id,
                    DiscoveredStudy.repository == repository
                ).first()

                if existing:
                    continue

                # Check if already imported as experiment
                imported_exp = db.query(Experiment).filter(
                    Experiment.name.ilike(f"%{result.study_id}%")
                ).first()

                discovered = DiscoveredStudy(
                    job_id=job.id,
                    study_id=result.study_id,
                    repository=repository,
                    title=result.title[:500] if result.title else None,
                    description=result.description[:2000] if result.description else None,
                    organism=result.organism,
                    omics_type=result.omics_type,
                    status="imported" if imported_exp else "new",
                    imported_experiment_id=imported_exp.id if imported_exp else None
                )
                db.add(discovered)
                studies_found += 1

            job.studies_found = studies_found
            job.status = "completed"
            job.completed_at = datetime.utcnow()
            db.commit()

            logger.info(f"[DISCOVERY] Job {job_id} completed: {studies_found} studies found")

        except Exception as e:
            job.status = "failed"
            job.error_message = str(e)[:1000]
            job.completed_at = datetime.utcnow()
            db.commit()
            logger.error(f"[DISCOVERY] Job {job_id} failed: {e}")

        return job_id

    finally:
        db_gen.close()


def get_pending_studies(limit: int = 100) -> List[dict]:
    """Get studies pending review."""
    db_gen = get_db()
    db = next(db_gen)
    try:
        studies = db.query(DiscoveredStudy).filter(
            DiscoveredStudy.status == "new"
        ).order_by(DiscoveredStudy.discovered_at.desc()).limit(limit).all()

        return [
            {
                "id": str(s.id),
                "study_id": s.study_id,
                "repository": s.repository,
                "title": s.title,
                "organism": s.organism,
                "omics_type": s.omics_type,
                "discovered_at": s.discovered_at.isoformat() if s.discovered_at else None
            }
            for s in studies
        ]
    finally:
        db_gen.close()


def import_discovered_study(discovered_study_id: str) -> Optional[UUID]:
    """Import a discovered study as an experiment."""
    from amprenta_rag.ingestion.design_integration import create_experiment_from_study

    db_gen = get_db()
    db = next(db_gen)
    try:
        study = db.query(DiscoveredStudy).filter(DiscoveredStudy.id == discovered_study_id).first()
        if not study:
            return None

        # Fetch full metadata and create experiment
        metadata = fetch_study_metadata(str(study.study_id), str(study.repository))
        if not metadata:
            return None

        exp_id = create_experiment_from_study(metadata)
        if exp_id:
            study.status = "imported"
            setattr(study, "imported_experiment_id", exp_id)
            db.commit()
            logger.info(f"[DISCOVERY] Imported study {study.study_id} as experiment {exp_id}")

            # Fire workflow trigger
            from amprenta_rag.automation.engine import fire_trigger
            fire_trigger("discovery_imported", {
                "study_id": str(study.id),
                "title": study.title or ""
            }, db)

        return exp_id
    finally:
        db_gen.close()
