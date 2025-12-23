"""Nextflow pipeline execution backend."""

from __future__ import annotations

import subprocess
import threading
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.database.models import NextflowJob
from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def check_nextflow_installed() -> Dict[str, Any]:
    """Check if Nextflow and Java are installed.

    Returns:
        {"installed": bool, "version": str, "java_version": str, "error": str}
    """
    result: Dict[str, Any] = {"installed": False, "version": None, "java_version": None, "error": None}

    # Check Java
    try:
        java_proc = subprocess.run(
            ["java", "-version"],
            capture_output=True,
            text=True,
            check=False,
        )
        if java_proc.returncode == 0:
            # Java version is in stderr
            result["java_version"] = (java_proc.stderr or "").split("\n")[0]
    except FileNotFoundError:
        result["error"] = "Java not found. Install: conda install -c conda-forge openjdk=17"
        return result

    # Check Nextflow
    try:
        nf_proc = subprocess.run(
            ["nextflow", "-version"],
            capture_output=True,
            text=True,
            check=False,
        )
        if nf_proc.returncode == 0:
            result["installed"] = True
            result["version"] = (nf_proc.stdout or "").strip()
    except FileNotFoundError:
        result["error"] = "Nextflow not found. Install: curl -s https://get.nextflow.io | bash"
        return result

    return result


def list_nfcore_pipelines() -> List[Dict[str, str]]:
    """Return list of supported nf-core pipelines."""
    return [
        {"name": "nf-core/rnaseq", "version": "3.14.0", "description": "RNA-seq analysis"},
        # Future: add more pipelines
    ]


def generate_sample_sheet(fastq_files: List[Dict[str, str]], output_path: Path) -> Path:
    """Generate nf-core sample sheet CSV.

    Args:
        fastq_files: List of {"sample": str, "fastq_1": str, "fastq_2": str (optional)}
        output_path: Where to save CSV

    Returns:
        Path to generated sample_sheet.csv
    """
    import csv

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
        for entry in fastq_files:
            writer.writerow(
                [
                    entry["sample"],
                    entry["fastq_1"],
                    entry.get("fastq_2", ""),
                    entry.get("strandedness", "auto"),
                ]
            )

    return output_path


def submit_nextflow_job(
    pipeline: str,
    sample_sheet: Path,
    genome: str,
    output_dir: Path,
    user_id: str,
    params: Optional[Dict[str, Any]] = None,
    profile: str = "docker",
) -> UUID:
    """Submit Nextflow job for background execution.

    Returns:
        job_id (UUID)
    """
    params = params or {}
    params.setdefault("profile", profile)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create job record
    with db_session() as db:
        job = NextflowJob(
            pipeline_name=pipeline,
            pipeline_version=(params.get("version", "latest") if isinstance(params, dict) else "latest"),
            status="pending",
            sample_sheet_path=str(sample_sheet),
            genome=genome,
            output_dir=str(output_dir),
            work_dir=str(output_dir / "work"),
            nextflow_log=str(output_dir / ".nextflow.log"),
            created_by=user_id,
            created_at=datetime.utcnow(),
            params=params,
        )
        db.add(job)
        db.flush()
        db.refresh(job)
        job_id = job.id

    # Start background thread
    thread = threading.Thread(
        target=_run_nextflow_async,
        args=(job_id,),
        daemon=True,
    )
    thread.start()

    return job_id


def _run_nextflow_async(job_id: UUID) -> None:
    """Background thread to execute Nextflow pipeline."""
    try:
        with db_session() as db:
            job = db.query(NextflowJob).filter_by(id=job_id).first()
            if not job:
                return

            job.status = "running"
            job.started_at = datetime.utcnow()
            job.progress_percent = 5
            db.add(job)

            # Build command
            profile = (job.params or {}).get("profile", "docker") if isinstance(job.params, dict) else "docker"
            cmd = [
                "nextflow",
                "run",
                job.pipeline_name,
                "--input",
                job.sample_sheet_path,
                "--genome",
                job.genome,
                "--outdir",
                job.output_dir,
                "-profile",
                profile,
                "-work-dir",
                job.work_dir,
            ]

            if isinstance(job.params, dict) and job.params.get("resume"):
                cmd.append("-resume")

            logger.info("[NEXTFLOW] Running: %s", " ".join(cmd))

            # Execute (capture output; Nextflow itself writes .nextflow.log under cwd)
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
                cwd=job.output_dir,
            )

            # Update progress based on log if available
            try:
                job.progress_percent = parse_nextflow_progress(Path(job.nextflow_log))
            except Exception:
                pass

            if proc.returncode == 0:
                job.status = "complete"
                job.progress_percent = 100
                # Find MultiQC report
                multiqc = Path(job.output_dir) / "multiqc" / "multiqc_report.html"
                if multiqc.exists():
                    job.multiqc_report = str(multiqc)
            else:
                job.status = "failed"
                job.error_message = (proc.stderr or "")[:2000] if proc.stderr else "Unknown error"

            job.completed_at = datetime.utcnow()
            db.add(job)

    except Exception as e:  # noqa: BLE001
        logger.exception("[NEXTFLOW] Job %s failed", job_id)
        with db_session() as db:
            job = db.query(NextflowJob).filter_by(id=job_id).first()
            if job:
                job.status = "failed"
                job.error_message = str(e)[:2000]
                job.completed_at = datetime.utcnow()
                db.add(job)


def get_nextflow_job_status(job_id: UUID) -> Optional[NextflowJob]:
    """Get job status from DB."""
    with db_session() as db:
        job = db.query(NextflowJob).filter_by(id=job_id).first()
        if job:
            db.expunge(job)
        return job


def list_nextflow_jobs(user_id: Optional[str] = None, limit: int = 20) -> List[NextflowJob]:
    """List recent Nextflow jobs."""
    with db_session() as db:
        q = db.query(NextflowJob).order_by(NextflowJob.created_at.desc())
        if user_id:
            q = q.filter_by(created_by=user_id)
        jobs = q.limit(limit).all()
        for j in jobs:
            db.expunge(j)
        return jobs


def parse_nextflow_progress(log_path: Path) -> int:
    """Parse .nextflow.log for progress percentage.

    Looks for lines like: [12/34] process > NFCORE_RNASEQ:...
    """
    if not log_path.exists():
        return 0

    try:
        completed = 0
        total = 0
        with open(log_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if "] process >" in line:
                    # Parse [completed/total]
                    bracket = line.split("]")[0].strip("[")
                    parts = bracket.split("/")
                    if len(parts) == 2:
                        try:
                            completed = int(parts[0])
                            total = int(parts[1])
                        except ValueError:
                            pass

        if total > 0:
            return int((completed / total) * 100)
        return 0
    except Exception:
        return 0


__all__ = [
    "check_nextflow_installed",
    "list_nfcore_pipelines",
    "generate_sample_sheet",
    "submit_nextflow_job",
    "get_nextflow_job_status",
    "list_nextflow_jobs",
    "parse_nextflow_progress",
]


