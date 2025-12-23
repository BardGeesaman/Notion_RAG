"""Nextflow Orchestrator page for the Streamlit dashboard.

Provides:
- Run Pipeline: select nf-core pipeline, upload/build sample sheet, submit job
- Job Monitor: monitor active jobs, progress bars, log viewer
- Results: view completed jobs, render MultiQC report, download artifacts
"""

from __future__ import annotations

import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import streamlit as st
import streamlit.components.v1 as components

from amprenta_rag.auth.session import get_current_user
from amprenta_rag.ingestion.genomics.nextflow_runner import (
    check_nextflow_installed,
    generate_sample_sheet,
    list_nfcore_pipelines,
    list_nextflow_jobs,
    parse_nextflow_progress,
    submit_nextflow_job,
)


def _user_id() -> str:
    user = get_current_user() or {}
    return str(user.get("id") or "00000000-0000-0000-0000-000000000001")


def _ensure_dir(primary: Path, fallback: Path) -> Path:
    try:
        primary.mkdir(parents=True, exist_ok=True)
        return primary
    except Exception:
        fallback.mkdir(parents=True, exist_ok=True)
        return fallback


RUNS_DIR = _ensure_dir(Path("/data/nextflow/runs"), Path("data/nextflow/runs"))


def _init_state() -> None:
    st.session_state.setdefault("nf_selected_job_id", None)
    st.session_state.setdefault("nf_uploaded_sample_sheet", None)
    st.session_state.setdefault("nf_built_sample_sheet_path", None)


def _render_nextflow_missing(info: Dict[str, Any]) -> None:
    st.error("Nextflow is not installed or not available in this environment.")
    if info.get("error"):
        st.caption(f"Details: {info['error']}")
    st.markdown("**Install hints:**")
    st.code(
        "\n".join(
            [
                "conda install -c conda-forge openjdk=17",
                "curl -s https://get.nextflow.io | bash",
                "sudo mv nextflow /usr/local/bin/",
            ]
        )
    )


def _save_upload(uploaded, out_dir: Path) -> Optional[Path]:
    if uploaded is None:
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / uploaded.name
    try:
        p.write_bytes(uploaded.getvalue())
        return p
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to save upload: {e}")
        return None


def _tab_run_pipeline() -> None:
    st.subheader("Run Pipeline")

    nf_info = check_nextflow_installed()
    if not nf_info.get("installed"):
        _render_nextflow_missing(nf_info)
        return

    st.success(f"Nextflow available: {nf_info.get('version')}")
    if nf_info.get("java_version"):
        st.caption(f"Java: {nf_info.get('java_version')}")

    pipelines = list_nfcore_pipelines()
    if not pipelines:
        st.warning("No pipelines configured.")
        return

    st.markdown("### 1) Pipeline")
    labels = [f"{p['name']} ({p.get('version','')}) - {p.get('description','')}" for p in pipelines]
    sel = st.selectbox("Pipeline", labels, index=0)
    pipeline = pipelines[labels.index(sel)]["name"]
    pipeline_version = pipelines[labels.index(sel)].get("version", "latest")

    st.markdown("### 2) Sample Sheet")
    ss_tab1, ss_tab2 = st.tabs(["Upload CSV", "Build CSV"])

    run_dir = RUNS_DIR / datetime.utcnow().strftime("%Y%m%d_%H%M%S")

    sample_sheet_path: Optional[Path] = None
    with ss_tab1:
        up = st.file_uploader("Upload nf-core sample sheet CSV", type=["csv"])
        if up is not None:
            p = _save_upload(up, run_dir)
            if p:
                st.session_state["nf_uploaded_sample_sheet"] = str(p)
                st.success(f"Saved sample sheet: {p}")

        cur = st.session_state.get("nf_uploaded_sample_sheet")
        if cur:
            st.caption(f"Using sample sheet: {cur}")
            sample_sheet_path = Path(cur)

    with ss_tab2:
        st.caption("Provide FASTQ file paths; strandedness defaults to 'auto'.")
        rows = st.data_editor(
            [
                {"sample": "sample1", "fastq_1": "/path/to/R1.fastq.gz", "fastq_2": "", "strandedness": "auto"},
            ],
            num_rows="dynamic",
            use_container_width=True,
            key="nf_sample_sheet_rows",
        )
        out_name = st.text_input("Output filename", value="sample_sheet.csv")
        if st.button("Generate sample sheet", type="secondary"):
            try:
                out_path = run_dir / out_name
                p2 = generate_sample_sheet(fastq_files=rows, output_path=out_path)  # type: ignore[arg-type]
                st.session_state["nf_built_sample_sheet_path"] = str(p2)
                st.success(f"Generated: {p2}")
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to generate sample sheet: {e}")

        built = st.session_state.get("nf_built_sample_sheet_path")
        if built:
            st.caption(f"Built sample sheet: {built}")
            sample_sheet_path = Path(built)

    st.markdown("### 3) Parameters")
    genome = st.selectbox("Genome", ["GRCh38", "GRCm39", "GRCh37"], index=0)
    profile = st.selectbox("Profile", ["docker", "conda", "singularity"], index=0)
    resume = st.checkbox("Resume", value=False)

    out_dir = st.text_input("Output directory", value=str(run_dir))

    st.markdown("### 4) Execute")
    disabled = sample_sheet_path is None or not sample_sheet_path.exists()
    if st.button("Submit Nextflow Job", type="primary", disabled=disabled, use_container_width=True):
        if sample_sheet_path is None or not sample_sheet_path.exists():
            st.error("Please provide a valid sample sheet.")
        else:
            params: Dict[str, Any] = {"profile": profile, "resume": resume, "version": pipeline_version}
            try:
                job_id = submit_nextflow_job(
                    pipeline=pipeline,
                    sample_sheet=sample_sheet_path,
                    genome=genome,
                    output_dir=Path(out_dir),
                    user_id=_user_id(),
                    params=params,
                    profile=profile,
                )
                st.session_state["nf_selected_job_id"] = str(job_id)
                st.success(f"Submitted job: {job_id}")
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to submit job: {e}")


def _tail_file(path: Path, max_chars: int = 8000) -> str:
    try:
        data = path.read_text(encoding="utf-8", errors="replace")
        if len(data) > max_chars:
            return data[-max_chars:]
        return data
    except Exception:
        return ""


def _tab_job_monitor() -> None:
    st.subheader("Job Monitor")

    nf_info = check_nextflow_installed()
    if not nf_info.get("installed"):
        _render_nextflow_missing(nf_info)
        return

    jobs = list_nextflow_jobs(user_id=_user_id(), limit=20)
    active = [j for j in jobs if (j.status or "").lower() in ("pending", "running")]
    if not active:
        st.info("No active jobs.")
    for job in active:
        title = f"{str(job.id)[:8]} | {job.pipeline_name} | {job.status} | {job.progress_percent or 0}%"
        with st.expander(title, expanded=False):
            pct = int(job.progress_percent or 0)
            # Try to parse progress from log file if available.
            try:
                if job.nextflow_log:
                    pct = max(pct, parse_nextflow_progress(Path(job.nextflow_log)))
            except Exception:
                pass
            st.progress(max(0, min(100, pct)))
            st.write(
                {
                    "id": str(job.id),
                    "pipeline": job.pipeline_name,
                    "status": job.status,
                    "genome": job.genome,
                    "output_dir": job.output_dir,
                    "work_dir": job.work_dir,
                    "nextflow_log": job.nextflow_log,
                    "multiqc_report": job.multiqc_report,
                    "error_message": job.error_message,
                    "started_at": str(job.started_at) if job.started_at else None,
                }
            )
            if job.nextflow_log:
                log_tail = _tail_file(Path(job.nextflow_log))
                if log_tail:
                    st.markdown("**.nextflow.log (tail)**")
                    st.code(log_tail, language="text")

    c1, c2 = st.columns(2)
    with c1:
        if st.button("Refresh"):
            st.rerun()
    with c2:
        auto = st.checkbox("Auto-refresh (2s)", value=False, key="nf_autorefresh")
    if auto:
        time.sleep(2.0)
        st.rerun()


def _tab_results() -> None:
    st.subheader("Results")

    nf_info = check_nextflow_installed()
    if not nf_info.get("installed"):
        _render_nextflow_missing(nf_info)
        return

    jobs = list_nextflow_jobs(user_id=_user_id(), limit=20)
    done = [j for j in jobs if (j.status or "").lower() in ("complete", "failed")]
    if not done:
        st.info("No completed jobs yet.")
        return

    opts = [f"{str(j.id)[:8]} | {j.pipeline_name} | {j.status}" for j in done]
    sel = st.selectbox("Select job", opts, index=0)
    job = done[opts.index(sel)]

    st.write(
        {
            "id": str(job.id),
            "status": job.status,
            "pipeline": job.pipeline_name,
            "output_dir": job.output_dir,
            "multiqc_report": job.multiqc_report,
            "error_message": job.error_message,
            "completed_at": str(job.completed_at) if job.completed_at else None,
        }
    )
    if job.error_message and (job.status or "").lower() == "failed":
        st.error(job.error_message)

    if job.multiqc_report:
        p = Path(job.multiqc_report)
        if p.exists():
            st.markdown("### MultiQC Report")
            try:
                html = p.read_text(encoding="utf-8", errors="replace")
                components.html(html, height=800, scrolling=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to render MultiQC report: {e}")

            try:
                st.download_button(
                    "Download MultiQC report",
                    data=p.read_bytes(),
                    file_name=p.name,
                    mime="text/html",
                    use_container_width=True,
                )
            except Exception:
                pass

    if job.nextflow_log:
        p2 = Path(job.nextflow_log)
        if p2.exists():
            st.markdown("### Nextflow Log")
            st.code(_tail_file(p2), language="text")
            try:
                st.download_button(
                    "Download .nextflow.log",
                    data=p2.read_bytes(),
                    file_name=p2.name,
                    mime="text/plain",
                    use_container_width=True,
                )
            except Exception:
                pass

    st.divider()
    st.subheader("Ingest Results")
    if st.button("Ingest as Dataset", key="ingest_nextflow", use_container_width=True):
        try:
            from amprenta_rag.database.models import Dataset
            from amprenta_rag.database.session import db_session

            out_dir = Path(job.output_dir or "")
            if not out_dir:
                st.error("Job output directory is missing.")
                return

            with db_session() as db:
                ds = Dataset(
                    name=f"Nextflow {job.pipeline_name} {str(job.id)[:8]}",
                    omics_type="transcriptomics",
                    description=f"Nextflow results for job {job.id} ({job.pipeline_name}).",
                    file_paths=[str(out_dir)],
                    dataset_source_type="Nextflow Orchestrator",
                    data_origin="Internal – Pipeline",
                )
                db.add(ds)
                db.flush()
                db.refresh(ds)
                st.success("Dataset created from Nextflow results")
                st.caption(f"Dataset ID: {ds.id}")
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to create dataset: {e}")


def render_nextflow_orchestrator_page() -> None:
    """Render the Nextflow Orchestrator page."""
    from scripts.dashboard.auth import require_auth

    require_auth()

    _init_state()

    st.header("🧪 Nextflow Orchestrator")
    st.caption("Submit and monitor nf-core Nextflow pipelines.")

    tab1, tab2, tab3 = st.tabs(["Run Pipeline", "Job Monitor", "Results"])
    with tab1:
        _tab_run_pipeline()
    with tab2:
        _tab_job_monitor()
    with tab3:
        _tab_results()


__all__ = ["render_nextflow_orchestrator_page"]


