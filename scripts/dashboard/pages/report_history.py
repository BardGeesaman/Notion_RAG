from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from typing import Optional

import pandas as pd
import streamlit as st

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import ReportArtifact, ReportSchedule
from amprenta_rag.reports.generator import generate_report
from amprenta_rag.reports.artifact_registry import save_artifact
from amprenta_rag.reports.scheduler import add_report_schedule, remove_report_schedule


def _params_hash(entity_type: str, entity_id: Optional[str], fmt: str) -> str:
    payload = {"entity_type": entity_type, "entity_id": entity_id or "", "format": fmt}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _render_artifacts():
    st.subheader("Recent Report Artifacts")
    with db_session() as db:
        artifacts = (
            db.query(ReportArtifact)
            .order_by(ReportArtifact.created_at.desc())
            .limit(100)
            .all()
        )
    if not artifacts:
        st.info("No report artifacts found.")
        return

    rows = []
    for a in artifacts:
        rows.append(
            {
                "ID": str(a.id),
                "Entity Type": a.entity_type,
                "Entity ID": str(a.entity_id),
                "Format": a.format,
                "Created At": a.created_at.isoformat() if a.created_at else "",
                "File Path": a.file_path,
            }
        )
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)

    with st.expander("Download artifacts"):
        for a in artifacts:
            if a.file_path and os.path.exists(a.file_path):
                with open(a.file_path, "rb") as f:
                    data = f.read()
                st.download_button(
                    label=f"⬇️ Download {a.format.upper()} ({a.entity_type} {a.entity_id})",
                    data=data,
                    file_name=os.path.basename(a.file_path),
                    mime="application/octet-stream",
                    key=f"dl_{a.id}",
                )
            else:
                st.text(f"{a.file_path or 'N/A'} (missing)")


def _render_manual_generate():
    st.subheader("Generate Report")
    with st.form("manual_generate"):
        entity_type = st.selectbox("Entity type", ["program", "dataset", "signature", "feature"])
        entity_id = st.text_input("Entity ID (UUID or identifier)")
        fmt = st.selectbox("Format", ["html", "pdf"])
        submitted = st.form_submit_button("Generate")
        if submitted:
            entity_id_val = entity_id.strip()
            if not entity_id_val:
                st.error("Entity ID is required.")
                return
            try:
                from uuid import UUID as _UUID
                entity_uuid = _UUID(entity_id_val)
            except Exception:
                st.error("Invalid entity ID format (must be UUID).")
                return
            try:
                params = {"entity_type": entity_type, "entity_id": str(entity_uuid), "format": fmt}
                report_path = generate_report(
                    template_name="narrative_report.ipynb",
                    params=params,
                    format=fmt,
                )
                phash = _params_hash(entity_type, str(entity_uuid), fmt)
                artifact = save_artifact(
                    entity_type=entity_type,
                    entity_id=str(entity_uuid),
                    format=fmt,
                    file_path=report_path,
                    params_hash=phash,
                    user_id=None,
                )
                st.success(f"Report generated: {report_path}")
                st.write(f"Artifact ID: {artifact.id}")
                if os.path.exists(report_path):
                    with open(report_path, "rb") as f:
                        data = f.read()
                    st.download_button(
                        label="Download report",
                        data=data,
                        file_name=os.path.basename(report_path),
                        mime="application/octet-stream",
                        key="manual_download",
                    )
            except Exception as exc:
                st.error(f"Failed to generate report: {exc}")


def _render_schedule_management():
    st.subheader("Report Schedules")
    with db_session() as db:
        schedules = db.query(ReportSchedule).order_by(ReportSchedule.created_at.desc()).all()

    if schedules:
        sched_rows = []
        for s in schedules:
            sched_rows.append(
                {
                    "ID": str(s.id),
                    "Name": s.name,
                    "Entity Type": s.entity_type,
                    "Entity ID": str(s.entity_id) if s.entity_id else "All",
                    "Format": s.format,
                    "Cron": s.cron_expression,
                    "Enabled": s.enabled,
                    "Last Run": s.last_run_at.isoformat() if s.last_run_at else "",
                }
            )
        st.dataframe(pd.DataFrame(sched_rows), hide_index=True, use_container_width=True)

        for s in schedules:
            col1, col2 = st.columns([3, 1])
            with col1:
                st.markdown(f"**{s.name}** — {s.entity_type} ({s.entity_id or 'all'}) @ `{s.cron_expression}`")
            with col2:
                new_enabled = st.checkbox("Enabled", value=s.enabled, key=f"toggle_{s.id}")
                if new_enabled != s.enabled:
                    with db_session() as db:
                        db_obj = db.query(ReportSchedule).filter(ReportSchedule.id == s.id).first()
                        if db_obj:
                            db_obj.enabled = new_enabled
                            db.commit()
                    if new_enabled:
                        add_report_schedule(s.id, s.cron_expression)
                    else:
                        remove_report_schedule(s.id)
                    st.experimental_rerun()
    else:
        st.info("No report schedules found.")

    st.markdown("---")
    st.subheader("Add Schedule")
    with st.form("add_schedule"):
        name = st.text_input("Name")
        entity_type = st.selectbox("Entity type", ["program", "dataset", "signature"])
        entity_id_input = st.text_input("Entity ID (leave blank for all)")
        fmt = st.selectbox("Format", ["html", "pdf"], key="sched_fmt")
        cron = st.text_input("Cron expression", value="0 8 * * MON")
        create = st.form_submit_button("Create schedule")
        if create:
            if not name.strip():
                st.error("Name is required.")
            elif not cron.strip():
                st.error("Cron expression is required.")
            else:
                with db_session() as db:
                    sched = ReportSchedule(
                        name=name.strip(),
                        entity_type=entity_type,
                        entity_id=entity_id_input.strip() or None,
                        format=fmt,
                        cron_expression=cron.strip(),
                        enabled=True,
                        created_at=datetime.utcnow(),
                    )
                    db.add(sched)
                    db.commit()
                    db.refresh(sched)
                    add_report_schedule(sched.id, cron.strip())
                    st.success(f"Schedule created: {sched.name}")
                    st.rerun()


def render_report_history_page() -> None:
    st.title("📄 Report History")
    _render_artifacts()
    st.markdown("---")
    _render_manual_generate()
    st.markdown("---")
    _render_schedule_management()

