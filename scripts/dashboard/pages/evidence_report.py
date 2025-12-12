import streamlit as st

from amprenta_rag.database.models import Dataset, Program, Signature
from amprenta_rag.reporting.evidence_report import (
    generate_dataset_report,
    generate_program_report,
    generate_signature_report,
)
from amprenta_rag.reporting.pdf_export import export_evidence_report_to_pdf
from scripts.dashboard.db_session import db_session


def render_evidence_report_page():
    st.header("Evidence Report Viewer")
    entity_type = st.selectbox("Select Entity Type", ["Program", "Dataset", "Signature"])
    
    with db_session() as db:
        if entity_type == "Program":
            options = [(p.id, p.name) for p in db.query(Program).all()]
        elif entity_type == "Dataset":
            options = [(d.id, d.name) for d in db.query(Dataset).all()]
        else:
            options = [(s.id, s.name) for s in db.query(Signature).all()]
    
    selected_id = st.selectbox("Select Entity", options, format_func=lambda x: x[1])
    if st.button("Generate Report"):
        if entity_type == "Program":
            rep = generate_program_report(selected_id[0])
        elif entity_type == "Dataset":
            rep = generate_dataset_report(selected_id[0])
        else:
            rep = generate_signature_report(selected_id[0])
        md = f"# Evidence Report: {entity_type} {rep.entity_id}\n_Generated at: {rep.generated_at}_\n\n"
        for section in rep.sections:
            md += f"## {section.title}\n\n{section.summary_text}\n\n"
            if section.supporting_datasets:
                md += f"**Supporting Datasets:** {section.supporting_datasets}\n\n"
            if section.key_features:
                md += f"**Key Features:** {section.key_features}\n\n"
            if section.signatures:
                md += f"**Signatures:** {section.signatures}\n\n"
            if section.references:
                md += f"**References:** {section.references}\n\n"
        st.markdown(md)

        # Download buttons
        col1, col2 = st.columns(2)
        fname_base = f"evidence_report_{entity_type.lower()}_{rep.entity_id}"

        with col1:
            st.download_button(
                "📄 Download Markdown",
                md,
                file_name=f"{fname_base}.md",
                mime="text/markdown",
            )

        with col2:
            try:
                pdf_bytes = export_evidence_report_to_pdf(rep)
                st.download_button(
                    "📕 Download PDF",
                    pdf_bytes,
                    file_name=f"{fname_base}.pdf",
                    mime="application/pdf",
                )
            except Exception as e:
                st.warning(f"PDF export unavailable: {e}")
