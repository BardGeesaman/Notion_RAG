import traceback

import streamlit as st


def render_ingest_error(exc: Exception):
    msg = str(exc)
    if any(x in msg.lower() for x in ["no feature column", "missing required column"]):
        st.error("❌ Ingestion failed: Could not find the necessary columns in your file.")
        st.info(
            "Check that required columns exist and are correctly named. See the ingestion format guide for details."
        )
    elif "parse" in msg.lower() or "delimiter" in msg.lower():
        st.error("❌ Ingestion failed: File could not be parsed as CSV/TSV. ")
        st.info("Check the file format, encoding, and delimiter.")
    elif "validation" in msg.lower():
        st.error(f"❌ Ingestion failed: {msg}")
        st.info("Please review all required metadata fields for validity.")
    else:
        st.error(f"❌ Ingestion failed: {msg}")
        st.info("See server logs for details or contact support.")
    # (Optionally: log full trace for debugging, stdout or to server logs)
    traceback.print_exc()
