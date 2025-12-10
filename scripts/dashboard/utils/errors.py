import traceback

import streamlit as st
from requests.exceptions import ConnectionError, HTTPError, Timeout
from sqlalchemy.exc import IntegrityError, OperationalError, ProgrammingError


def render_ingest_error(exc: Exception):
    msg = str(exc)
    if isinstance(exc, OperationalError):
        st.error("❌ Database connection failed.")
        st.info("Check DB_URL, credentials, and that Postgres is running.")
    elif isinstance(exc, ProgrammingError):
        st.error("❌ Database schema error during ingestion.")
        st.info("Run migrations (alembic upgrade head) and retry.")
    elif isinstance(exc, IntegrityError):
        st.error("❌ Duplicate or constraint violation detected.")
        st.info("Check for existing records or conflicting keys before re-ingesting.")
    elif isinstance(exc, (ConnectionError, Timeout)):
        st.error("❌ Network/API unreachable during ingestion.")
        st.info("Verify internet connectivity and upstream API status, then retry.")
    elif any(x in msg.lower() for x in ["no feature column", "missing required column"]):
        st.error("❌ Ingestion failed: Could not find the necessary columns in your file.")
        st.info(
            "Check that required columns exist and are correctly named. See the ingestion format guide for details."
        )
    elif "parse" in msg.lower() or "delimiter" in msg.lower():
        st.error("❌ Ingestion failed: File could not be parsed as CSV/TSV.")
        st.info("Check the file format, encoding, and delimiter.")
    elif "validation" in msg.lower():
        st.error(f"❌ Ingestion failed: {msg}")
        st.info("Please review all required metadata fields for validity.")
    else:
        st.error(f"❌ Ingestion failed: {msg}")
        st.info("See server logs for details or contact support.")
    traceback.print_exc()


def render_api_error(exc: Exception, service: str = "API"):
    if isinstance(exc, (ConnectionError, Timeout)):
        st.error(f"❌ {service} unreachable.")
        st.info("Check internet connectivity, VPN/proxy settings, and service status page.")
    elif isinstance(exc, HTTPError):
        st.error(f"❌ {service} returned an HTTP error: {exc}")
        st.info("Inspect the status code and response; ensure credentials and payload are valid.")
    else:
        st.error(f"❌ {service} request failed: {exc}")
        st.info("Retry the request or contact support if the issue persists.")
    traceback.print_exc()
