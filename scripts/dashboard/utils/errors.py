import traceback
from typing import Callable, Any, Optional

import streamlit as st
import httpx
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


def render_error_with_action(
    title: str,
    message: str,
    action_label: str = "Retry",
    action_callback: Optional[Callable] = None
):
    """Render error with actionable button.
    
    Args:
        title: Error title to display
        message: Detailed error message
        action_label: Label for the action button
        action_callback: Function to call when button is clicked
    """
    st.error(f"❌ {title}")
    st.info(message)
    if action_callback:
        if st.button(action_label, key=f"error_action_{hash(title)}"):
            action_callback()
            st.rerun()


def render_connection_error(service: str = "API"):
    """Standardized connection error display.
    
    Args:
        service: Name of the service that couldn't be reached
    """
    st.warning(f"⚠️ Unable to connect to {service}")
    st.info("Please check that the backend server is running and try again.")
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("🔄 Retry", key=f"retry_{service}"):
            st.rerun()
    with col2:
        error_details = f"Connection failed to {service} at {st.session_state.get('last_error_time', 'unknown time')}"
        st.button("📋 Copy Error Details", 
                 disabled=True, 
                 help=error_details,
                 key=f"copy_error_{service}")


def render_empty_state(
    entity: str,
    message: Optional[str] = None,
    action_label: Optional[str] = None,
    action_callback: Optional[Callable] = None
):
    """Render empty state for no data scenarios.
    
    Args:
        entity: Type of entity (e.g., "datasets", "experiments")
        message: Custom message to display
        action_label: Label for the action button
        action_callback: Function to call when button is clicked
    """
    default_message = f"No {entity} found."
    st.info(f"ℹ️ {message or default_message}")
    
    if action_label and action_callback:
        if st.button(action_label, key=f"empty_action_{entity}"):
            action_callback()


def render_loading_error(
    entity: str,
    error: Exception,
    retry_callback: Optional[Callable] = None
):
    """Render loading error with retry option.
    
    Args:
        entity: Type of entity that failed to load
        error: The exception that occurred
        retry_callback: Function to call for retry
    """
    st.error(f"❌ Failed to load {entity}")
    
    # Provide specific error context based on exception type
    if isinstance(error, httpx.ConnectError):
        st.info("The server appears to be unavailable. Please check the connection and try again.")
    elif isinstance(error, httpx.TimeoutException):
        st.info("The request timed out. The server may be experiencing high load.")
    elif isinstance(error, httpx.HTTPStatusError):
        st.info(f"Server returned error {error.response.status_code}. Please try again later.")
    else:
        st.info(f"An unexpected error occurred: {str(error)}")
    
    if retry_callback:
        if st.button("🔄 Retry", key=f"retry_loading_{entity}"):
            retry_callback()
            st.rerun()


def safe_api_call(
    func: Callable, 
    fallback: Any = None, 
    error_message: str = "Operation failed",
    show_error: bool = True
):
    """Wrapper for API calls with error handling.
    
    Args:
        func: Function to execute
        fallback: Value to return on error
        error_message: Error message prefix
        show_error: Whether to display error in UI
    
    Returns:
        Function result or fallback value on error
    """
    try:
        return func()
    except httpx.ConnectError:
        if show_error:
            render_connection_error()
        return fallback
    except httpx.HTTPStatusError as e:
        if show_error:
            st.error(f"❌ {error_message}: HTTP {e.response.status_code}")
        return fallback
    except httpx.TimeoutException:
        if show_error:
            st.error(f"❌ {error_message}: Request timed out")
        return fallback
    except Exception as e:
        if show_error:
            st.error(f"❌ {error_message}: {str(e)}")
        return fallback


def with_error_handling(
    error_message: str = "Operation failed",
    fallback: Any = None,
    show_traceback: bool = False
):
    """Decorator for functions that need error handling.
    
    Args:
        error_message: Message to display on error
        fallback: Value to return on error
        show_traceback: Whether to show full traceback
    
    Usage:
        @with_error_handling("Failed to load data", fallback=[])
        def load_data():
            return api_call()
    """
    def decorator(func: Callable) -> Callable:
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                st.error(f"❌ {error_message}: {str(e)}")
                if show_traceback:
                    st.code(traceback.format_exc())
                return fallback
        return wrapper
    return decorator


def render_validation_error(field: str, message: str):
    """Render validation error for form fields.
    
    Args:
        field: Name of the field with validation error
        message: Validation error message
    """
    st.error(f"❌ {field}: {message}")


def render_success_message(message: str, details: Optional[str] = None):
    """Render success message with optional details.
    
    Args:
        message: Success message
        details: Optional additional details
    """
    st.success(f"✅ {message}")
    if details:
        st.info(details)


def render_warning_message(message: str, details: Optional[str] = None):
    """Render warning message with optional details.
    
    Args:
        message: Warning message
        details: Optional additional details
    """
    st.warning(f"⚠️ {message}")
    if details:
        st.info(details)
