"""Security headers component for Streamlit."""
import streamlit.components.v1 as components


def inject_security_meta_tags():
    """Inject security-related meta tags into the page."""
    meta_tags = """
    <meta http-equiv="X-Content-Type-Options" content="nosniff">
    <meta http-equiv="X-Frame-Options" content="SAMEORIGIN">
    <meta name="referrer" content="strict-origin-when-cross-origin">
    """
    components.html(meta_tags, height=0, width=0)
