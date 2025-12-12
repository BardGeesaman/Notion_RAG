"""System Health page for monitoring platform status."""
from __future__ import annotations

import streamlit as st

from amprenta_rag.utils.health import get_db_stats, check_api_connectivity, get_system_info
from scripts.dashboard.db_session import db_session
from amprenta_rag.auth.session import get_current_user


def render_system_health_page() -> None:
    """Render the System Health page (admin only)."""
    user = get_current_user()
    
    # Admin only
    if not user or user.get("role") != "admin":
        st.error("Access denied. Only administrators can view system health.")
        return
    
    st.header("🏥 System Health")
    st.markdown("Monitor platform status, database statistics, and API connectivity.")
    
    with db_session() as db:
        # Database Statistics
        st.subheader("📊 Database Statistics")
        stats = get_db_stats(db)
        
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            st.metric("Users", stats["users"])
        with col2:
            st.metric("Experiments", stats["experiments"])
        with col3:
            st.metric("Compounds", stats["compounds"])
        with col4:
            st.metric("Signatures", stats["signatures"])
        with col5:
            st.metric("Datasets", stats["datasets"])
        
        st.markdown("---")
        
        # API Connectivity
        st.subheader("🔌 API Connectivity")
        with st.spinner("Checking API connectivity..."):
            api_status = check_api_connectivity()
        
        col1, col2 = st.columns(2)
        with col1:
            if api_status["openai"]:
                st.success("✅ OpenAI API: Connected")
            else:
                st.error("❌ OpenAI API: Not connected")
        
        with col2:
            if api_status["pinecone"]:
                st.success("✅ Pinecone API: Connected")
            else:
                st.error("❌ Pinecone API: Not connected")
        
        st.markdown("---")
        
        # System Information
        st.subheader("💻 System Information")
        sys_info = get_system_info()
        
        st.markdown(f"**Python Version:** `{sys_info['python_version'].split()[0]}`")
        st.markdown(f"**Platform:** `{sys_info['platform']}`")
        st.markdown(f"**Processor:** `{sys_info['processor'] or 'Unknown'}`")
        st.markdown(f"**Architecture:** `{sys_info['architecture']}`")
        
        if sys_info["memory"].get("available"):
            st.markdown("**Memory:**")
            mem = sys_info["memory"]
            st.markdown(f"- Total: {mem['total_gb']} GB")
            st.markdown(f"- Available: {mem['available_gb']} GB")
            st.markdown(f"- Used: {mem['used_gb']} GB ({mem['percent']}%)")
        else:
            st.info("Memory information not available (psutil not installed)")


if __name__ == "__main__":
    render_system_health_page()
