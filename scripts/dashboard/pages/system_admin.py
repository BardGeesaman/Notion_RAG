"""System Administration page."""

from __future__ import annotations

import os
from typing import Any

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: dict, *, timeout: int = 30) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_system_admin_page() -> None:
    """Render the System Administration page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("⚙️ System Administration")
    st.caption("Monitor system health, manage caches, and view queue status.")
    
    tab1, tab2, tab3, tab4 = st.tabs(["System Health", "Cache Management", "Queue Health", "Connections"])
    
    with tab1:
        render_system_health_tab()
    
    with tab2:
        render_cache_management_tab()
    
    with tab3:
        render_queue_health_tab()
    
    with tab4:
        render_connections_tab()


def render_system_health_tab() -> None:
    """Render system health monitoring tab."""
    st.subheader("System Health")
    
    # Admin check
    if st.session_state.get("user", {}).get("role") != "admin":
        st.error("⚠️ Admin access required")
        return
    
    col1, col2 = st.columns([3, 1])
    
    with col2:
        if st.button("Refresh", type="secondary"):
            # Clear cached data to force refresh
            if "system_health" in st.session_state:
                del st.session_state["system_health"]
    
    if st.button("Load System Health", type="primary") or "system_health" in st.session_state:
        try:
            if "system_health" not in st.session_state:
                health_data = _api_get("/api/v1/admin/health/system")
                st.session_state["system_health"] = health_data
            else:
                health_data = st.session_state["system_health"]
            
            # Display metrics in columns
            col1, col2, col3 = st.columns(3)
            
            with col1:
                cpu_usage = health_data.get("cpu_percent", 0)
                st.metric("CPU Usage", f"{cpu_usage:.1f}%")
                if cpu_usage > 80:
                    st.error("⚠️ High CPU usage")
                elif cpu_usage > 60:
                    st.warning("⚠️ Moderate CPU usage")
                else:
                    st.success("✅ CPU usage normal")
            
            with col2:
                memory_data = health_data.get("memory", {})
                memory_used_gb = memory_data.get("used", 0) / (1024**3)  # Convert to GB
                memory_total_gb = memory_data.get("total", 0) / (1024**3)
                memory_percent = memory_data.get("percent", 0)
                
                st.metric("Memory Usage", f"{memory_used_gb:.1f}GB / {memory_total_gb:.1f}GB")
                st.metric("Memory %", f"{memory_percent:.1f}%")
                
                if memory_percent > 85:
                    st.error("⚠️ High memory usage")
                elif memory_percent > 70:
                    st.warning("⚠️ Moderate memory usage")
                else:
                    st.success("✅ Memory usage normal")
            
            with col3:
                disk_data = health_data.get("disk", {})
                disk_used_gb = disk_data.get("used", 0) / (1024**3)
                disk_total_gb = disk_data.get("total", 0) / (1024**3)
                disk_percent = disk_data.get("percent", 0)
                
                st.metric("Disk Usage", f"{disk_used_gb:.1f}GB / {disk_total_gb:.1f}GB")
                st.metric("Disk %", f"{disk_percent:.1f}%")
                
                if disk_percent > 90:
                    st.error("⚠️ High disk usage")
                elif disk_percent > 75:
                    st.warning("⚠️ Moderate disk usage")
                else:
                    st.success("✅ Disk usage normal")
            
            # Timestamp
            timestamp = health_data.get("timestamp", "Unknown")
            st.caption(f"Last updated: {timestamp}")
            
        except Exception as e:
            st.error(f"Failed to load system health: {e}")


def render_cache_management_tab() -> None:
    """Render cache management tab."""
    st.subheader("Cache Management")
    
    # Admin check
    if st.session_state.get("user", {}).get("role") != "admin":
        st.error("⚠️ Admin access required")
        return
    
    if st.button("Load Cache Stats", type="primary"):
        try:
            # Get cache summary
            summary = _api_get("/api/v1/admin/caches/summary")
            st.session_state["cache_summary"] = summary
            
            # Get detailed cache stats
            caches = _api_get("/api/v1/admin/caches")
            st.session_state["cache_stats"] = caches
            
        except Exception as e:
            st.error(f"Failed to load cache stats: {e}")
    
    # Display cache summary
    summary = st.session_state.get("cache_summary")
    if summary:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Total Caches", summary.get("total_caches", 0))
        
        with col2:
            st.metric("Total Entries", summary.get("total_entries", 0))
        
        with col3:
            error_count = summary.get("caches_with_errors", 0)
            st.metric("Caches with Errors", error_count)
            if error_count > 0:
                st.warning(f"⚠️ {error_count} cache(s) have errors")
    
    # Display cache stats table
    cache_stats = st.session_state.get("cache_stats")
    if cache_stats:
        st.markdown("### Cache Statistics")
        
        # Convert to DataFrame
        df_data = []
        for cache_name, stats in cache_stats.items():
            df_data.append({
                "Cache": cache_name,
                "Entries": stats.get("entries", 0),
                "Hits": stats.get("hits", 0),
                "Misses": stats.get("misses", 0),
                "Evictions": stats.get("evictions", 0),
                "TTL": stats.get("ttl", "N/A"),
            })
        
        if df_data:
            df = pd.DataFrame(df_data)
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            # Cache management buttons
            st.markdown("### Cache Actions")
            
            col1, col2 = st.columns(2)
            
            with col1:
                cache_to_clear = st.selectbox("Select cache to clear", [row["Cache"] for row in df_data])
                
                if st.button("Clear Selected Cache", type="secondary"):
                    try:
                        _api_post(f"/api/v1/admin/caches/{cache_to_clear}/clear", {})
                        st.success(f"✅ Cache '{cache_to_clear}' cleared successfully")
                        # Clear cached data to force refresh
                        if "cache_stats" in st.session_state:
                            del st.session_state["cache_stats"]
                        if "cache_summary" in st.session_state:
                            del st.session_state["cache_summary"]
                    except Exception as e:
                        st.error(f"Failed to clear cache: {e}")
            
            with col2:
                st.markdown("**Clear All Caches**")
                confirm_clear_all = st.checkbox("I confirm clearing ALL caches")
                
                clear_all_button = st.button(
                    "Clear All Caches", 
                    type="primary",
                    disabled=not confirm_clear_all
                )
                
                if clear_all_button:
                    try:
                        _api_post("/api/v1/admin/caches/clear-all", {})
                        st.success("✅ All caches cleared successfully")
                        # Clear cached data
                        if "cache_stats" in st.session_state:
                            del st.session_state["cache_stats"]
                        if "cache_summary" in st.session_state:
                            del st.session_state["cache_summary"]
                    except Exception as e:
                        st.error(f"Failed to clear all caches: {e}")


def render_queue_health_tab() -> None:
    """Render queue health monitoring tab."""
    st.subheader("Queue Health")
    
    # Admin check
    if st.session_state.get("user", {}).get("role") != "admin":
        st.error("⚠️ Admin access required")
        return
    
    if st.button("Load Queue Status", type="primary"):
        try:
            queue_data = _api_get("/api/v1/admin/health/queues")
            st.session_state["queue_health"] = queue_data
        except Exception as e:
            st.error(f"Failed to load queue status: {e}")
    
    queue_data = st.session_state.get("queue_health")
    if queue_data:
        # Task counts
        st.markdown("### Task Counts")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            active_tasks = queue_data.get("active_tasks", 0)
            st.metric("Active Tasks", active_tasks)
        
        with col2:
            reserved_tasks = queue_data.get("reserved_tasks", 0)
            st.metric("Reserved Tasks", reserved_tasks)
        
        with col3:
            scheduled_tasks = queue_data.get("scheduled_tasks", 0)
            st.metric("Scheduled Tasks", scheduled_tasks)
        
        # Worker information
        st.markdown("### Workers")
        workers = queue_data.get("workers", [])
        worker_count = len(workers)
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Worker Count", worker_count)
        
        with col2:
            queue_status = queue_data.get("queue_status", "unknown")
            if queue_status == "connected":
                st.success("✅ Queue connected")
            else:
                st.error(f"❌ Queue status: {queue_status}")
        
        # Worker details
        if workers:
            worker_data = []
            for worker in workers:
                worker_data.append({
                    "Worker ID": worker.get("id", "Unknown"),
                    "Status": worker.get("status", "Unknown"),
                    "Last Seen": worker.get("last_seen", "Unknown"),
                })
            
            df = pd.DataFrame(worker_data)
            st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.info("No workers currently active")


def render_connections_tab() -> None:
    """Render connection status tab."""
    st.subheader("Connection Status")
    
    # Admin check
    if st.session_state.get("user", {}).get("role") != "admin":
        st.error("⚠️ Admin access required")
        return
    
    if st.button("Check Connections", type="primary"):
        try:
            conn_data = _api_get("/api/v1/admin/health/connections")
            st.session_state["connection_health"] = conn_data
        except Exception as e:
            st.error(f"Failed to check connections: {e}")
    
    conn_data = st.session_state.get("connection_health")
    if conn_data:
        # PostgreSQL connection
        st.markdown("### Database Connections")
        
        postgres_status = conn_data.get("postgresql", {})
        postgres_connected = postgres_status.get("connected", False)
        postgres_error = postgres_status.get("error")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**PostgreSQL**")
            if postgres_connected:
                st.success("✅ Connected")
            else:
                st.error("❌ Disconnected")
                if postgres_error:
                    st.error(f"Error: {postgres_error}")
        
        # Redis connection
        redis_status = conn_data.get("redis", {})
        redis_connected = redis_status.get("connected", False)
        redis_error = redis_status.get("error")
        
        with col2:
            st.markdown("**Redis**")
            if redis_connected:
                st.success("✅ Connected")
            else:
                st.error("❌ Disconnected")
                if redis_error:
                    st.error(f"Error: {redis_error}")
        
        # Connection timestamps
        st.markdown("### Connection Details")
        postgres_timestamp = postgres_status.get("timestamp", "Unknown")
        redis_timestamp = redis_status.get("timestamp", "Unknown")
        
        col1, col2 = st.columns(2)
        with col1:
            st.caption(f"PostgreSQL last checked: {postgres_timestamp}")
        
        with col2:
            st.caption(f"Redis last checked: {redis_timestamp}")


if __name__ == "__main__":
    render_system_admin_page()
