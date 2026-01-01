"""Streamlit caching utilities with cache invalidation support."""

import streamlit as st
import httpx
from typing import Any, Optional
from functools import wraps

# Default TTL for cached data (5 minutes)
DEFAULT_TTL = 300


def get_api_base() -> str:
    """Get API base URL from environment or default."""
    import os
    return os.environ.get("API_BASE", "http://localhost:8000")


@st.cache_data(ttl=DEFAULT_TTL)
def fetch_programs() -> list:
    """Fetch programs list with caching.
    
    Cache invalidation: Call clear_all_caches() after program mutations.
    """
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{get_api_base()}/api/v1/programs")
            response.raise_for_status()
            return response.json()
    except Exception:
        return []


@st.cache_data(ttl=DEFAULT_TTL)
def fetch_datasets(limit: int = 50) -> list:
    """Fetch datasets list with caching.
    
    Cache invalidation: Call clear_all_caches() after dataset mutations.
    """
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{get_api_base()}/api/v1/datasets", params={"limit": limit})
            response.raise_for_status()
            return response.json()
    except Exception:
        return []


@st.cache_data(ttl=DEFAULT_TTL)
def fetch_experiments(limit: int = 50) -> list:
    """Fetch experiments list with caching.
    
    Cache invalidation: Call clear_all_caches() after experiment mutations.
    """
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{get_api_base()}/api/v1/experiments", params={"limit": limit})
            response.raise_for_status()
            return response.json()
    except Exception:
        return []


@st.cache_data(ttl=DEFAULT_TTL)
def fetch_compounds(limit: int = 50) -> list:
    """Fetch compounds list with caching.
    
    Cache invalidation: Call clear_all_caches() after compound mutations.
    """
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{get_api_base()}/api/v1/compounds", params={"limit": limit})
            response.raise_for_status()
            return response.json()
    except Exception:
        return []


@st.cache_data(ttl=DEFAULT_TTL)
def fetch_structures(limit: int = 50) -> list:
    """Fetch protein structures list with caching.
    
    Cache invalidation: Call clear_all_caches() after structure mutations.
    """
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{get_api_base()}/api/v1/structures", params={"limit": limit})
            response.raise_for_status()
            return response.json()
    except Exception:
        return []


@st.cache_data(ttl=DEFAULT_TTL)
def fetch_docking_runs(limit: int = 50) -> list:
    """Fetch docking runs list with caching.
    
    Cache invalidation: Call clear_all_caches() after docking mutations.
    """
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{get_api_base()}/api/v1/docking/runs", params={"limit": limit})
            response.raise_for_status()
            return response.json()
    except Exception:
        return []


@st.cache_data(ttl=DEFAULT_TTL)
def fetch_jobs(limit: int = 50) -> list:
    """Fetch jobs list with caching.
    
    Cache invalidation: Call clear_all_caches() after job mutations.
    """
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{get_api_base()}/api/v1/jobs", params={"limit": limit})
            response.raise_for_status()
            return response.json().get("jobs", [])
    except Exception:
        return []


def clear_all_caches():
    """Clear all Streamlit caches. Call after mutations.
    
    This function should be called after any create/update/delete operations
    to ensure fresh data is fetched on the next request.
    """
    fetch_programs.clear()
    fetch_datasets.clear()
    fetch_experiments.clear()
    fetch_compounds.clear()
    fetch_structures.clear()
    fetch_docking_runs.clear()
    fetch_jobs.clear()
    st.cache_data.clear()


def clear_specific_cache(cache_name: str):
    """Clear a specific cache by name.
    
    Args:
        cache_name: Name of the cache function (e.g., 'programs', 'datasets')
    """
    cache_map = {
        'programs': fetch_programs,
        'datasets': fetch_datasets,
        'experiments': fetch_experiments,
        'compounds': fetch_compounds,
        'structures': fetch_structures,
        'docking_runs': fetch_docking_runs,
        'jobs': fetch_jobs,
    }
    
    if cache_name in cache_map:
        cache_map[cache_name].clear()


def with_refresh_button(cache_func):
    """Decorator to add refresh button that clears specific cache.
    
    Usage:
        @with_refresh_button
        def my_cached_function():
            return fetch_programs()
    """
    @wraps(cache_func)
    def wrapper(*args, **kwargs):
        col1, col2 = st.columns([10, 1])
        with col2:
            if st.button("🔄", key=f"refresh_{cache_func.__name__}", help="Refresh data"):
                cache_func.clear()
                st.rerun()
        return cache_func(*args, **kwargs)
    return wrapper


def cache_metrics() -> dict:
    """Get cache metrics for monitoring.
    
    Returns:
        Dictionary with cache statistics
    """
    # Note: Streamlit doesn't expose cache hit/miss stats directly
    # This is a placeholder for future implementation
    return {
        "cache_functions": 7,
        "ttl_seconds": DEFAULT_TTL,
        "status": "operational"
    }
