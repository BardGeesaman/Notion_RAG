"""Standardized loading state components for dashboard."""

import streamlit as st
from typing import Callable, Any, Optional
from functools import wraps


def loading_spinner(message: str = "Loading..."):
    """Decorator for consistent loading states on functions.
    
    Args:
        message: Loading message to display (e.g., "Loading datasets...")
    
    Usage:
        @loading_spinner("Loading experiments...")
        def fetch_experiments():
            return api_call()
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            with st.spinner(message):
                return func(*args, **kwargs)
        return wrapper
    return decorator


def render_loading_placeholder(
    message: str = "Loading data...",
    height: int = 200
):
    """Render a loading placeholder with consistent styling.
    
    Args:
        message: Loading message to display
        height: Height of the placeholder in pixels
    """
    st.markdown(
        f"""
        <div style="
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            height: {height}px;
            background: linear-gradient(90deg, #f0f0f0 25%, #e0e0e0 50%, #f0f0f0 75%);
            background-size: 200% 100%;
            animation: shimmer 1.5s infinite;
            border-radius: 8px;
            margin: 10px 0;
        ">
            <p style="color: #666; margin: 0;">⏳ {message}</p>
        </div>
        <style>
            @keyframes shimmer {{
                0% {{ background-position: -200% 0; }}
                100% {{ background-position: 200% 0; }}
            }}
        </style>
        """,
        unsafe_allow_html=True
    )


def render_skeleton_table(rows: int = 5, cols: int = 4):
    """Render a skeleton loader for data tables.
    
    Args:
        rows: Number of skeleton rows
        cols: Number of skeleton columns
    """
    placeholder_data = [["▓" * 10 for _ in range(cols)] for _ in range(rows)]
    st.dataframe(
        placeholder_data,
        use_container_width=True,
        hide_index=True
    )


def render_skeleton_metrics(count: int = 4):
    """Render skeleton loaders for metric cards.
    
    Args:
        count: Number of metric placeholders to show
    """
    cols = st.columns(count)
    for col in cols:
        with col:
            st.metric(label="Loading...", value="—", delta=None)


def render_skeleton_chart(height: int = 300):
    """Render a skeleton loader for charts.
    
    Args:
        height: Height of the chart placeholder
    """
    st.markdown(
        f"""
        <div style="
            height: {height}px;
            background: linear-gradient(45deg, #f8f9fa 25%, transparent 25%), 
                        linear-gradient(-45deg, #f8f9fa 25%, transparent 25%), 
                        linear-gradient(45deg, transparent 75%, #f8f9fa 75%), 
                        linear-gradient(-45deg, transparent 75%, #f8f9fa 75%);
            background-size: 20px 20px;
            background-position: 0 0, 0 10px, 10px -10px, -10px 0px;
            border-radius: 8px;
            border: 1px solid #e9ecef;
            display: flex;
            align-items: center;
            justify-content: center;
            margin: 10px 0;
        ">
            <p style="color: #6c757d; margin: 0;">📊 Loading chart...</p>
        </div>
        """,
        unsafe_allow_html=True
    )


class LoadingContext:
    """Context manager for loading states with progress tracking.
    
    Usage:
        with LoadingContext("Processing files...", show_progress=True) as loader:
            for i, file in enumerate(files):
                process_file(file)
                loader.update_progress(i / len(files), f"Processing {file}...")
    """
    
    def __init__(self, message: str = "Processing...", show_progress: bool = False):
        self.message = message
        self.show_progress = show_progress
        self._progress_bar = None
        self._spinner_container = None
        
    def __enter__(self):
        if self.show_progress:
            self._progress_bar = st.progress(0, text=self.message)
        else:
            self._spinner_container = st.spinner(self.message)
            self._spinner_container.__enter__()
        return self
    
    def update_progress(self, value: float, text: Optional[str] = None):
        """Update progress bar (0.0 to 1.0).
        
        Args:
            value: Progress value between 0.0 and 1.0
            text: Optional text to display with progress
        """
        if self._progress_bar:
            # Ensure value is within bounds
            value = max(0.0, min(1.0, value))
            self._progress_bar.progress(value, text=text or self.message)
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._progress_bar:
            self._progress_bar.empty()
        elif self._spinner_container:
            self._spinner_container.__exit__(exc_type, exc_val, exc_tb)
        return False


def with_loading_state(
    loading_message: str = "Loading...",
    empty_message: str = "No data found",
    error_message: str = "Failed to load data"
):
    """Decorator that adds loading, empty, and error states to functions.
    
    Args:
        loading_message: Message to show while loading
        empty_message: Message to show when result is empty
        error_message: Message to show on error
    
    Usage:
        @with_loading_state("Loading experiments...", "No experiments found")
        def render_experiments():
            return fetch_experiments()
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            try:
                with st.spinner(loading_message):
                    result = func(*args, **kwargs)
                
                # Handle empty results
                if not result:
                    st.info(f"ℹ️ {empty_message}")
                    return None
                
                return result
                
            except Exception as e:
                st.error(f"❌ {error_message}: {str(e)}")
                return None
        return wrapper
    return decorator


# Standardized loading messages for common operations
LOADING_MESSAGES = {
    "datasets": "Loading datasets...",
    "experiments": "Loading experiments...", 
    "programs": "Loading programs...",
    "compounds": "Loading compounds...",
    "structures": "Loading protein structures...",
    "jobs": "Loading jobs...",
    "results": "Loading results...",
    "analysis": "Running analysis...",
    "processing": "Processing data...",
    "saving": "Saving changes...",
    "uploading": "Uploading file...",
    "downloading": "Preparing download...",
    "generating": "Generating report...",
    "validating": "Validating data...",
    "importing": "Importing data...",
    "exporting": "Exporting data...",
    "refreshing": "Refreshing cache...",
    "connecting": "Connecting to service...",
}


def get_loading_message(operation: str, entity: str = None) -> str:
    """Get standardized loading message for common operations.
    
    Args:
        operation: Type of operation (e.g., "loading", "processing", "saving")
        entity: Entity being operated on (e.g., "datasets", "experiments")
    
    Returns:
        Standardized loading message
    
    Examples:
        get_loading_message("loading", "datasets") -> "Loading datasets..."
        get_loading_message("processing") -> "Processing data..."
    """
    if entity and f"{operation}_{entity}" in LOADING_MESSAGES:
        return LOADING_MESSAGES[f"{operation}_{entity}"]
    elif entity:
        return f"{operation.title()} {entity}..."
    elif operation in LOADING_MESSAGES:
        return LOADING_MESSAGES[operation]
    else:
        return f"{operation.title()}..."
