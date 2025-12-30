"""Sidebar rendering and navigation."""

from __future__ import annotations

import os
import streamlit as st
from typing import Iterable, Any, Dict, List, Optional

from amprenta_rag.database.session import db_session
from amprenta_rag.auth.audit import log_logout
from amprenta_rag.notifications.service import get_unread_count, get_unread_notifications, mark_as_read, mark_all_read
from amprenta_rag.utils.global_search import global_search
from scripts.dashboard.themes import apply_theme, get_theme_names
from scripts.dashboard.shortcuts import inject_shortcuts_js, handle_shortcut_action, render_shortcuts_help
from scripts.dashboard.help_content import get_help, search_help
from scripts.dashboard.help_chat import render_help_chat
from scripts.dashboard.components.command_palette import render_command_palette
from scripts.dashboard.core.favorites import (
    render_favorites_section,
    render_bookmarks_section,
    toggle_favorite,
    get_user_favorites,
    update_recent_pages,
)
# Import moved inside function to avoid circular import
from scripts.dashboard.core.jupyter_auth import get_jupyterhub_url, get_voila_url
from scripts.dashboard.components.alerts_bell import render_alerts_bell


def render_user_info(user: Dict[str, Any] | None) -> None:
    if not user:
        return
    username = str(user.get("username", ""))
    role = str(user.get("role", ""))
    st.markdown(f"**👤 {username}** ({role})")
    if st.button("Logout", key="logout_btn"):
        log_logout(str(user.get("id") or ""), username)
        from amprenta_rag.auth.session import clear_session

        clear_session()
        st.rerun()
    st.divider()


def render_notifications(user: Dict[str, Any] | None) -> None:
    if not user:
        return
    user_id_raw = user.get("id")
    if user_id_raw in (None, "00000000-0000-0000-0000-000000000001"):
        return
    user_id = str(user_id_raw)
    try:
        with db_session() as db:
            unread_count = get_unread_count(user_id, db)
            notifications = get_unread_notifications(user_id, db, limit=10)

            if unread_count > 0:
                st.markdown(f"**🔔 Notifications ({unread_count})**")
            else:
                st.markdown("**🔔 Notifications**")

            with st.expander("View Notifications", expanded=False):
                if notifications:
                    for notif in notifications:
                        col1, col2 = st.columns([3, 1])
                        with col1:
                            st.markdown(f"**{notif.title}**")
                            if notif.message:
                                st.caption(notif.message)
                            ts = notif.created_at.strftime("%Y-%m-%d %H:%M") if notif.created_at else ""
                            st.caption(f"{ts} • {notif.notification_type}")
                        with col2:
                            if st.button("✓", key=f"mark_read_{notif.id}", help="Mark as read"):
                                mark_as_read(str(notif.id), db)
                                st.rerun()
                        st.divider()
                    if st.button("Mark All Read", key="mark_all_read", use_container_width=True):
                        mark_all_read(user_id, db)
                        st.rerun()
                else:
                    st.info("No unread notifications")
    except Exception as e:
        st.error(f"Error loading notifications: {e}")
    st.divider()


def render_theme_selector() -> None:
    if "theme" not in st.session_state:
        st.session_state["theme"] = "dark"
    apply_theme(st.session_state["theme"])
    theme_options = get_theme_names()
    selected_theme = st.selectbox(
        "Theme",
        theme_options,
        index=theme_options.index(st.session_state["theme"]) if st.session_state["theme"] in theme_options else 0,
        key="theme_selector",
    )
    if selected_theme != st.session_state["theme"]:
        st.session_state["theme"] = selected_theme
        st.rerun()
    st.divider()


def render_search() -> None:
    search_query = st.text_input("🔍 Search...", key="global_search", placeholder="Search experiments, compounds, datasets...")
    if not search_query:
        return
    with st.expander("🔍 Search Results", expanded=True):
        try:
            with db_session() as db:
                results = global_search(search_query, db, limit=5)
                if results["experiments"]:
                    st.markdown("**Experiments**")
                    for exp in results["experiments"]:
                        if st.button(f"🔬 {exp['name']}", key=f"search_exp_{exp['id']}", use_container_width=True):
                            st.session_state["selected_page"] = "Experiments"
                            st.session_state["selected_experiment_id"] = exp["id"]
                            st.rerun()
                if results["compounds"]:
                    st.markdown("**Compounds**")
                    for comp in results["compounds"]:
                        if st.button(f"⚗️ {comp['compound_id']}", key=f"search_comp_{comp['id']}", use_container_width=True):
                            st.session_state["selected_page"] = "Chemistry"
                            st.session_state["selected_compound_id"] = comp["id"]
                            st.rerun()
                if results["signatures"]:
                    st.markdown("**Signatures**")
                    for sig in results["signatures"]:
                        if st.button(f"📊 {sig['name']}", key=f"search_sig_{sig['id']}", use_container_width=True):
                            st.session_state["selected_page"] = "Signatures"
                            st.session_state["selected_signature_id"] = sig["id"]
                            st.rerun()
                if results["datasets"]:
                    st.markdown("**Datasets**")
                    for ds in results["datasets"]:
                        if st.button(f"📁 {ds['name']}", key=f"search_ds_{ds['id']}", use_container_width=True):
                            st.session_state["selected_page"] = "Datasets"
                            st.session_state["selected_dataset_id"] = ds["id"]
                            st.rerun()
                if not any([results["experiments"], results["compounds"], results["signatures"], results["datasets"]]):
                    st.info("No results found")
        except Exception as e:
            st.error(f"Search error: {e}")


# DEPRECATED: render_navigation function removed
# This legacy function has been replaced by render_grouped_sidebar
# which uses the new PAGE_GROUPS structure for better organization


def render_recent() -> None:
    if st.session_state.get("recent_pages"):
        with st.expander("🕐 Recent", expanded=False):
            for recent_page in st.session_state["recent_pages"][:5]:
                if st.button(f"🕐 {recent_page}", key=f"recent_{recent_page}", use_container_width=True):
                    st.session_state["selected_page"] = recent_page
                    update_recent_pages(recent_page)
                    st.rerun()
        st.divider()


def render_help(page: str):
    with st.expander("❓ Help", expanded=False):
        help_query = st.text_input("Search help", key="help_search")
        if help_query:
            results = search_help(help_query)
            for r in results[:3]:
                st.markdown(f"**{r['title']}**: {r['description']}")
        else:
            current_help = get_help(page)
            if current_help:
                st.markdown(f"### {current_help['title']}")
                st.write(current_help['description'])
                if current_help.get("tips"):
                    st.markdown("**Tips:**")
                    for tip in current_help["tips"]:
                        st.markdown(f"• {tip}")

    with st.expander("🤖 Ask Assistant", expanded=False):
        render_help_chat()


def render_sidebar(user: Dict[str, Any] | None, visible_pages: Iterable[str], groups: Any) -> str:
    with st.sidebar:
        st.title("Navigation")

        render_user_info(user)
        render_alerts_bell()
        render_notifications(user)
        render_theme_selector()

        if user and user.get("role") == "admin":
            if st.button("➕ Register User", key="register_btn"):

                st.session_state["show_register"] = True
                st.rerun()

        render_search()
        st.divider()

        render_favorites_section(user, update_recent_pages)
        render_bookmarks_section(user)
        render_recent()

        if "selected_page" not in st.session_state:
            st.session_state["selected_page"] = "Overview"

        # Allow deep-linking via URL query param, e.g. /?page=Cytoscape%20Demo
        # This is especially useful for smoke tests and bookmarking.
    page_param: Optional[str] = None

    # Debug: show raw query params if enabled
    if os.environ.get("AMPRENTA_DEBUG_NAV", "").lower() in ("1", "true", "yes"):
        try:
            st.sidebar.caption(f"raw_query_params: {dict(st.query_params)}")
        except Exception as e:  # noqa: BLE001
            st.sidebar.caption(f"query_params error: {e}")

    # Streamlit 1.30+ query param reading
    try:
        if hasattr(st, "query_params") and "page" in st.query_params:
            raw = st.query_params["page"]
            page_param = raw[0] if isinstance(raw, list) else (str(raw) if raw else None)
    except Exception:
        pass

    # Fallback for older Streamlit
    if page_param is None:
        try:
            raw2 = st.experimental_get_query_params().get("page")
            page_param = raw2[0] if isinstance(raw2, list) and raw2 else None
        except Exception:
            pass

    # Normalize query param (some setups return URL-encoded values).
    if page_param:
        try:
            from urllib.parse import unquote_plus

            page_param = unquote_plus(str(page_param)).strip()
        except Exception:
            page_param = str(page_param).strip()

    if os.environ.get("AMPRENTA_DEBUG_NAV", "").lower() in ("1", "true", "yes"):
        in_visible = bool(page_param and page_param in set(visible_pages))
        st.sidebar.caption(
            f"nav_debug: page_param={page_param!r} "
            f"selected={st.session_state.get('selected_page')!r} "
            f"in_visible={in_visible}"
        )
        if page_param and not in_visible:
            st.sidebar.warning(f"nav_debug: unknown page param {page_param!r}")

    if page_param and page_param in set(visible_pages):
        if st.session_state.get("selected_page") != page_param:
            st.session_state["selected_page"] = page_param
            update_recent_pages(page_param)
            st.rerun()

        inject_shortcuts_js()
        shortcut_actions = [
            "navigate_home",
            "navigate_experiments",
            "navigate_chemistry",
            "navigate_analysis",
            "navigate_search",
            "navigate_ingestion",
            "navigate_protocols",
            "navigate_qa",
        ]
        for action in shortcut_actions:
            if st.session_state.get(f"shortcut_{action}", False):
                handle_shortcut_action(action)
                st.session_state.pop(f"shortcut_{action}", None)
                break

        if st.session_state.get("show_shortcuts", False):
            render_shortcuts_help()
            if st.button("Close", key="close_shortcuts"):
                st.session_state.pop("show_shortcuts", None)
                st.rerun()

        render_command_palette()

    # Use enhanced grouped navigation
    current_page = st.session_state.get("selected_page", "Overview")
    
    # Import here to avoid circular import
    from scripts.dashboard.components.sidebar_nav import render_grouped_sidebar
    selected = render_grouped_sidebar(current_page=current_page)
    
    if selected:
        page = selected
        st.session_state["selected_page"] = selected
        update_recent_pages(selected)
        st.rerun()
    else:
        page = current_page

    page = st.session_state.get("selected_page", page)

    if user and user.get("id") != "00000000-0000-0000-0000-000000000001":
        user_id_str = str(user.get("id"))
        favorites = get_user_favorites(user_id_str)
        is_fav = page in favorites
        star_label = "⭐" if is_fav else "☆"
        if st.button(f"{star_label} Favorite", key="toggle_favorite", use_container_width=True):
            toggle_favorite(user_id_str, page)
            st.rerun()

    render_help(page)
    update_recent_pages(page)

    st.sidebar.divider()
    st.sidebar.subheader("Jupyter Notebooks")

    # Get username from user object or default
    username = user.get("username", "scientist") if user else "scientist"

    # Get current page context from session state
    experiment_id = st.session_state.get("current_experiment_id")
    dataset_id = st.session_state.get("current_dataset_id")
    compound_id = st.session_state.get("current_compound_id")

    jupyter_url = get_jupyterhub_url(
        username,
        experiment_id=experiment_id,
        dataset_id=dataset_id,
        compound_id=compound_id
    )
    st.sidebar.link_button(
        "Open in JupyterLab",
        jupyter_url,
        use_container_width=True
    )
    st.sidebar.caption("Launch notebooks with API access")

    if st.sidebar.button("Generate Report", use_container_width=True):
        st.info("Report generation coming soon")

    # Determine dashboard based on context
    current_page = st.session_state.get("selected_page", "")
    campaign_id = st.session_state.get("current_campaign_id")

    if "HTS" in current_page or "Screening" in current_page or campaign_id:
        dashboard = "hts_plate_viewer.ipynb"
        context = {"campaign_id": campaign_id} if campaign_id else {}
    elif experiment_id:
        dashboard = "experiment_dashboard.ipynb"
        context = {"experiment_id": str(experiment_id)} if experiment_id else {}
    else:
        dashboard = "experiment_dashboard.ipynb"
        context = {}

    # Extract base URL for Voila
    jupyter_base_url = jupyter_url.split("/hub/login")[0] if "/hub/login" in jupyter_url else "http://localhost:8888"
    voila_url = get_voila_url(username, dashboard, context if context else None, jupyter_base_url)
    st.sidebar.link_button(
        "View as Dashboard",
        voila_url,
        use_container_width=True
    )
    st.sidebar.caption("Login to JupyterHub first")

    return page

