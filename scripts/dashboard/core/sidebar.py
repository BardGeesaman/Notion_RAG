"""Sidebar rendering and navigation."""

from __future__ import annotations

import streamlit as st
from typing import Iterable

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
from scripts.dashboard.core.jupyter_auth import get_jupyterhub_url


def render_user_info(user: dict | None):
    if not user:
        return
    st.markdown(f"**👤 {user['username']}** ({user['role']})")
    if st.button("Logout", key="logout_btn"):
        log_logout(user.get("id"), user.get("username"))
        from amprenta_rag.auth.session import clear_session

        clear_session()
        st.rerun()
    st.divider()


def render_notifications(user: dict | None):
    if not user:
        return
    user_id = user.get("id")
    if user_id in (None, "00000000-0000-0000-0000-000000000001"):
        return
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


def render_theme_selector():
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


def render_search():
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


def render_navigation(user: dict | None, visible_pages: Iterable[str], *, groups) -> str:
    page = st.session_state.get("selected_page", "Overview")

    filtered_discovery_pages = [p for p in groups["DISCOVERY_PAGES"] if p in visible_pages]
    filtered_analysis_pages = [p for p in groups["ANALYSIS_PAGES"] if p in visible_pages]
    filtered_eln_pages = [p for p in groups["ELN_PAGES"] if p in visible_pages]
    filtered_admin_pages = [p for p in groups["ADMIN_PAGES"] if p in visible_pages]

    if filtered_discovery_pages:
        with st.expander("🔍 Discovery", expanded=True):
            for p in filtered_discovery_pages:
                if st.button(p, key=f"btn_{p}", use_container_width=True):
                    page = p
                    st.session_state["selected_page"] = p
                    update_recent_pages(p)
                    st.rerun()

    if filtered_analysis_pages:
        with st.expander("📊 Analysis", expanded=False):
            for p in filtered_analysis_pages:
                if st.button(p, key=f"btn_{p}", use_container_width=True):
                    page = p
                    st.session_state["selected_page"] = p
                    update_recent_pages(p)
                    st.rerun()

    if filtered_eln_pages:
        with st.expander("📋 ELN", expanded=False):
            for p in filtered_eln_pages:
                if st.button(p, key=f"btn_{p}", use_container_width=True):
                    page = p
                    st.session_state["selected_page"] = p
                    update_recent_pages(p)
                    st.rerun()

    if user and user.get("role") == "admin" and filtered_admin_pages:
        with st.expander("⚙️ Admin", expanded=False):
            for p in filtered_admin_pages:
                if st.button(p, key=f"btn_{p}", use_container_width=True):
                    page = p
                    st.session_state["selected_page"] = p
                    update_recent_pages(p)
                    st.rerun()

    other_pages = [p for p in groups["ALL_PAGES"] if p not in groups["DISCOVERY_PAGES"] + groups["ANALYSIS_PAGES"] + groups["ELN_PAGES"] + groups["ADMIN_PAGES"]]
    filtered_other_pages = [p for p in other_pages if p in visible_pages]
    if filtered_other_pages:
        with st.expander("📚 Other Pages", expanded=False):
            for p in filtered_other_pages:
                if st.button(p, key=f"btn_{p}", use_container_width=True):
                    page = p
                    st.session_state["selected_page"] = p
                    update_recent_pages(p)
                    st.rerun()

    return page


def render_recent():
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


def render_sidebar(user: dict | None, visible_pages: Iterable[str], groups) -> str:
    with st.sidebar:
        st.title("Navigation")

        render_user_info(user)
        render_notifications(user)
        render_theme_selector()

        if user and user.get("role") == "admin":
            if st.button("➕ Register User", key="register_btn"):
                from scripts.dashboard.pages.auth.register import render_register_page

                st.session_state["show_register"] = True
                st.rerun()

        render_search()
        st.divider()

        render_favorites_section(user, update_recent_pages)
        render_bookmarks_section(user)
        render_recent()

        if "selected_page" not in st.session_state:
            st.session_state["selected_page"] = "Overview"

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

        page = render_navigation(
            user,
            visible_pages,
            groups={
                "DISCOVERY_PAGES": groups.DISCOVERY_PAGES,
                "ANALYSIS_PAGES": groups.ANALYSIS_PAGES,
                "ELN_PAGES": groups.ELN_PAGES,
                "ADMIN_PAGES": groups.ADMIN_PAGES,
                "ALL_PAGES": groups.ALL_PAGES,
            },
        )

        page = st.session_state.get("selected_page", page)

        if user and user.get("id") != "00000000-0000-0000-0000-000000000001":
            favorites = get_user_favorites(user.get("id"))
            is_fav = page in favorites
            star_label = "⭐" if is_fav else "☆"
            if st.button(f"{star_label} Favorite", key="toggle_favorite", use_container_width=True):
                toggle_favorite(user.get("id"), page)
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
        
        # Extract base URL for Voila
        jupyter_base_url = jupyter_url.split("/hub/login")[0] if "/hub/login" in jupyter_url else "http://localhost:8888"
        voila_url = f"{jupyter_base_url}/user/{username}/voila/render/templates/experiment_dashboard.ipynb"
        st.sidebar.link_button(
            "View as Dashboard",
            voila_url,
            use_container_width=True
        )
        st.sidebar.caption("Login to JupyterHub first")

        return page

