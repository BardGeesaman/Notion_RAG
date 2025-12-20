"""Favorites and bookmarks helpers."""

from __future__ import annotations

from uuid import UUID as UUID_cls

import streamlit as st

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import UserFavorite, Experiment, Compound, Signature
from amprenta_rag.utils.bookmarks import get_user_bookmarks


def get_user_favorites(user_id: str) -> list[str]:
    if not user_id or user_id == "00000000-0000-0000-0000-000000000001":
        return []
    try:
        with db_session() as db:
            favorites = (
                db.query(UserFavorite)
                .filter(UserFavorite.user_id == UUID_cls(user_id))
                .order_by(UserFavorite.created_at.desc())
                .all()
            )
            return [f.page_name for f in favorites if f.page_name]
    except Exception:
        return []


def toggle_favorite(user_id: str, page_name: str) -> None:
    if not user_id or user_id == "00000000-0000-0000-0000-000000000001":
        return
    try:
        with db_session() as db:
            existing = (
                db.query(UserFavorite)
                .filter(UserFavorite.user_id == UUID_cls(user_id), UserFavorite.page_name == page_name)
                .first()
            )
            if existing:
                db.delete(existing)
            else:
                favorite = UserFavorite(user_id=UUID_cls(user_id), page_name=page_name)  # type: ignore[arg-type]
                db.add(favorite)
            db.commit()
    except Exception:
        pass


def update_recent_pages(page_name: str) -> None:
    if "recent_pages" not in st.session_state:
        st.session_state["recent_pages"] = []
    if page_name not in st.session_state["recent_pages"]:
        st.session_state["recent_pages"].insert(0, page_name)
        st.session_state["recent_pages"] = st.session_state["recent_pages"][:5]
    else:
        st.session_state["recent_pages"].remove(page_name)
        st.session_state["recent_pages"].insert(0, page_name)


def render_favorites_section(user: dict | None, update_recent_pages_fn) -> None:
    if not user or user.get("id") == "00000000-0000-0000-0000-000000000001":
        return
    favorites = get_user_favorites(str(user.get("id")))
    if favorites:
        with st.sidebar.expander("⭐ Favorites", expanded=True):
            for fav_page in favorites:
                if st.button(f"📌 {fav_page}", key=f"fav_{fav_page}", use_container_width=True):
                    st.session_state["selected_page"] = fav_page
                    update_recent_pages_fn(fav_page)
                    st.rerun()
        st.sidebar.divider()


def render_bookmarks_section(user: dict | None) -> None:
    if not user or user.get("id") == "00000000-0000-0000-0000-000000000001":
        return
    try:
        with db_session() as db:
            bookmarks = get_user_bookmarks(str(user.get("id")), db)
            if not bookmarks:
                return
            with st.sidebar.expander("📌 Bookmarks", expanded=False):
                experiments = [b for b in bookmarks if b.entity_type == "experiment"]
                compounds = [b for b in bookmarks if b.entity_type == "compound"]
                signatures = [b for b in bookmarks if b.entity_type == "signature"]

                if experiments:
                    st.markdown("**Experiments**")
                    for bm in experiments:
                        exp = db.query(Experiment).filter(Experiment.id == bm.entity_id).first()
                        if exp and st.button(f"🔬 {exp.name}", key=f"bm_exp_{bm.id}", use_container_width=True):
                            st.session_state["selected_page"] = "Experiments"
                            st.session_state["selected_experiment_id"] = str(bm.entity_id)
                            st.rerun()

                if compounds:
                    st.markdown("**Compounds**")
                    for bm in compounds:
                        comp = db.query(Compound).filter(Compound.id == bm.entity_id).first()
                        if comp and st.button(f"⚗️ {comp.compound_id}", key=f"bm_comp_{bm.id}", use_container_width=True):
                            st.session_state["selected_page"] = "Chemistry"
                            st.session_state["selected_compound_id"] = str(bm.entity_id)
                            st.rerun()

                if signatures:
                    st.markdown("**Signatures**")
                    for bm in signatures:
                        sig = db.query(Signature).filter(Signature.id == bm.entity_id).first()
                        if sig and st.button(f"📊 {sig.name}", key=f"bm_sig_{bm.id}", use_container_width=True):
                            st.session_state["selected_page"] = "Signatures"
                            st.session_state["selected_signature_id"] = str(bm.entity_id)
                            st.rerun()
    except Exception:
        pass

