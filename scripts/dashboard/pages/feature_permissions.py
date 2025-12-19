"""Feature Permissions management page (admin only)."""
from __future__ import annotations

import streamlit as st

from amprenta_rag.auth.session import get_current_user
from amprenta_rag.auth.feature_permissions import (
    get_all_permissions,
    set_page_visibility,
    initialize_default_permissions,
)
from scripts.dashboard.db_session import db_session


def render_feature_permissions_page() -> None:
    """Render the Feature Permissions management page."""
    st.header("🔐 Feature Permissions")
    st.markdown("Control which pages are visible to each user role")

    user = get_current_user()

    # Admin only
    if not user or user.get("role") != "admin":
        st.error("Access denied. Only administrators can manage feature permissions.")
        return

    with db_session() as db:
        # Get current permissions
        permissions = get_all_permissions(db)

        # Get all pages
        from scripts.run_dashboard import ALL_PAGES

        st.markdown("---")

        # Reset to defaults button
        if st.button("🔄 Reset to Defaults", type="secondary", help="Reset all permissions to default values"):
            try:
                # Delete all existing permissions
                from amprenta_rag.database.models import FeaturePermission
                db.query(FeaturePermission).delete()
                db.commit()

                # Initialize defaults
                initialize_default_permissions(db)
                st.success("Permissions reset to defaults!")
                st.rerun()
            except Exception as e:
                st.error(f"Error resetting permissions: {e}")

        st.markdown("---")

        # Display permissions matrix for each role
        roles = ["admin", "researcher", "viewer"]

        for role in roles:
            st.subheader(f"Role: {role.title()}")

            # Get current permissions for this role
            role_perms = permissions.get(role, {})

            # Create columns for checkboxes (4 columns)
            num_cols = 4
            cols = st.columns(num_cols)

            # Group pages into columns
            for i, page_name in enumerate(ALL_PAGES):
                col_idx = i % num_cols
                with cols[col_idx]:
                    # Get current visibility (default True if not set)
                    current_visible = role_perms.get(page_name, True)

                    # Create checkbox
                    new_visible = st.checkbox(
                        page_name,
                        value=current_visible,
                        key=f"perm_{role}_{page_name}",
                        help=f"Toggle visibility of '{page_name}' for {role} role",
                    )

                    # Update if changed
                    if new_visible != current_visible:
                        try:
                            success = set_page_visibility(role, page_name, new_visible, db)
                            if success:
                                st.success(f"Updated {page_name} for {role}")
                                st.rerun()
                            else:
                                st.error(f"Failed to update {page_name} for {role}")
                        except Exception as e:
                            st.error(f"Error updating permission: {e}")

            st.markdown("---")

        # Summary statistics
        st.subheader("Summary")
        col1, col2, col3 = st.columns(3)

        for i, role in enumerate(roles):
            role_perms = permissions.get(role, {})
            visible_count = sum(1 for v in role_perms.values() if v)
            total_count = len(role_perms) if role_perms else len(ALL_PAGES)

            with [col1, col2, col3][i]:
                st.metric(
                    f"{role.title()}",
                    f"{visible_count}/{total_count} visible",
                )

