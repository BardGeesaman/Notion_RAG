"""User registration page (admin only)."""
import streamlit as st
from amprenta_rag.database.models import User
from amprenta_rag.database.session import db_session
from amprenta_rag.auth.password import hash_password
from amprenta_rag.auth.session import get_current_user


def render_register_page():
    user = get_current_user()

    # Only admins can register new users
    if not user or user.get("role") != "admin":
        st.error("Access denied. Only administrators can register new users.")
        return

    st.title("👥 Register New User")
    st.markdown("Create a new user account")

    with st.form("register_form"):
        username = st.text_input("Username", max_chars=100)
        email = st.text_input("Email")
        password = st.text_input("Password", type="password")
        confirm_password = st.text_input("Confirm Password", type="password")
        role = st.selectbox("Role", ["researcher", "viewer", "admin"])
        submitted = st.form_submit_button("Create User", type="primary")

        if submitted:
            if not username or not email or not password:
                st.error("All fields are required")
                return
            if password != confirm_password:
                st.error("Passwords do not match")
                return
            if len(password) < 8:
                st.error("Password must be at least 8 characters")
                return

            with db_session() as db:
                # Check for existing user
                existing = db.query(User).filter(
                    (User.username == username) | (User.email == email)
                ).first()
                if existing:
                    st.error("Username or email already exists")
                    return

                new_user = User(
                    username=username,
                    email=email,
                    password_hash=hash_password(password),
                    role=role,
                    is_active=True
                )
                db.add(new_user)
                db.commit()
                st.success(f"User '{username}' created successfully!")


if __name__ == "__main__":
    render_register_page()
