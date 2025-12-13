"""Login page for Amprenta dashboard."""
import streamlit as st
from amprenta_rag.database.models import User
from amprenta_rag.database.session import db_session
from amprenta_rag.auth.password import verify_password
from amprenta_rag.auth.session import set_current_user
from amprenta_rag.auth.audit import log_login
from datetime import datetime


def render_login_page():
    st.title("🔐 Login")
    st.markdown("Sign in to access the Amprenta Multi-Omics Platform")

    with st.form("login_form"):
        username = st.text_input("Username")
        password = st.text_input("Password", type="password")
        submitted = st.form_submit_button("Login", type="primary")

        if submitted:
            if not username or not password:
                st.error("Please enter username and password")
                return

            with db_session() as db:
                user = db.query(User).filter(User.username == username).first()
                if user and user.is_active and verify_password(password, user.password_hash):
                    set_current_user({
                        "id": str(user.id),
                        "username": user.username,
                        "email": user.email,
                        "role": user.role
                    })
                    user.last_login = datetime.utcnow()
                    db.commit()
                    log_login(str(user.id), user.username)
                    st.success("Login successful!")
                    st.rerun()
                else:
                    st.error("Invalid username or password")


if __name__ == "__main__":
    render_login_page()
