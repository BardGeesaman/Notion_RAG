import os

import requests
import streamlit as st


def show_login_page():
    st.title("Login")
    username = st.text_input("Username")
    password = st.text_input("Password", type="password")
    if st.button("Log in"):
        # For demo, local check; normally should POST to /api/v1/token
        resp = requests.post(
            f"{os.getenv('API_URL', 'http://localhost:8000')}/api/v1/token",
            data={"username": username, "password": password},
        )
        if resp.status_code == 200:
            data = resp.json()
            st.session_state["auth_user"] = {"username": username}
            st.session_state["access_token"] = data["access_token"]
            st.success("Login successful!")
            st.rerun()
        else:
            st.warning(f"Login failed: {resp.json().get('detail', resp.status_code)}")


def require_auth():
    if not os.getenv("AUTH_ENABLED", "0").isdigit() or int(os.getenv("AUTH_ENABLED", "0")) == 0:
        # When auth is disabled, return admin user for local development
        return {"username": "dev_user", "role": "admin"}
    user = st.session_state.get("auth_user")
    if not user:
        show_login_page()
        st.stop()
    return user


def require_admin(user=None):
    user = user or require_auth()
    if user.get("role") != "admin":
        st.error("Admin access required.")
        st.stop()
    return user
