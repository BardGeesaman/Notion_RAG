
"""Login page for Amprenta dashboard."""
import streamlit as st
from amprenta_rag.database.models import User
from amprenta_rag.database.session import db_session
from amprenta_rag.auth.password import verify_password
from amprenta_rag.auth.session import set_current_user
from amprenta_rag.auth.audit import log_login
from datetime import datetime

from scripts.dashboard.utils.accessibility import (
    render_skip_link,
    accessible_text_input,
    accessible_button,
    announce_to_screen_reader,
    add_heading_structure,
    ensure_minimum_contrast
)


def render_login_page():
    # Add accessibility features
    render_skip_link("login-form")
    ensure_minimum_contrast()
    
    # Structured heading for screen readers
    add_heading_structure("🔐 Login", level=1, id="login-heading")
    st.markdown("Sign in to access the Amprenta Multi-Omics Platform")
    
    # Add form role and ARIA labels
    st.markdown(
        """
        <div role="form" aria-labelledby="login-heading" aria-describedby="login-description">
            <div id="login-description" style="margin-bottom: 1rem; color: #6c757d;">
                Enter your username and password to access the platform
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )

    with st.form("login_form", clear_on_submit=False):
        # Accessible form inputs with proper ARIA labels
        username = accessible_text_input(
            "Username",
            key="username_input",
            aria_label="Username for login",
            aria_describedby="username-help",
            help="Enter your registered username"
        )
        
        # Add hidden help text for screen readers
        st.markdown(
            '<div id="username-help" class="sr-only">Enter the username you registered with</div>',
            unsafe_allow_html=True
        )
        
        password = accessible_text_input(
            "Password", 
            key="password_input",
            type="password",
            aria_label="Password for login",
            aria_describedby="password-help",
            help="Enter your password"
        )
        
        # Add hidden help text for screen readers
        st.markdown(
            '<div id="password-help" class="sr-only">Enter your account password</div>',
            unsafe_allow_html=True
        )
        
        submitted = accessible_button(
            "Login", 
            key="login_submit",
            aria_label="Submit login form",
            help="Click to sign in to the platform"
        )

        if submitted:
            # Announce form submission to screen readers
            announce_to_screen_reader("Processing login request", priority="assertive")
            
            if not username or not password:
                error_msg = "Please enter both username and password"
                st.error(f"❌ {error_msg}")
                # Add error announcement with proper ARIA
                st.markdown(
                    f"""
                    <div role="alert" aria-live="assertive" id="login-error">
                        {error_msg}
                    </div>
                    """,
                    unsafe_allow_html=True
                )
                announce_to_screen_reader(error_msg, priority="assertive")
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
                    
                    success_msg = f"Login successful! Welcome back, {user.username}"
                    st.success(f"✅ {success_msg}")
                    
                    # Add success announcement with proper ARIA
                    st.markdown(
                        f"""
                        <div role="status" aria-live="polite" id="login-success">
                            {success_msg}
                        </div>
                        """,
                        unsafe_allow_html=True
                    )
                    announce_to_screen_reader(success_msg, priority="assertive")
                    st.rerun()
                else:
                    error_msg = "Invalid username or password. Please check your credentials and try again."
                    st.error(f"❌ {error_msg}")
                    
                    # Add error announcement with proper ARIA
                    st.markdown(
                        f"""
                        <div role="alert" aria-live="assertive" id="login-error">
                            {error_msg}
                        </div>
                        """,
                        unsafe_allow_html=True
                    )
                    announce_to_screen_reader(error_msg, priority="assertive")


if __name__ == "__main__":
    render_login_page()
