"""Email Settings page for managing email subscriptions."""
from __future__ import annotations

import streamlit as st
from uuid import UUID

from amprenta_rag.database.models import EmailSubscription, User
from amprenta_rag.auth.session import get_current_user
from amprenta_rag.notifications.email_service import is_email_configured
from scripts.dashboard.db_session import db_session


def render_email_settings_page() -> None:
    """Render the Email Settings page."""
    st.header("📧 Email Settings")
    st.markdown("Manage your email notification preferences")
    
    user = get_current_user()
    if not user or user.get("id") == "test":
        st.warning("Email settings require a logged-in user account.")
        st.info("Currently in test/dev mode. Log in with a real account to manage email preferences.")
        return
    
    user_id = user.get("id")
    if not user_id:
        st.error("Invalid user session.")
        return
    
    # Check if email is configured
    if not is_email_configured():
        st.warning("⚠️ Email service is not configured. Contact your administrator to enable email notifications.")
        return
    
    with db_session() as db:
        # Get user email
        user_obj = db.query(User).filter(User.id == UUID(user_id)).first()
        if not user_obj:
            st.error("User not found.")
            return
        
        st.info(f"📬 Notifications will be sent to: **{user_obj.email}**")
        st.markdown("---")
        
        # Get existing subscriptions
        subscriptions = db.query(EmailSubscription).filter(
            EmailSubscription.user_id == UUID(user_id)
        ).all()
        
        # Create subscription dict for easy lookup
        subscription_dict = {sub.subscription_type: sub for sub in subscriptions}
        
        # Digest subscription
        st.subheader("📊 Activity Digest")
        digest_sub = subscription_dict.get("digest")
        digest_frequency = digest_sub.frequency if digest_sub and digest_sub.is_active else "off"
        
        digest_options = ["off", "daily", "weekly"]
        selected_digest = st.radio(
            "Receive activity digest",
            digest_options,
            index=digest_options.index(digest_frequency) if digest_frequency in digest_options else 0,
            key="digest_frequency",
            help="Daily or weekly summary of platform activity"
        )
        
        if selected_digest != digest_frequency:
            if selected_digest == "off":
                # Deactivate or delete subscription
                if digest_sub:
                    digest_sub.is_active = False
                    db.commit()
                    st.success("Digest subscription disabled")
                    st.rerun()
            else:
                # Create or update subscription
                if digest_sub:
                    digest_sub.frequency = selected_digest
                    digest_sub.is_active = True
                else:
                    digest_sub = EmailSubscription(
                        user_id=UUID(user_id),
                        subscription_type="digest",
                        frequency=selected_digest,
                        is_active=True,
                    )
                    db.add(digest_sub)
                db.commit()
                st.success(f"Digest subscription set to {selected_digest}")
                st.rerun()
        
        st.markdown("---")
        
        # Alert notifications
        st.subheader("🔔 Alert Notifications")
        alert_sub = subscription_dict.get("alerts")
        alerts_enabled = alert_sub.is_active if alert_sub else False
        
        alerts_toggle = st.checkbox(
            "Enable immediate alert notifications",
            value=alerts_enabled,
            key="alerts_toggle",
            help="Receive immediate email notifications for important events"
        )
        
        if alerts_toggle != alerts_enabled:
            if alerts_toggle:
                # Enable alerts
                if alert_sub:
                    alert_sub.is_active = True
                    alert_sub.frequency = "immediate"
                else:
                    alert_sub = EmailSubscription(
                        user_id=UUID(user_id),
                        subscription_type="alerts",
                        frequency="immediate",
                        is_active=True,
                    )
                    db.add(alert_sub)
                db.commit()
                st.success("Alert notifications enabled")
                st.rerun()
            else:
                # Disable alerts
                if alert_sub:
                    alert_sub.is_active = False
                    db.commit()
                    st.success("Alert notifications disabled")
                    st.rerun()
        
        st.markdown("---")
        
        # Share notifications
        st.subheader("🔗 Share Notifications")
        share_sub = subscription_dict.get("shares")
        shares_enabled = share_sub.is_active if share_sub else True  # Default to enabled
        
        shares_toggle = st.checkbox(
            "Enable share notifications",
            value=shares_enabled,
            key="shares_toggle",
            help="Receive email notifications when someone shares content with you"
        )
        
        if shares_toggle != shares_enabled:
            if shares_toggle:
                # Enable shares
                if share_sub:
                    share_sub.is_active = True
                    share_sub.frequency = "immediate"
                else:
                    share_sub = EmailSubscription(
                        user_id=UUID(user_id),
                        subscription_type="shares",
                        frequency="immediate",
                        is_active=True,
                    )
                    db.add(share_sub)
                db.commit()
                st.success("Share notifications enabled")
                st.rerun()
            else:
                # Disable shares
                if share_sub:
                    share_sub.is_active = False
                    db.commit()
                    st.success("Share notifications disabled")
                    st.rerun()
        
        st.markdown("---")
        
        # Show current subscriptions summary
        st.subheader("Current Subscriptions")
        active_subs = db.query(EmailSubscription).filter(
            EmailSubscription.user_id == UUID(user_id),
            EmailSubscription.is_active == True
        ).all()
        
        if active_subs:
            for sub in active_subs:
                st.write(f"✅ **{sub.subscription_type.title()}**: {sub.frequency}")
        else:
            st.info("No active email subscriptions")

