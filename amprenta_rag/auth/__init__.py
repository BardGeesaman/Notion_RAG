"""Authentication utilities."""
from amprenta_rag.auth.password import hash_password, verify_password
from amprenta_rag.auth.session import get_current_user, set_current_user, clear_session, require_auth
from amprenta_rag.auth.audit import log_action, log_login, log_logout, log_create, log_update, log_delete

__all__ = [
    "hash_password",
    "verify_password",
    "get_current_user",
    "set_current_user",
    "clear_session",
    "require_auth",
    "log_action",
    "log_login",
    "log_logout",
    "log_create",
    "log_update",
    "log_delete",
]
