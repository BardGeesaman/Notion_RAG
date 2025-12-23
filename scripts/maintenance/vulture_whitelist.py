"""
Vulture whitelist for known false positives.

Vulture works by flagging "unused" code. In modern Python apps, many call sites are
dynamic (framework callbacks, decorators, reflection), so we list a set of common
entrypoints and symbols here to reduce noise.

This file should be included in vulture runs:
  vulture ... scripts/maintenance/vulture_whitelist.py
"""

# __init__/dunder methods
__all__ = []


def __str__():  # type: ignore[no-redef]
    ...


def __repr__():  # type: ignore[no-redef]
    ...


# Alembic migrations
def upgrade():
    ...


def downgrade():
    ...


# FastAPI dependency patterns
def get_db():
    ...


def get_current_user():
    ...


def get_current_user_with_company():
    ...


def set_company_context():
    ...


def extract_subdomain():
    ...


# Pydantic validators / model hooks (common names)
def validator():
    ...


def root_validator():
    ...


def field_validator():
    ...


def model_validator():
    ...


# Pytest fixtures / hooks (common names)
def pytest_addoption():
    ...


def pytest_configure():
    ...


def client():
    ...


def db_session():
    ...


def page():
    ...


def base_url():
    ...


def streamlit_server():
    ...


