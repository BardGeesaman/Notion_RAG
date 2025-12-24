"""Shared utilities for seeding scripts."""

from __future__ import annotations

import argparse
from typing import Any, Iterable, Optional, Type, TypeVar

from sqlalchemy import inspect
from sqlalchemy.exc import OperationalError, ProgrammingError
from sqlalchemy.orm import Session

T = TypeVar("T")


def get_common_parser(description: str) -> argparse.ArgumentParser:
    """Create a reusable argparse parser for seeders."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--size", choices=["small", "medium", "large"], default="small")
    parser.add_argument("--reset", action="store_true")
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def with_progress(iterable: Iterable[T], desc: str) -> Iterable[T]:
    """Wrap an iterable with tqdm if available (graceful fallback if not installed)."""
    try:
        from tqdm import tqdm  # type: ignore

        return tqdm(iterable, desc=desc)
    except Exception:
        return iterable


def _table_exists(session: Session, tablename: str) -> bool:
    try:
        insp = inspect(session.get_bind())
        return bool(insp.has_table(tablename))
    except Exception:
        return False


def safe_delete(session: Session, model: Type[Any]) -> int:
    """Delete all rows in a table if it exists; otherwise no-op."""
    tablename = getattr(model, "__tablename__", None)
    if not tablename:
        return 0
    if not _table_exists(session, tablename):
        return 0
    try:
        return int(session.query(model).delete(synchronize_session=False))
    except (OperationalError, ProgrammingError):
        # Table may not exist in older schemas
        return 0


def validate_counts(session: Session, model: Type[Any], expected: int) -> int:
    """Return row count for model and raise if it doesn't match expected."""
    tablename = getattr(model, "__tablename__", None)
    if not tablename or not _table_exists(session, tablename):
        raise RuntimeError(f"Table missing for {model}: {tablename}")
    actual = int(session.query(model).count())
    if actual != expected:
        raise AssertionError(f"Count mismatch for {tablename}: expected={expected}, actual={actual}")
    return actual


def safe_count(session: Session, model: Type[Any]) -> Optional[int]:
    """Return row count if table exists; otherwise None."""
    tablename = getattr(model, "__tablename__", None)
    if not tablename or not _table_exists(session, tablename):
        return None
    try:
        return int(session.query(model).count())
    except (OperationalError, ProgrammingError):
        return None


