"""Account lockout tracking for brute force protection.

Uses database for persistence (P1 fix: survives restarts).
Uses datetime.now(timezone.utc) per Python 3.12+ (P1 fix).
"""
from datetime import datetime, timedelta, timezone
from typing import Optional, Tuple
import logging
import os

from sqlalchemy import Column, Integer, String, DateTime, func
from sqlalchemy.orm import Session

from amprenta_rag.database.base import Base
from amprenta_rag.database.session import db_session

logger = logging.getLogger(__name__)

MAX_FAILED_ATTEMPTS = int(os.getenv("LOCKOUT_MAX_ATTEMPTS", "5"))
LOCKOUT_DURATION = timedelta(minutes=int(os.getenv("LOCKOUT_DURATION_MINUTES", "15")))
ATTEMPT_WINDOW = timedelta(minutes=5)


class AuthLockoutAttempt(Base):
    """Track failed authentication attempts."""
    __tablename__ = 'auth_lockout_attempts'
    
    id = Column(Integer, primary_key=True)
    identifier = Column(String(255), nullable=False, index=True)
    attempt_time = Column(DateTime(timezone=True), nullable=False)


def record_failed_attempt(identifier: str, db: Optional[Session] = None) -> int:
    """Record a failed authentication attempt.
    
    Args:
        identifier: Unique identifier (e.g., "sign:192.168.1.1")
        db: Optional session (creates new if not provided)
    
    Returns:
        Current count of failed attempts in window
    """
    now = datetime.now(timezone.utc)
    
    def _record(session: Session) -> int:
        # Add new attempt
        attempt = AuthLockoutAttempt(identifier=identifier, attempt_time=now)
        session.add(attempt)
        session.commit()
        
        # Clean old attempts
        cutoff = now - ATTEMPT_WINDOW
        session.query(AuthLockoutAttempt).filter(
            AuthLockoutAttempt.identifier == identifier,
            AuthLockoutAttempt.attempt_time < cutoff
        ).delete()
        session.commit()
        
        # Get current count
        count = session.query(AuthLockoutAttempt).filter(
            AuthLockoutAttempt.identifier == identifier,
            AuthLockoutAttempt.attempt_time >= cutoff
        ).count()
        
        logger.warning(f"Failed auth attempt for {identifier}: {count}/{MAX_FAILED_ATTEMPTS}")
        return count
    
    if db:
        return _record(db)
    else:
        with db_session() as session:
            return _record(session)


def is_locked_out(identifier: str, db: Optional[Session] = None) -> Tuple[bool, Optional[int]]:
    """Check if identifier is locked out.
    
    Args:
        identifier: Unique identifier to check
        db: Optional session
    
    Returns:
        (is_locked, seconds_remaining) - seconds_remaining is None if not locked
    """
    now = datetime.now(timezone.utc)
    cutoff = now - ATTEMPT_WINDOW
    
    def _check(session: Session) -> Tuple[bool, Optional[int]]:
        count = session.query(AuthLockoutAttempt).filter(
            AuthLockoutAttempt.identifier == identifier,
            AuthLockoutAttempt.attempt_time >= cutoff
        ).count()
        
        if count >= MAX_FAILED_ATTEMPTS:
            # Get last attempt time
            last = session.query(func.max(AuthLockoutAttempt.attempt_time)).filter(
                AuthLockoutAttempt.identifier == identifier
            ).scalar()
            
            if last:
                # Handle timezone-aware comparison
                if last.tzinfo is None:
                    last = last.replace(tzinfo=timezone.utc)
                lockout_end = last + LOCKOUT_DURATION
                if now < lockout_end:
                    remaining = int((lockout_end - now).total_seconds())
                    return True, remaining
        
        return False, None
    
    if db:
        return _check(db)
    else:
        with db_session() as session:
            return _check(session)


def clear_failed_attempts(identifier: str, db: Optional[Session] = None) -> None:
    """Clear failed attempts after successful auth.
    
    Args:
        identifier: Identifier to clear
        db: Optional session
    """
    def _clear(session: Session):
        deleted = session.query(AuthLockoutAttempt).filter(
            AuthLockoutAttempt.identifier == identifier
        ).delete()
        session.commit()
        if deleted:
            logger.info(f"Cleared {deleted} failed attempts for {identifier}")
    
    if db:
        _clear(db)
    else:
        with db_session() as session:
            _clear(session)
