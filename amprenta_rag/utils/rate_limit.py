"""Rate limiting utilities."""
from __future__ import annotations

import time
from typing import Dict, List, Tuple
from collections import defaultdict

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class RateLimitError(Exception):
    """Raised when rate limit is exceeded."""
    pass


class RateLimiter:
    """In-memory rate limiter using sliding window."""

    def __init__(self):
        # Structure: {user_id: {action: [(timestamp, count)]}}
        self._records: Dict[str, Dict[str, List[Tuple[float, int]]]] = defaultdict(lambda: defaultdict(list))
        self._cleanup_interval = 300  # Clean up old records every 5 minutes
        self._last_cleanup = time.time()

    def _cleanup_old_records(self, window_seconds: int) -> None:
        """Remove records older than the window."""
        current_time = time.time()
        if current_time - self._last_cleanup < self._cleanup_interval:
            return

        cutoff_time = current_time - window_seconds

        for user_id in list(self._records.keys()):
            for action in list(self._records[user_id].keys()):
                self._records[user_id][action] = [
                    (ts, count) for ts, count in self._records[user_id][action]
                    if ts > cutoff_time
                ]
                # Remove empty action dicts
                if not self._records[user_id][action]:
                    del self._records[user_id][action]
            # Remove empty user dicts
            if not self._records[user_id]:
                del self._records[user_id]

        self._last_cleanup = current_time

    def check_rate_limit(
        self,
        user_id: str,
        action: str,
        limit: int,
        window_seconds: int,
    ) -> bool:
        """
        Check if an action is allowed under rate limit.

        Args:
            user_id: User identifier
            action: Action identifier (e.g., "query_rag", "import_data")
            limit: Maximum number of actions allowed in the window
            window_seconds: Time window in seconds

        Returns:
            True if allowed, False if rate limit exceeded
        """
        current_time = time.time()
        cutoff_time = current_time - window_seconds

        # Clean up old records periodically
        self._cleanup_old_records(window_seconds)

        # Get records for this user/action within the window
        records = self._records[user_id][action]
        recent_records = [(ts, count) for ts, count in records if ts > cutoff_time]

        # Count total actions in window
        total_count = sum(count for _, count in recent_records)

        if total_count >= limit:
            logger.warning("[RATE_LIMIT] Rate limit exceeded for user %s, action %s (%d/%d)",
                          user_id, action, total_count, limit)
            return False

        # Record this action
        recent_records.append((current_time, 1))
        self._records[user_id][action] = recent_records

        return True

    def get_remaining(
        self,
        user_id: str,
        action: str,
        limit: int,
        window_seconds: int,
    ) -> int:
        """
        Get remaining actions allowed in the current window.

        Args:
            user_id: User identifier
            action: Action identifier
            limit: Maximum number of actions allowed
            window_seconds: Time window in seconds

        Returns:
            Number of remaining actions
        """
        current_time = time.time()
        cutoff_time = current_time - window_seconds

        records = self._records[user_id][action]
        recent_records = [(ts, count) for ts, count in records if ts > cutoff_time]

        total_count = sum(count for _, count in recent_records)
        remaining = max(0, limit - total_count)

        return remaining


# Singleton instance
_rate_limiter: RateLimiter | None = None


def get_rate_limiter() -> RateLimiter:
    """Get the singleton RateLimiter instance."""
    global _rate_limiter
    if _rate_limiter is None:
        _rate_limiter = RateLimiter()
    return _rate_limiter
