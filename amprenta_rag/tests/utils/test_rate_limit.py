from __future__ import annotations

import time

from amprenta_rag.utils import rate_limit


def test_rate_limiter_allows_then_blocks(monkeypatch):
    rl = rate_limit.RateLimiter()
    base = time.time()
    times = [base, base + 1, base + 2, base + 3]
    monkeypatch.setattr(time, "time", lambda: times.pop(0))

    assert rl.check_rate_limit("u", "a", limit=2, window_seconds=5) is True
    assert rl.check_rate_limit("u", "a", limit=2, window_seconds=5) is True
    assert rl.check_rate_limit("u", "a", limit=2, window_seconds=5) is False


def test_rate_limiter_remaining(monkeypatch):
    rl = rate_limit.RateLimiter()
    t0 = time.time()
    monkeypatch.setattr(time, "time", lambda: t0)
    rl.check_rate_limit("u", "a", limit=3, window_seconds=10)
    monkeypatch.setattr(time, "time", lambda: t0 + 1)
    remaining = rl.get_remaining("u", "a", limit=3, window_seconds=10)
    assert remaining == 2

