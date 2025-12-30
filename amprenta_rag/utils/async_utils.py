"""Async utilities for wrapping sync functions."""

import asyncio
from functools import wraps
from typing import Any, Callable, Coroutine, List, TypeVar

T = TypeVar('T')


def run_sync(func: Callable[..., T]) -> Callable[..., Coroutine[Any, Any, T]]:
    """Decorator to run sync function in thread pool."""
    @wraps(func)
    async def wrapper(*args, **kwargs):
        return await asyncio.to_thread(func, *args, **kwargs)
    return wrapper


async def gather_with_limit(n: int, *coros: Coroutine) -> List[Any]:
    """Run coroutines with max concurrency limit."""
    semaphore = asyncio.Semaphore(n)
    
    async def sem_coro(coro):
        async with semaphore:
            return await coro
    
    return await asyncio.gather(*(sem_coro(c) for c in coros))
