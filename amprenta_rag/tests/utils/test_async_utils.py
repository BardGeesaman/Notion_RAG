"""Tests for async utilities."""

import asyncio
import time
from unittest.mock import AsyncMock, MagicMock

import pytest

from amprenta_rag.utils.async_utils import run_sync, gather_with_limit


class TestRunSync:
    """Test run_sync decorator."""

    @pytest.mark.asyncio
    async def test_run_sync_basic(self):
        """Test decorator converts sync function to async."""
        def sync_func():
            return "hello"

        async_func = run_sync(sync_func)
        
        # Verify it returns a coroutine
        coro = async_func()
        assert asyncio.iscoroutine(coro)
        
        # Verify result is correct
        result = await coro
        assert result == "hello"

    @pytest.mark.asyncio
    async def test_run_sync_with_args(self):
        """Test arguments pass through correctly."""
        def sync_func(a, b, c=None):
            return f"{a}-{b}-{c}"

        async_func = run_sync(sync_func)
        
        # Test positional and keyword args
        result = await async_func("pos1", "pos2", c="keyword")
        assert result == "pos1-pos2-keyword"

    @pytest.mark.asyncio
    async def test_run_sync_exception(self):
        """Test exceptions propagate properly."""
        def sync_func():
            raise ValueError("test error")

        async_func = run_sync(sync_func)
        
        # Verify exception propagates
        with pytest.raises(ValueError, match="test error"):
            await async_func()

    @pytest.mark.asyncio
    async def test_run_sync_preserves_function_metadata(self):
        """Test decorator preserves function metadata."""
        def sync_func():
            """Original docstring."""
            return "result"

        async_func = run_sync(sync_func)
        
        # Verify metadata is preserved
        assert async_func.__name__ == "sync_func"
        assert async_func.__doc__ == "Original docstring."


class TestGatherWithLimit:
    """Test gather_with_limit function."""

    @pytest.mark.asyncio
    async def test_gather_with_limit_concurrency(self):
        """Test limits concurrent execution."""
        execution_times = []
        
        async def slow_task(task_id: int, duration: float):
            start = time.time()
            execution_times.append(("start", task_id, start))
            await asyncio.sleep(duration)
            end = time.time()
            execution_times.append(("end", task_id, end))
            return f"task_{task_id}"

        # Create 4 tasks with limit of 2
        tasks = [
            slow_task(1, 0.1),
            slow_task(2, 0.1), 
            slow_task(3, 0.1),
            slow_task(4, 0.1)
        ]
        
        start_time = time.time()
        results = await gather_with_limit(2, *tasks)
        total_time = time.time() - start_time
        
        # Verify results
        assert results == ["task_1", "task_2", "task_3", "task_4"]
        
        # Verify concurrency limit (should take ~0.2s, not ~0.1s)
        # With limit of 2, tasks 1&2 run first (0.1s), then 3&4 (another 0.1s)
        assert total_time >= 0.15  # Allow some overhead
        assert total_time < 0.5   # But not sequential execution

    @pytest.mark.asyncio
    async def test_gather_with_limit_ordering(self):
        """Test results maintain order."""
        async def numbered_task(n: int):
            # Add small random delay to ensure ordering isn't by completion time
            await asyncio.sleep(0.01 * (5 - n))  # Reverse order delays
            return n * 10

        tasks = [numbered_task(i) for i in range(1, 6)]
        results = await gather_with_limit(3, *tasks)
        
        # Results should maintain original order despite completion order
        assert results == [10, 20, 30, 40, 50]

    @pytest.mark.asyncio
    async def test_gather_with_limit_exception_handling(self):
        """Test exception handling in limited gather."""
        async def failing_task():
            await asyncio.sleep(0.01)
            raise ValueError("task failed")
        
        async def success_task():
            await asyncio.sleep(0.01)
            return "success"

        tasks = [success_task(), failing_task(), success_task()]
        
        # Exception should propagate and stop execution
        with pytest.raises(ValueError, match="task failed"):
            await gather_with_limit(2, *tasks)

    @pytest.mark.asyncio
    async def test_gather_with_limit_single_task(self):
        """Test with single task and limit of 1."""
        async def single_task():
            return "single_result"

        result = await gather_with_limit(1, single_task())
        assert result == ["single_result"]

    @pytest.mark.asyncio
    async def test_gather_with_limit_no_tasks(self):
        """Test with no tasks."""
        result = await gather_with_limit(5)
        assert result == []
