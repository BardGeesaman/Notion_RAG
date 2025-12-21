from __future__ import annotations

from amprenta_rag.utils import parallel


def test_parallel_map_success_and_failure():
    items = [1, 2, 0]

    def func(x):
        if x == 0:
            raise ValueError("boom")
        return x * 2

    results = parallel.parallel_map(func, items, max_workers=2)
    assert results[:2] == [2, 4]
    assert results[2] is None


def test_chunked_parallel_chunks():
    items = list(range(5))
    results = parallel.chunked_parallel(lambda x: x + 1, items, chunk_size=2, max_workers=2)
    assert results == [1, 2, 3, 4, 5]

