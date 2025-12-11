from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Callable, Iterable, List, Optional

from tqdm import tqdm

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def parallel_map(
    func: Callable[[Any], Any],
    items: Iterable[Any],
    max_workers: int = 4,
    desc: Optional[str] = None,
) -> List[Any]:
    items_list = list(items)
    results: List[Any] = [None] * len(items_list)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {executor.submit(func, item): idx for idx, item in enumerate(items_list)}
        iterator = as_completed(future_map)
        if desc:
            iterator = tqdm(iterator, total=len(items_list), desc=desc)
        for fut in iterator:
            idx = future_map[fut]
            try:
                results[idx] = fut.result()
            except Exception as exc:
                logger.warning("[PARALLEL] Task %d failed: %r", idx, exc)
                results[idx] = None
    return results


def chunked_parallel(
    func: Callable[[Any], Any],
    items: Iterable[Any],
    chunk_size: int = 100,
    max_workers: int = 4,
    desc: Optional[str] = None,
) -> List[Any]:
    """Process items in chunks, applying func to each element within the chunk.

    This prevents overwhelming the thread pool with a huge number of tasks while
    still parallelizing within each chunk.
    """
    items_list = list(items)
    chunks = [items_list[i : i + chunk_size] for i in range(0, len(items_list), chunk_size)]
    all_results: List[Any] = []

    for chunk_idx, chunk in enumerate(tqdm(chunks, desc=desc) if desc else chunks):
        try:
            chunk_results = parallel_map(func, chunk, max_workers=max_workers)
        except Exception as exc:
            logger.warning("[PARALLEL] Chunk %d failed: %r", chunk_idx, exc)
            chunk_results = [None] * len(chunk)
        all_results.extend(chunk_results)

    return all_results


__all__ = ["parallel_map", "chunked_parallel"]
