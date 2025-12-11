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
    func: Callable[[List[Any]], Any],
    items: Iterable[Any],
    chunk_size: int = 100,
    max_workers: int = 4,
    desc: Optional[str] = None,
) -> List[Any]:
    items_list = list(items)
    chunks = [items_list[i : i + chunk_size] for i in range(0, len(items_list), chunk_size)]
    results: List[Any] = [None] * len(chunks)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {executor.submit(func, chunk): idx for idx, chunk in enumerate(chunks)}
        iterator = as_completed(future_map)
        if desc:
            iterator = tqdm(iterator, total=len(chunks), desc=desc)
        for fut in iterator:
            idx = future_map[fut]
            try:
                results[idx] = fut.result()
            except Exception as exc:
                logger.warning("[PARALLEL] Chunk %d failed: %r", idx, exc)
                results[idx] = None
    return results


__all__ = ["parallel_map", "chunked_parallel"]
