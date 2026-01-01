"""Spectrum matching utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Tuple


@dataclass(frozen=True)
class MatchResult:
    spectral_reference_id: str
    lipid_name: str
    score: float
    mz_error_ppm: float
    matched_peaks: int
    total_peaks: int


def _ppm_error(obs: float, ref: float) -> float:
    if ref == 0:
        return 0.0
    return (obs - ref) / ref * 1e6


def _normalize(peaks: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    if not peaks:
        return []
    norm = sum((i * i) for _, i in peaks) ** 0.5
    if norm <= 0:
        return [(mz, 0.0) for mz, _ in peaks]
    return [(mz, float(i) / norm) for mz, i in peaks]


def _greedy_cosine(
    query_peaks: List[Tuple[float, float]],
    ref_peaks: List[Tuple[float, float]],
    fragment_tolerance_da: float,
) -> tuple[float, int]:
    """Internal greedy cosine similarity (fallback if matchms not installed)."""
    qn = _normalize(query_peaks)
    rn = _normalize(ref_peaks)
    used = set()
    score = 0.0
    matched = 0
    for q_mz, q_i in qn:
        best_j = None
        best_prod = 0.0
        for j, (r_mz, r_i) in enumerate(rn):
            if j in used:
                continue
            if abs(q_mz - r_mz) <= fragment_tolerance_da:
                prod = q_i * r_i
                if prod > best_prod:
                    best_prod = prod
                    best_j = j
        if best_j is not None and best_prod > 0:
            used.add(best_j)
            score += best_prod
            matched += 1
    return float(score), matched


def match_spectrum(
    query_mz: float,
    query_peaks: List[Tuple[float, float]],
    references: List[Dict[str, Any]],
    *,
    precursor_tolerance_ppm: float = 10.0,
    fragment_tolerance_da: float = 0.01,
    top_k: int = 10,
) -> List[MatchResult]:
    """Match a query spectrum against a list of reference spectra dicts.

    Each reference dict must include:
    - id (UUID or str)
    - lipid_name (str)
    - precursor_mz (float)
    - spectrum {"mz": [...], "intensity": [...]}
    """
    if not query_peaks:
        return []

    # Precursor filter
    cand = []
    for r in references:
        try:
            precursor_mz = r.get("precursor_mz")
            if precursor_mz is None:
                continue
            pmz = float(precursor_mz)
        except Exception:
            continue
        ppm = abs(_ppm_error(query_mz, pmz))
        if ppm <= float(precursor_tolerance_ppm):
            cand.append((ppm, r))

    results: List[MatchResult] = []
    for ppm, r in cand:
        spec = r.get("spectrum") or {}
        mzs = spec.get("mz") or []
        intens = spec.get("intensity") or []
        ref_peaks = []
        for m, i in zip(mzs, intens):
            try:
                ref_peaks.append((float(m), float(i)))
            except Exception:
                continue
        score, matched = _greedy_cosine(query_peaks, ref_peaks, fragment_tolerance_da=float(fragment_tolerance_da))
        results.append(
            MatchResult(
                spectral_reference_id=str(r.get("id")),
                lipid_name=str(r.get("lipid_name") or "unknown"),
                score=float(score),
                mz_error_ppm=float(ppm),
                matched_peaks=int(matched),
                total_peaks=int(len(query_peaks)),
            )
        )

    results.sort(key=lambda x: x.score, reverse=True)
    return results[: int(top_k)]


__all__ = ["MatchResult", "match_spectrum"]


