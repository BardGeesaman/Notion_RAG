"""Spectral matching service for lipid annotations."""

from __future__ import annotations

from typing import Any, Dict, List, Tuple
from uuid import UUID

from amprenta_rag.database.models import Feature, LipidAnnotation, SpectralReference
from amprenta_rag.spectral.matcher import MatchResult, match_spectrum


def _load_feature_spectrum(feature: Feature) -> tuple[float, List[Tuple[float, float]]]:
    """Load query spectrum from Feature.external_ids (MVP convention).

    Expected structure:
    feature.external_ids["spectrum"] = {
      "precursor_mz": float,
      "peaks": [[mz, intensity], ...]
    }
    """
    ext = getattr(feature, "external_ids", None) or {}
    if not isinstance(ext, dict):
        raise ValueError("Feature.external_ids missing (expected dict with spectrum)")
    spec = ext.get("spectrum") or {}
    if not isinstance(spec, dict):
        raise ValueError("Feature.external_ids.spectrum missing/invalid")
    precursor_mz = float(spec.get("precursor_mz"))
    peaks_raw = spec.get("peaks") or []
    peaks: List[Tuple[float, float]] = []
    for p in peaks_raw:
        try:
            peaks.append((float(p[0]), float(p[1])))
        except Exception:
            continue
    if not peaks:
        raise ValueError("Feature spectrum peaks empty")
    return precursor_mz, peaks


def match_feature(
    feature_id: UUID,
    db,
    *,
    precursor_tolerance_ppm: float = 10.0,
    fragment_tolerance_da: float = 0.01,
    confident_score: float = 0.7,
    confident_ppm: float = 10.0,
) -> List[LipidAnnotation]:
    """Match a feature spectrum to reference spectra and store LipidAnnotation rows."""
    feat = db.query(Feature).filter(Feature.id == feature_id).first()
    if not feat:
        raise ValueError("Feature not found")

    query_mz, query_peaks = _load_feature_spectrum(feat)

    refs = db.query(SpectralReference).all()
    ref_dicts: List[Dict[str, Any]] = []
    for r in refs:
        ref_dicts.append(
            {
                "id": str(r.id),
                "lipid_name": r.lipid_name,
                "precursor_mz": r.precursor_mz,
                "spectrum": r.spectrum,
            }
        )

    matches: List[MatchResult] = match_spectrum(
        query_mz,
        query_peaks,
        ref_dicts,
        precursor_tolerance_ppm=precursor_tolerance_ppm,
        fragment_tolerance_da=fragment_tolerance_da,
        top_k=10,
    )

    # Remove existing annotations (recompute)
    db.query(LipidAnnotation).filter(LipidAnnotation.feature_id == feature_id).delete()

    out: List[LipidAnnotation] = []
    # Determine ambiguity: multiple confident hits close in score
    confident = [m for m in matches if (m.score >= confident_score and m.mz_error_ppm <= confident_ppm)]
    ambiguous = False
    if len(confident) >= 2:
        if abs(confident[0].score - confident[1].score) <= 0.05:
            ambiguous = True

    for idx, m in enumerate(matches, start=1):
        is_conf = bool(m.score >= confident_score and m.mz_error_ppm <= confident_ppm)
        ann = LipidAnnotation(
            feature_id=feature_id,
            spectral_reference_id=UUID(m.spectral_reference_id),
            spectral_score=float(m.score),
            mz_error_ppm=float(m.mz_error_ppm),
            matched_peaks=int(m.matched_peaks),
            total_peaks=int(m.total_peaks),
            is_confident=is_conf,
            is_ambiguous=bool(ambiguous and is_conf),
            rank=int(idx),
            manually_reviewed=False,
            review_status=None,
            reviewed_at=None,
        )
        db.add(ann)
        out.append(ann)

    db.commit()
    return out


__all__ = ["match_feature"]


