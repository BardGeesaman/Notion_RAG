"""Lipidomics spectral matching API endpoints."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, ConfigDict, Field

from amprenta_rag.database.models import LipidAnnotation, SpectralLibrary, SpectralReference
from amprenta_rag.database.session import db_session
from amprenta_rag.spectral.library_loader import load_library
from amprenta_rag.spectral.matching_service import match_feature


router = APIRouter(prefix="/spectral", tags=["Spectral"])


class SpectralLibraryResponse(BaseModel):
    id: UUID
    name: str
    version: Optional[str] = None
    source_url: Optional[str] = None
    file_path: Optional[str] = None
    n_spectra: Optional[int] = None

    model_config = ConfigDict(from_attributes=True)


class IngestLibraryRequest(BaseModel):
    mgf_path: str
    name: str
    version: Optional[str] = None
    source_url: Optional[str] = None


class AnnotationResponse(BaseModel):
    id: UUID
    feature_id: UUID
    spectral_reference_id: UUID
    rank: Optional[int] = None
    spectral_score: float
    mz_error_ppm: Optional[float] = None
    matched_peaks: Optional[int] = None
    total_peaks: Optional[int] = None
    is_confident: bool
    is_ambiguous: bool
    manually_reviewed: bool
    review_status: Optional[str] = None
    reviewed_at: Optional[str] = None

    lipid_name: Optional[str] = None
    lipid_class: Optional[str] = None


class ReviewRequest(BaseModel):
    status: str = Field(..., description="confirmed|rejected")


@router.get("/libraries", response_model=List[SpectralLibraryResponse])
def list_libraries(limit: int = 100) -> List[SpectralLibraryResponse]:
    with db_session() as db:
        rows = db.query(SpectralLibrary).order_by(SpectralLibrary.created_at.desc()).limit(limit).all()
        return [SpectralLibraryResponse.model_validate(r) for r in rows]


@router.post("/libraries/ingest", response_model=SpectralLibraryResponse)
def ingest_library(payload: IngestLibraryRequest) -> SpectralLibraryResponse:
    lib = load_library(payload.mgf_path, payload.name, version=payload.version, source_url=payload.source_url)
    return SpectralLibraryResponse.model_validate(lib)


@router.post("/match/{feature_id}")
def run_match(feature_id: UUID) -> dict:
    with db_session() as db:
        try:
            anns = match_feature(feature_id, db)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
        return {"matched": len(anns)}


@router.get("/annotations/{feature_id}", response_model=List[AnnotationResponse])
def get_annotations(feature_id: UUID) -> List[AnnotationResponse]:
    with db_session() as db:
        rows = (
            db.query(LipidAnnotation, SpectralReference)
            .join(SpectralReference, SpectralReference.id == LipidAnnotation.spectral_reference_id)
            .filter(LipidAnnotation.feature_id == feature_id)
            .order_by(LipidAnnotation.rank.asc().nullslast(), LipidAnnotation.spectral_score.desc())
            .all()
        )
        out: List[AnnotationResponse] = []
        for ann, ref in rows:
            out.append(
                AnnotationResponse(
                    id=ann.id,
                    feature_id=ann.feature_id,
                    spectral_reference_id=ann.spectral_reference_id,
                    rank=ann.rank,
                    spectral_score=ann.spectral_score,
                    mz_error_ppm=ann.mz_error_ppm,
                    matched_peaks=ann.matched_peaks,
                    total_peaks=ann.total_peaks,
                    is_confident=bool(ann.is_confident),
                    is_ambiguous=bool(ann.is_ambiguous),
                    manually_reviewed=bool(ann.manually_reviewed),
                    review_status=ann.review_status,
                    reviewed_at=ann.reviewed_at.isoformat() if ann.reviewed_at else None,
                    lipid_name=ref.lipid_name,
                    lipid_class=ref.lipid_class,
                )
            )
        return out


@router.post("/annotations/{annotation_id}/review", response_model=AnnotationResponse)
def review_annotation(annotation_id: UUID, payload: ReviewRequest) -> AnnotationResponse:
    status = payload.status.strip().lower()
    if status not in ("confirmed", "rejected"):
        raise HTTPException(status_code=400, detail="status must be confirmed|rejected")

    with db_session() as db:
        ann = db.query(LipidAnnotation).filter_by(id=annotation_id).first()
        if not ann:
            raise HTTPException(status_code=404, detail="Annotation not found")
        ann.manually_reviewed = True
        ann.review_status = status
        ann.reviewed_at = datetime.now(timezone.utc)
        db.add(ann)

        ref = db.query(SpectralReference).filter_by(id=ann.spectral_reference_id).first()
        return AnnotationResponse(
            id=ann.id,
            feature_id=ann.feature_id,
            spectral_reference_id=ann.spectral_reference_id,
            rank=ann.rank,
            spectral_score=ann.spectral_score,
            mz_error_ppm=ann.mz_error_ppm,
            matched_peaks=ann.matched_peaks,
            total_peaks=ann.total_peaks,
            is_confident=bool(ann.is_confident),
            is_ambiguous=bool(ann.is_ambiguous),
            manually_reviewed=bool(ann.manually_reviewed),
            review_status=ann.review_status,
            reviewed_at=ann.reviewed_at.isoformat() if ann.reviewed_at else None,
            lipid_name=ref.lipid_name if ref else None,
            lipid_class=ref.lipid_class if ref else None,
        )


__all__ = ["router"]


