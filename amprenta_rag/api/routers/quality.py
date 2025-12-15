from __future__ import annotations

from typing import List

from fastapi import APIRouter, Query

from amprenta_rag.analysis.quality_watcher import (
    get_low_quality_datasets,
    get_quality_summary,
    scan_all_datasets,
)
from amprenta_rag.api import schemas

router = APIRouter()


@router.get(
    "/quality/summary",
    summary="Dataset quality summary",
    response_model=schemas.QualitySummary,
)
def quality_summary() -> schemas.QualitySummary:
    summary = get_quality_summary()
    return schemas.QualitySummary(**summary)


@router.get(
    "/quality/datasets",
    summary="Quality scores for all datasets",
    response_model=List[schemas.DatasetQualityReport],
)
def quality_datasets() -> List[schemas.DatasetQualityReport]:
    reports = scan_all_datasets()
    return [
        schemas.DatasetQualityReport(
            dataset_id=report.dataset_id,
            dataset_name=report.dataset_name,
            score=report.score,
            status=report.status,
            issues=report.issues,
            metrics=report.metrics,
        )
        for report in reports
    ]


@router.get(
    "/quality/alerts",
    summary="Low quality dataset alerts",
    response_model=List[schemas.DatasetQualityReport],
)
def quality_alerts(
    threshold: float = Query(50.0, ge=0, le=100, description="Score threshold"),
) -> List[schemas.DatasetQualityReport]:
    reports = get_low_quality_datasets(threshold=threshold)
    return [
        schemas.DatasetQualityReport(
            dataset_id=report.dataset_id,
            dataset_name=report.dataset_name,
            score=report.score,
            status=report.status,
            issues=report.issues,
            metrics=report.metrics,
        )
        for report in reports
    ]

