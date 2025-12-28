"""Docking (Vina) API endpoints."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional
from uuid import UUID, uuid4

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field

from amprenta_rag.database.models import BindingSite, Compound, DockingPose, DockingRun, ProteinStructure
from amprenta_rag.database.session import db_session
from amprenta_rag.structural.docking_service import DockingService


router = APIRouter(prefix="/docking", tags=["Docking"])

DOCKING_SERVICE = DockingService(max_workers=4)


class CreateDockingRunRequest(BaseModel):
    structure_id: UUID
    binding_site_id: Optional[UUID] = None
    compound_ids: List[str] = Field(default_factory=list, description="List of Compound UUIDs or compound_id strings")
    exhaustiveness: int = Field(default=8, ge=1, le=64)


class DockingRunResponse(BaseModel):
    id: UUID
    structure_id: UUID
    binding_site_id: Optional[UUID] = None
    center_x: Optional[float] = None
    center_y: Optional[float] = None
    center_z: Optional[float] = None
    size_x: Optional[float] = None
    size_y: Optional[float] = None
    size_z: Optional[float] = None
    status: str
    total_compounds: Optional[int] = None
    completed_compounds: int
    progress: float
    error_log: Optional[str] = None

    model_config = ConfigDict(from_attributes=True)


class DockingPoseResponse(BaseModel):
    id: UUID
    docking_run_id: UUID
    compound_id: UUID
    binding_affinity: Optional[float] = None
    pose_rank: Optional[int] = None
    rmsd_lb: Optional[float] = None
    rmsd_ub: Optional[float] = None
    pose_pdbqt_path: Optional[str] = None

    model_config = ConfigDict(from_attributes=True)


def _progress(run: DockingRun) -> float:
    total = int(run.total_compounds or 0)
    done = int(run.completed_compounds or 0)
    if total <= 0:
        return 0.0
    return max(0.0, min(1.0, done / total))


@router.post("/runs", response_model=DockingRunResponse)
def create_run(payload: CreateDockingRunRequest) -> DockingRunResponse:
    with db_session() as db:
        s = db.query(ProteinStructure).filter_by(id=payload.structure_id).first()
        if not s:
            raise HTTPException(status_code=404, detail="Structure not found")

        pocket = None
        if payload.binding_site_id:
            pocket = db.query(BindingSite).filter_by(id=payload.binding_site_id).first()
            if not pocket:
                raise HTTPException(status_code=404, detail="Binding site not found")
            if pocket.structure_id != payload.structure_id:
                raise HTTPException(status_code=400, detail="Binding site does not belong to structure")

        # Resolve compound identifiers (UUID or compound_id string)
        resolved: List[UUID] = []
        missing: List[str] = []
        for cid in payload.compound_ids or []:
            cid_s = str(cid).strip()
            if not cid_s:
                continue
            try:
                resolved.append(UUID(cid_s))
                continue
            except Exception:
                pass
            c = db.query(Compound).filter(Compound.compound_id == cid_s).first()
            if c:
                resolved.append(c.id)
            else:
                missing.append(cid_s)

        if missing:
            raise HTTPException(status_code=400, detail=f"Unknown compounds: {missing[:10]}")
        if not resolved:
            raise HTTPException(status_code=400, detail="compound_ids must be non-empty")

        run = DockingRun(
            id=uuid4(),
            structure_id=payload.structure_id,
            binding_site_id=payload.binding_site_id,
            center_x=(pocket.center_x if pocket else None),
            center_y=(pocket.center_y if pocket else None),
            center_z=(pocket.center_z if pocket else None),
            size_x=20.0,
            size_y=20.0,
            size_z=20.0,
            status="pending",
            total_compounds=len(resolved),
            completed_compounds=0,
            compound_ids=[str(x) for x in resolved],
            exhaustiveness=int(payload.exhaustiveness or 8),
        )
        db.add(run)
        db.commit()
        db.refresh(run)

    DOCKING_SERVICE.start_run(run.id)

    with db_session() as db:
        rr = db.query(DockingRun).filter_by(id=run.id).first()
        if not rr:
            raise HTTPException(status_code=404, detail="Run not found")
        return DockingRunResponse(
            id=rr.id,
            structure_id=rr.structure_id,
            binding_site_id=rr.binding_site_id,
            center_x=rr.center_x,
            center_y=rr.center_y,
            center_z=rr.center_z,
            size_x=rr.size_x,
            size_y=rr.size_y,
            size_z=rr.size_z,
            status=rr.status,
            total_compounds=rr.total_compounds,
            completed_compounds=int(rr.completed_compounds or 0),
            progress=_progress(rr),
            error_log=rr.error_log,
        )


@router.get("/runs", response_model=List[DockingRunResponse])
def list_runs(
    status: Optional[str] = Query(None),
    limit: int = Query(100, ge=1, le=1000),
) -> List[DockingRunResponse]:
    with db_session() as db:
        q = db.query(DockingRun).order_by(DockingRun.started_at.desc().nullslast())
        if status:
            q = q.filter(DockingRun.status == status)
        rows = q.limit(limit).all()
        return [
            DockingRunResponse(
                id=r.id,
                structure_id=r.structure_id,
                binding_site_id=r.binding_site_id,
                center_x=r.center_x,
                center_y=r.center_y,
                center_z=r.center_z,
                size_x=r.size_x,
                size_y=r.size_y,
                size_z=r.size_z,
                status=r.status,
                total_compounds=r.total_compounds,
                completed_compounds=int(r.completed_compounds or 0),
                progress=_progress(r),
                error_log=r.error_log,
            )
            for r in rows
        ]


@router.get("/runs/{run_id}", response_model=DockingRunResponse)
def get_run(run_id: UUID) -> DockingRunResponse:
    with db_session() as db:
        r = db.query(DockingRun).filter_by(id=run_id).first()
        if not r:
            raise HTTPException(status_code=404, detail="Run not found")
        return DockingRunResponse(
            id=r.id,
            structure_id=r.structure_id,
            binding_site_id=r.binding_site_id,
            center_x=r.center_x,
            center_y=r.center_y,
            center_z=r.center_z,
            size_x=r.size_x,
            size_y=r.size_y,
            size_z=r.size_z,
            status=r.status,
            total_compounds=r.total_compounds,
            completed_compounds=int(r.completed_compounds or 0),
            progress=_progress(r),
            error_log=r.error_log,
        )


@router.get("/runs/{run_id}/poses", response_model=List[DockingPoseResponse])
def list_run_poses(run_id: UUID) -> List[DockingPoseResponse]:
    with db_session() as db:
        rows = (
            db.query(DockingPose)
            .filter(DockingPose.docking_run_id == run_id)
            .order_by(DockingPose.binding_affinity.asc().nullslast())
            .all()
        )
        return [DockingPoseResponse.model_validate(r) for r in rows]


@router.get("/poses/{pose_id}/download")
def download_pose(pose_id: UUID):
    with db_session() as db:
        p = db.query(DockingPose).filter_by(id=pose_id).first()
        if not p:
            raise HTTPException(status_code=404, detail="Pose not found")
        if not p.pose_pdbqt_path:
            raise HTTPException(status_code=404, detail="Pose file not available")
        path = Path(p.pose_pdbqt_path)
        if not path.exists():
            raise HTTPException(status_code=404, detail="Pose file missing on disk")
        return FileResponse(str(path), filename=path.name, media_type="chemical/x-pdbqt")


__all__ = ["router"]


