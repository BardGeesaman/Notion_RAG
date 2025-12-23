"""Binding site detection (fpocket) API endpoints."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel

from amprenta_rag.database.models import BindingSite, ProteinStructure
from amprenta_rag.database.session import db_session
from amprenta_rag.structural.fpocket_parser import parse_fpocket_output
from amprenta_rag.structural.fpocket_runner import run_fpocket


router = APIRouter(prefix="/pockets", tags=["Pockets"])


class DetectPocketsResponse(BaseModel):
    ok: bool
    pockets_created: int


class BindingSiteResponse(BaseModel):
    id: UUID
    structure_id: UUID
    pocket_rank: int
    score: Optional[float] = None
    volume: Optional[float] = None
    center_x: Optional[float] = None
    center_y: Optional[float] = None
    center_z: Optional[float] = None
    residues: Optional[List[str]] = None
    pocket_pdb_path: Optional[str] = None
    detection_method: str

    class Config:
        from_attributes = True


@router.post("/structures/{structure_id}/detect-pockets", response_model=DetectPocketsResponse)
def detect_pockets(structure_id: UUID) -> DetectPocketsResponse:
    with db_session() as db:
        s = db.query(ProteinStructure).filter(ProteinStructure.id == structure_id).first()
        if not s:
            raise HTTPException(status_code=404, detail="Structure not found")

        raw = next((f for f in (s.files or []) if f.file_type in ("pdb", "prepared")), None)
        if raw is None:
            raise HTTPException(status_code=400, detail="No PDB file available for structure")

        in_path = raw.file_path
        out_dir = str(Path(in_path).parent / "fpocket_out")
        ok = run_fpocket(in_path, out_dir)
        if not ok:
            raise HTTPException(status_code=500, detail="fpocket execution failed")

        pockets = parse_fpocket_output(out_dir)
        created = 0
        for p in pockets:
            # Upsert by unique key (structure_id, method, rank)
            existing = (
                db.query(BindingSite)
                .filter_by(structure_id=structure_id, detection_method="fpocket", pocket_rank=p.pocket_rank)
                .first()
            )
            if existing is None:
                bs = BindingSite(
                    structure_id=structure_id,
                    pocket_rank=p.pocket_rank,
                    score=p.score,
                    volume=p.volume,
                    center_x=p.center_x,
                    center_y=p.center_y,
                    center_z=p.center_z,
                    residues=p.residues,
                    pocket_pdb_path=p.pocket_pdb_path,
                    detection_method="fpocket",
                )
                db.add(bs)
                created += 1
            else:
                existing.score = p.score
                existing.volume = p.volume
                existing.center_x = p.center_x
                existing.center_y = p.center_y
                existing.center_z = p.center_z
                existing.pocket_pdb_path = p.pocket_pdb_path
        db.commit()

        return DetectPocketsResponse(ok=True, pockets_created=created)


@router.get("", response_model=List[BindingSiteResponse])
def list_pockets(structure_id: UUID = Query(...)) -> List[BindingSiteResponse]:
    with db_session() as db:
        rows = (
            db.query(BindingSite)
            .filter(BindingSite.structure_id == structure_id)
            .order_by(BindingSite.pocket_rank.asc())
            .all()
        )
        return [BindingSiteResponse.model_validate(r) for r in rows]


@router.get("/{pocket_id}/download")
def download_pocket(pocket_id: UUID):
    with db_session() as db:
        p = db.query(BindingSite).filter(BindingSite.id == pocket_id).first()
        if not p:
            raise HTTPException(status_code=404, detail="Pocket not found")
        if not p.pocket_pdb_path:
            raise HTTPException(status_code=404, detail="Pocket PDB not available")
        path = Path(p.pocket_pdb_path)
        if not path.exists():
            raise HTTPException(status_code=404, detail="Pocket file missing on disk")
        return FileResponse(str(path), filename=path.name, media_type="chemical/x-pdb")


__all__ = ["router"]


