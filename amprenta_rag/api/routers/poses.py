"""Pose QC + interactions API endpoints."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, ConfigDict
from sqlalchemy.orm import selectinload, joinedload

from amprenta_rag.database.models import DockingPose, DockingRun, PoseInteraction, PoseQuality, ProteinStructure
from amprenta_rag.database.session import db_session
from amprenta_rag.structural.pose_qc import analyze_pose


router = APIRouter(prefix="/poses", tags=["Poses"])


class PoseQualityResponse(BaseModel):
    pose_id: UUID
    num_hbonds: int
    num_hydrophobic: int
    num_salt_bridges: int
    num_pi_stacking: int
    num_pi_cation: int
    num_halogen: int
    num_metal: int
    total_interactions: int
    has_clashes: bool
    ligand_efficiency: Optional[float] = None

    model_config = ConfigDict(from_attributes=True)


class PoseInteractionResponse(BaseModel):
    id: UUID
    pose_id: UUID
    interaction_type: str
    ligand_atom: Optional[str] = None
    protein_residue: Optional[str] = None
    distance: Optional[float] = None
    angle: Optional[float] = None

    model_config = ConfigDict(from_attributes=True)


def _get_receptor_pdb_path(structure: ProteinStructure) -> Optional[str]:
    files = list(structure.files or [])
    preferred = next((f for f in files if f.file_type == "prepared"), None) or next((f for f in files if f.file_type == "pdb"), None)
    return preferred.file_path if preferred else None


@router.post("/{pose_id}/analyze", response_model=PoseQualityResponse)
def analyze(pose_id: UUID) -> PoseQualityResponse:
    with db_session() as db:
        pose = db.query(DockingPose).options(
            joinedload(DockingPose.compound),
            joinedload(DockingPose.docking_run).joinedload(DockingRun.structure).selectinload(ProteinStructure.files)
        ).filter(DockingPose.id == pose_id).first()
        if not pose:
            raise HTTPException(status_code=404, detail="Pose not found")
        
        run = pose.docking_run
        if not run or not run.structure:
            raise HTTPException(status_code=400, detail="Pose missing docking run/structure")
        structure = run.structure

        receptor_pdb = _get_receptor_pdb_path(structure)
        if not receptor_pdb:
            raise HTTPException(status_code=400, detail="No receptor PDB available for structure")

        try:
            metrics, interactions = analyze_pose(pose, receptor_pdb=receptor_pdb)
        except RuntimeError as e:
            msg = str(e)
            if "obabel" in msg.lower():
                raise HTTPException(status_code=503, detail=msg)
            raise HTTPException(status_code=500, detail=msg)

        # Replace existing QC records
        db.query(PoseInteraction).filter(PoseInteraction.pose_id == pose_id).delete()
        db.query(PoseQuality).filter(PoseQuality.pose_id == pose_id).delete()

        q = PoseQuality(
            pose_id=pose_id,
            num_hbonds=metrics.num_hbonds,
            num_hydrophobic=metrics.num_hydrophobic,
            num_salt_bridges=metrics.num_salt_bridges,
            num_pi_stacking=metrics.num_pi_stacking,
            num_pi_cation=metrics.num_pi_cation,
            num_halogen=metrics.num_halogen,
            num_metal=metrics.num_metal,
            total_interactions=metrics.total_interactions,
            has_clashes=metrics.has_clashes,
            ligand_efficiency=metrics.ligand_efficiency,
        )
        db.add(q)
        for i in interactions:
            db.add(
                PoseInteraction(
                    pose_id=pose_id,
                    interaction_type=i.interaction_type,
                    ligand_atom=i.ligand_atom,
                    protein_residue=i.protein_residue,
                    distance=i.distance,
                    angle=i.angle,
                )
            )
        db.commit()
        db.refresh(q)

        return PoseQualityResponse(
            pose_id=q.pose_id,
            num_hbonds=q.num_hbonds,
            num_hydrophobic=q.num_hydrophobic,
            num_salt_bridges=q.num_salt_bridges,
            num_pi_stacking=q.num_pi_stacking,
            num_pi_cation=q.num_pi_cation,
            num_halogen=q.num_halogen,
            num_metal=q.num_metal,
            total_interactions=q.total_interactions,
            has_clashes=bool(q.has_clashes),
            ligand_efficiency=q.ligand_efficiency,
        )


@router.get("/{pose_id}/quality", response_model=PoseQualityResponse)
def get_quality(pose_id: UUID) -> PoseQualityResponse:
    with db_session() as db:
        q = db.query(PoseQuality).filter(PoseQuality.pose_id == pose_id).first()
        if not q:
            raise HTTPException(status_code=404, detail="Quality not found")
        return PoseQualityResponse(
            pose_id=q.pose_id,
            num_hbonds=q.num_hbonds,
            num_hydrophobic=q.num_hydrophobic,
            num_salt_bridges=q.num_salt_bridges,
            num_pi_stacking=q.num_pi_stacking,
            num_pi_cation=q.num_pi_cation,
            num_halogen=q.num_halogen,
            num_metal=q.num_metal,
            total_interactions=q.total_interactions,
            has_clashes=bool(q.has_clashes),
            ligand_efficiency=q.ligand_efficiency,
        )


@router.get("/{pose_id}/interactions", response_model=List[PoseInteractionResponse])
def get_interactions(pose_id: UUID) -> List[PoseInteractionResponse]:
    with db_session() as db:
        rows = (
            db.query(PoseInteraction)
            .filter(PoseInteraction.pose_id == pose_id)
            .order_by(PoseInteraction.interaction_type.asc())
            .all()
        )
        return [PoseInteractionResponse.model_validate(r) for r in rows]


__all__ = ["router"]


