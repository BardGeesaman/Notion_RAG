"""3D visualization API endpoints (conformers, overlays, protein PDB retrieval)."""

from __future__ import annotations

import asyncio
from pathlib import Path
from typing import Any, Dict
from uuid import UUID

from fastapi import APIRouter, HTTPException

from amprenta_rag.api.schemas import ConformerRequest, ConformerResponse, OverlayRequest, OverlayResponse
from amprenta_rag.chemistry.conformers import (
    RDKIT_AVAILABLE,
    align_molecules_to_reference,
    conformer_energies,
    conformer_to_pdb,
    generate_conformers,
    optimize_conformer,
)
from amprenta_rag.database.models import ProteinStructure, StructureFile
from amprenta_rag.database.session import db_session


router = APIRouter(prefix="/viz3d", tags=["viz3d"])


# Sync helper functions for compute-intensive operations
def _sync_get_conformers(payload: ConformerRequest) -> ConformerResponse:
    """Sync helper for RDKit conformer generation."""
    if not RDKIT_AVAILABLE:
        raise HTTPException(status_code=400, detail="RDKit not installed; 3D conformers unavailable")

    try:
        mols = generate_conformers(payload.smiles, n_conformers=int(payload.n_conformers), method="ETKDG")
        if payload.optimize:
            for m in mols:
                optimize_conformer(m, force_field="MMFF")
        pdbs = [conformer_to_pdb(m, conf_id=0) for m in mols]
        ens = [conformer_energies(m, prefer="MMFF")[0] for m in mols]
        return ConformerResponse(pdb_strings=pdbs, energies=ens)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except TimeoutError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=str(e))


def _sync_overlay(payload: OverlayRequest) -> OverlayResponse:
    """Sync helper for RDKit molecular overlay alignment."""
    if not RDKIT_AVAILABLE:
        raise HTTPException(status_code=400, detail="RDKit not installed; 3D overlays unavailable")
    if not payload.smiles_list:
        raise HTTPException(status_code=400, detail="smiles_list must be non-empty")

    try:
        mols = []
        for smi in payload.smiles_list:
            ms = generate_conformers(smi, n_conformers=1, method="ETKDG")
            m = ms[0]
            optimize_conformer(m, force_field="MMFF")
            mols.append(m)

        align_molecules_to_reference(mols, reference_idx=int(payload.reference_idx))
        pdbs = [conformer_to_pdb(m, conf_id=0) for m in mols]
        return OverlayResponse(aligned_pdb_strings=pdbs)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except TimeoutError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=str(e))


def _sync_get_protein_pdb(structure_id: UUID, file_type: str = "pdb") -> Dict[str, Any]:
    """Sync helper for protein PDB file I/O operations."""
    ft = str(file_type or "pdb").lower()
    with db_session() as db:
        stc = db.query(ProteinStructure).filter(ProteinStructure.id == structure_id).first()
        if not stc:
            raise HTTPException(status_code=404, detail="ProteinStructure not found")
        q = db.query(StructureFile).filter(StructureFile.structure_id == structure_id)
        if ft:
            q = q.filter(StructureFile.file_type == ft)
        sf = q.order_by(StructureFile.created_at.desc()).first()
        if not sf:
            # Fallback: any file
            sf = (
                db.query(StructureFile)
                .filter(StructureFile.structure_id == structure_id)
                .order_by(StructureFile.created_at.desc())
                .first()
            )
        if not sf:
            raise HTTPException(status_code=404, detail="No StructureFile found for structure")

        fp = Path(str(sf.file_path))
        if not fp.exists():
            raise HTTPException(status_code=404, detail=f"Structure file not found on disk: {sf.file_path}")
        try:
            pdb_str = fp.read_text(encoding="utf-8", errors="replace")
        except Exception as e:  # noqa: BLE001
            raise HTTPException(status_code=500, detail=f"Failed to read structure file: {e}")

        meta = {
            "structure_id": str(structure_id),
            "pdb_id": stc.pdb_id,
            "alphafold_uniprot_id": stc.alphafold_uniprot_id,
            "source": stc.source,
            "file_id": str(sf.id),
            "file_type": sf.file_type,
            "file_path": sf.file_path,
            "file_size_bytes": sf.file_size_bytes,
        }
        return {"pdb_string": pdb_str, "metadata": meta}


@router.post("/conformers", response_model=ConformerResponse)
async def get_conformers(payload: ConformerRequest) -> ConformerResponse:
    """Generate 3D conformers using async thread pool for RDKit computations."""
    return await asyncio.to_thread(_sync_get_conformers, payload)


@router.post("/overlay", response_model=OverlayResponse)
async def overlay(payload: OverlayRequest) -> OverlayResponse:
    """Perform molecular overlay alignment using async thread pool for RDKit computations."""
    return await asyncio.to_thread(_sync_overlay, payload)


@router.get("/protein/{structure_id}")
async def get_protein_pdb(structure_id: UUID, file_type: str = "pdb") -> Dict[str, Any]:
    """Return the PDB string for a ProteinStructure using async thread pool for file I/O."""
    return await asyncio.to_thread(_sync_get_protein_pdb, structure_id, file_type)


__all__ = ["router"]


