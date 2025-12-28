"""Protein structure store API endpoints."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict

from amprenta_rag.database.models import ProteinStructure, StructureFile
from amprenta_rag.database.session import db_session
from amprenta_rag.structural.alphafold_fetcher import fetch_from_alphafold
from amprenta_rag.structural.pdb_fetcher import fetch_from_pdb
from amprenta_rag.structural.pdb_parser import parse_pdb_metadata, to_fasta
from amprenta_rag.structural.prep import prepare_structure
from amprenta_rag.structural.storage import save_structure_file


router = APIRouter(prefix="/structures", tags=["Structures"])

HAS_MULTIPART = False
try:
    import multipart  # type: ignore[import-not-found]  # noqa: F401

    from fastapi import File, Form, UploadFile  # noqa: E402

    HAS_MULTIPART = True
except Exception:
    HAS_MULTIPART = False


def _md5(content: bytes) -> str:
    h = hashlib.md5()  # noqa: S324 - checksum only
    h.update(content)
    return h.hexdigest()


class StructureFileResponse(BaseModel):
    id: UUID
    file_type: str
    file_path: str
    file_size_bytes: Optional[int] = None
    md5_hash: Optional[str] = None

    model_config = ConfigDict(from_attributes=True)


class ProteinStructureResponse(BaseModel):
    id: UUID
    feature_id: Optional[UUID] = None
    pdb_id: Optional[str] = None
    alphafold_uniprot_id: Optional[str] = None
    source: str
    resolution: Optional[float] = None
    method: Optional[str] = None
    chain_ids: Optional[List[str]] = None
    prep_status: Optional[str] = None
    prep_log: Optional[str] = None
    files: List[StructureFileResponse] = []

    model_config = ConfigDict(from_attributes=True)


class FetchStructureRequest(BaseModel):
    source: str  # pdb | alphafold
    pdb_id: Optional[str] = None
    uniprot_id: Optional[str] = None
    feature_id: Optional[UUID] = None
    chain_id: Optional[str] = None


class PrepareRequest(BaseModel):
    chain_id: Optional[str] = None


@router.get("", response_model=List[ProteinStructureResponse])
def list_structures(
    source: Optional[str] = Query(None),
    prep_status: Optional[str] = Query(None),
    limit: int = Query(100, ge=1, le=1000),
) -> List[ProteinStructureResponse]:
    with db_session() as db:
        q = db.query(ProteinStructure).order_by(ProteinStructure.created_at.desc())
        if source:
            q = q.filter(ProteinStructure.source == source)
        if prep_status:
            q = q.filter(ProteinStructure.prep_status == prep_status)
        rows = q.limit(limit).all()
        # eager-load files
        for r in rows:
            _ = list(r.files or [])
        return [ProteinStructureResponse.model_validate(r) for r in rows]


@router.get("/{structure_id}", response_model=ProteinStructureResponse)
def get_structure(structure_id: UUID) -> ProteinStructureResponse:
    with db_session() as db:
        s = db.query(ProteinStructure).filter(ProteinStructure.id == structure_id).first()
        if not s:
            raise HTTPException(status_code=404, detail="Structure not found")
        _ = list(s.files or [])
        return ProteinStructureResponse.model_validate(s)


@router.post("/fetch", response_model=ProteinStructureResponse)
def fetch_structure(payload: FetchStructureRequest) -> ProteinStructureResponse:
    src = (payload.source or "").lower().strip()
    if src not in ("pdb", "alphafold"):
        raise HTTPException(status_code=400, detail="source must be 'pdb' or 'alphafold'")

    if src == "pdb":
        if not payload.pdb_id:
            raise HTTPException(status_code=400, detail="pdb_id required for source=pdb")
        content = fetch_from_pdb(payload.pdb_id)
        pdb_id = payload.pdb_id.strip()
        af_id = None
    else:
        if not payload.uniprot_id:
            raise HTTPException(status_code=400, detail="uniprot_id required for source=alphafold")
        content = fetch_from_alphafold(payload.uniprot_id)
        pdb_id = None
        af_id = payload.uniprot_id.strip()

    meta = parse_pdb_metadata(content)
    sequences = meta.get("sequences") if isinstance(meta, dict) else {}
    chain_ids = meta.get("chain_ids") if isinstance(meta, dict) else None
    resolution = meta.get("resolution") if isinstance(meta, dict) else None
    method = meta.get("method") if isinstance(meta, dict) else None

    if payload.chain_id:
        cid = payload.chain_id.strip()
        if cid:
            chain_ids = [cid]
            if isinstance(sequences, dict) and cid in sequences:
                sequences = {cid: sequences[cid]}

    with db_session() as db:
        s = ProteinStructure(
            feature_id=payload.feature_id,
            pdb_id=pdb_id,
            alphafold_uniprot_id=af_id,
            source=src,
            chain_ids=chain_ids if isinstance(chain_ids, list) and chain_ids else ([payload.chain_id] if payload.chain_id else None),
            resolution=float(resolution) if isinstance(resolution, (int, float)) else None,
            method=str(method) if method else None,
            prep_status="raw",
            prep_log=None,
        )
        db.add(s)
        db.commit()
        db.refresh(s)

        path = save_structure_file(s.id, content, "pdb")
        p = Path(path)
        f = StructureFile(
            structure_id=s.id,
            file_type="pdb",
            file_path=path,
            file_size_bytes=p.stat().st_size if p.exists() else None,
            md5_hash=_md5(content),
        )
        db.add(f)

        # Save FASTA from SEQRES if present.
        if isinstance(sequences, dict) and sequences:
            fasta = to_fasta(sequences, header_prefix=str(s.id))
            fasta_bytes = fasta.encode("utf-8")
            fasta_path = save_structure_file(s.id, fasta_bytes, "fasta")
            fp = Path(fasta_path)
            db.add(
                StructureFile(
                    structure_id=s.id,
                    file_type="fasta",
                    file_path=fasta_path,
                    file_size_bytes=fp.stat().st_size if fp.exists() else None,
                    md5_hash=_md5(fasta_bytes),
                )
            )

        db.commit()
        db.refresh(s)
        _ = list(s.files or [])
        return ProteinStructureResponse.model_validate(s)


if HAS_MULTIPART:

    @router.post("/upload", response_model=ProteinStructureResponse)
    async def upload_structure(
        file: UploadFile = File(...),
        feature_id: Optional[str] = Form(None),
        chain_id: Optional[str] = Form(None),
        file_type: Optional[str] = Form("pdb"),
    ) -> ProteinStructureResponse:
        """Upload a structure file and register it."""
        try:
            content = await file.read()
        except Exception as e:  # noqa: BLE001
            raise HTTPException(status_code=400, detail=f"Failed to read upload: {e}")

        ft = (file_type or "pdb").strip().lower()
        fid: Optional[UUID] = None
        if feature_id:
            try:
                fid = UUID(feature_id)
            except Exception:
                raise HTTPException(status_code=400, detail="feature_id must be a UUID")

        # Parse PDB metadata if file looks like PDB.
        resolution = None
        method = None
        chain_ids = [chain_id] if chain_id else None
        sequences: dict = {}
        if ft in ("pdb", "prepared"):
            meta = parse_pdb_metadata(content)
            if isinstance(meta, dict):
                sequences = meta.get("sequences") if isinstance(meta.get("sequences"), dict) else {}
                chain_ids = meta.get("chain_ids") if isinstance(meta.get("chain_ids"), list) else chain_ids
                resolution = meta.get("resolution")
                method = meta.get("method")
            if chain_id and sequences and chain_id in sequences:
                sequences = {chain_id: sequences[chain_id]}
                chain_ids = [chain_id]

        with db_session() as db:
            s = ProteinStructure(
                feature_id=fid,
                source="uploaded",
                chain_ids=chain_ids,
                resolution=float(resolution) if isinstance(resolution, (int, float)) else None,
                method=str(method) if method else None,
                prep_status="raw" if ft == "pdb" else "prepared" if ft == "prepared" else "raw",
            )
            db.add(s)
            db.commit()
            db.refresh(s)

            path = save_structure_file(s.id, content, ft)
            p = Path(path)
            db.add(
                StructureFile(
                    structure_id=s.id,
                    file_type=ft,
                    file_path=path,
                    file_size_bytes=p.stat().st_size if p.exists() else None,
                    md5_hash=_md5(content),
                )
            )

            if sequences:
                fasta = to_fasta(sequences, header_prefix=str(s.id))
                fasta_bytes = fasta.encode("utf-8")
                fasta_path = save_structure_file(s.id, fasta_bytes, "fasta")
                fp = Path(fasta_path)
                db.add(
                    StructureFile(
                        structure_id=s.id,
                        file_type="fasta",
                        file_path=fasta_path,
                        file_size_bytes=fp.stat().st_size if fp.exists() else None,
                        md5_hash=_md5(fasta_bytes),
                    )
                )

            db.commit()
            db.refresh(s)
            _ = list(s.files or [])
            return ProteinStructureResponse.model_validate(s)

else:

    @router.post("/upload")
    async def upload_structure_unavailable() -> dict:
        raise HTTPException(
            status_code=503,
            detail='File uploads require "python-multipart" to be installed on the server.',
        )


@router.post("/{structure_id}/prepare", response_model=ProteinStructureResponse)
def prepare_structure_endpoint(structure_id: UUID, payload: PrepareRequest) -> ProteinStructureResponse:
    with db_session() as db:
        s = db.query(ProteinStructure).filter(ProteinStructure.id == structure_id).first()
        if not s:
            raise HTTPException(status_code=404, detail="Structure not found")

        raw = next((f for f in (s.files or []) if f.file_type == "pdb"), None)
        if raw is None:
            raise HTTPException(status_code=400, detail="No raw PDB file found for structure")

        in_path = raw.file_path
        out_path = str(Path(in_path).with_name("prepared.pdb"))
        try:
            prepare_structure(in_path, out_path, chain_id=payload.chain_id)
            s.prep_status = "prepared"
            s.prep_log = "prepared ok"
        except Exception as e:  # noqa: BLE001
            s.prep_status = "failed"
            s.prep_log = str(e)
            db.commit()
            raise HTTPException(status_code=500, detail=f"Preparation failed: {e}")

        out_bytes = Path(out_path).read_bytes()
        f = StructureFile(
            structure_id=s.id,
            file_type="prepared",
            file_path=out_path,
            file_size_bytes=len(out_bytes),
            md5_hash=_md5(out_bytes),
        )
        db.add(f)
        db.commit()
        db.refresh(s)
        _ = list(s.files or [])
        return ProteinStructureResponse.model_validate(s)


@router.get("/{structure_id}/download")
def download_structure(structure_id: UUID, file_type: str = Query("prepared")):
    ft = (file_type or "").strip().lower()
    with db_session() as db:
        s = db.query(ProteinStructure).filter(ProteinStructure.id == structure_id).first()
        if not s:
            raise HTTPException(status_code=404, detail="Structure not found")
        f = next((x for x in (s.files or []) if x.file_type.lower() == ft), None)
        if f is None:
            raise HTTPException(status_code=404, detail="File not found for requested file_type")
        path = f.file_path
        if not path or not Path(path).exists():
            raise HTTPException(status_code=404, detail="File missing on disk")
        filename = Path(path).name
        return FileResponse(path, filename=filename, media_type="chemical/x-pdb")


__all__ = ["router"]


