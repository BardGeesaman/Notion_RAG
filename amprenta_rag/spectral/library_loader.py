"""Spectral library loader (MGF -> DB)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from amprenta_rag.database.models import SpectralLibrary, SpectralReference
from amprenta_rag.database.session import db_session
from amprenta_rag.spectral.library_parser import parse_mgf


def _infer_lipid_class(name: str) -> Optional[str]:
    # Simple heuristic: take prefix before first space or "("
    n = (name or "").strip()
    if not n:
        return None
    for sep in (" ", "(", "["):
        if sep in n:
            return n.split(sep, 1)[0]
    return n


def load_library(mgf_path: str, library_name: str, version: str | None = None, source_url: str | None = None) -> SpectralLibrary:
    """Create SpectralLibrary and bulk insert SpectralReference records."""
    p = Path(mgf_path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    spectra = parse_mgf(str(p))
    with db_session() as db:
        # Upsert library by (name, version)
        existing = db.query(SpectralLibrary).filter_by(name=library_name, version=version).first()
        if existing is None:
            lib = SpectralLibrary(
                name=library_name,
                version=version,
                source_url=source_url,
                file_path=str(p),
                n_spectra=len(spectra),
            )
            db.add(lib)
            db.commit()
            db.refresh(lib)
        else:
            existing.source_url = source_url or existing.source_url
            existing.file_path = str(p)
            existing.n_spectra = len(spectra)
            db.add(existing)
            db.commit()
            db.refresh(existing)
            lib = existing

        refs = []
        for s in spectra:
            mzs = [float(x[0]) for x in (s.get("peaks") or [])]
            intens = [float(x[1]) for x in (s.get("peaks") or [])]
            name = str(s.get("name") or "unknown")
            refs.append(
                SpectralReference(
                    library_id=lib.id,
                    lipid_name=name,
                    lipid_class=_infer_lipid_class(name),
                    smiles=None,
                    inchi_key=None,
                    precursor_mz=float(s.get("precursor_mz")),
                    precursor_type=s.get("precursor_type"),
                    collision_energy=s.get("collision_energy"),
                    spectrum={"mz": mzs, "intensity": intens},
                )
            )

        if refs:
            db.bulk_save_objects(refs)

        return lib


__all__ = ["load_library"]


