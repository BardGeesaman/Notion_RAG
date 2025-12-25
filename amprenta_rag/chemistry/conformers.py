"""RDKit-based conformer generation utilities (best-effort, optional dependency)."""

from __future__ import annotations

from contextlib import contextmanager
from typing import List, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS, rdMolAlign

    RDKIT_AVAILABLE = True
except Exception:  # noqa: BLE001
    RDKIT_AVAILABLE = False


def _require_rdkit() -> None:
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for conformer generation. Install rdkit to enable 3D features.")


@contextmanager
def _timeout(seconds: int):
    """Best-effort timeout guard (POSIX only)."""
    try:
        import signal
    except Exception:
        # No timeout available in this environment; proceed without guarding.
        yield
        return

    def _handler(signum, frame):  # noqa: ANN001
        raise TimeoutError("Conformer embedding timed out")

    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(int(seconds))
    try:
        yield
    finally:
        try:
            signal.alarm(0)
        except Exception:
            pass
        try:
            signal.signal(signal.SIGALRM, old)
        except Exception:
            pass


def _mol_from_smiles(smiles: str):
    _require_rdkit()
    m = Chem.MolFromSmiles(str(smiles or ""))
    if m is None:
        raise ValueError("Invalid SMILES")
    return m


def _copy_single_conformer(mol, conf_id: int):
    _require_rdkit()
    m2 = Chem.Mol(mol)
    conf = mol.GetConformer(int(conf_id))
    conf2 = Chem.Conformer(conf)
    m2.RemoveAllConformers()
    m2.AddConformer(conf2, assignId=True)
    return m2


def _mmff_energy(mol, conf_id: int) -> Optional[float]:
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if props is None:
        return None
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=int(conf_id))
    if ff is None:
        return None
    return float(ff.CalcEnergy())


def _uff_energy(mol, conf_id: int) -> Optional[float]:
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=int(conf_id))
    if ff is None:
        return None
    return float(ff.CalcEnergy())


def _optimize_mmff(mol, conf_id: int, max_iters: int = 200) -> Optional[float]:
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if props is None:
        return None
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=int(conf_id))
    if ff is None:
        return None
    ff.Minimize(maxIts=int(max_iters))
    return float(ff.CalcEnergy())


def _optimize_uff(mol, conf_id: int, max_iters: int = 200) -> Optional[float]:
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=int(conf_id))
    if ff is None:
        return None
    ff.Minimize(maxIts=int(max_iters))
    return float(ff.CalcEnergy())


def generate_conformers(smiles: str, n_conformers: int = 5, method: str = "ETKDG") -> List["Chem.Mol"]:
    """
    Generate multiple conformers for a SMILES (best-effort).

    Returns a list of RDKit Mol objects, each containing exactly one conformer.
    """
    _require_rdkit()
    n = int(n_conformers)
    if n < 1:
        raise ValueError("n_conformers must be >= 1")

    m = _mol_from_smiles(smiles)
    m = Chem.AddHs(m)

    meth = (method or "ETKDG").upper()
    if meth in ("ETKDG", "ETKDG3", "ETKDGv3".upper()):
        params = AllChem.ETKDGv3()
    else:
        params = AllChem.ETKDG()

    params.randomSeed = 0xF00D
    params.numThreads = 0
    # RDKit version differences: when using `params=...`, `maxAttempts` kwarg may be unsupported.
    try:
        params.maxAttempts = 50  # type: ignore[attr-defined]
    except Exception:
        pass

    # ETKDG can hang on pathological inputs; apply best-effort timeout guard.
    with _timeout(15):
        conf_ids = list(AllChem.EmbedMultipleConfs(m, numConfs=n, params=params))

    if not conf_ids:
        raise ValueError("Failed to embed conformers")

    return [_copy_single_conformer(m, cid) for cid in conf_ids]


def optimize_conformer(mol: "Chem.Mol", force_field: str = "MMFF") -> "Chem.Mol":
    """Optimize conformer geometry in-place; MMFF preferred with UFF fallback."""
    _require_rdkit()
    if mol is None or mol.GetNumConformers() < 1:
        raise ValueError("Mol has no conformers")

    ff = (force_field or "MMFF").upper()
    for cid in range(mol.GetNumConformers()):
        if ff == "MMFF":
            e = _optimize_mmff(mol, cid)
            if e is None:
                _ = _optimize_uff(mol, cid)
        else:
            e = _optimize_uff(mol, cid)
            if e is None:
                _ = _optimize_mmff(mol, cid)
    return mol


def get_lowest_energy_conformer(mol: "Chem.Mol") -> "Chem.Mol":
    """Select and return the single-conformer Mol with the minimum (MMFF/UFF) energy."""
    _require_rdkit()
    if mol is None or mol.GetNumConformers() < 1:
        raise ValueError("Mol has no conformers")

    best_idx = 0
    best_e: float | None = None
    for cid in range(mol.GetNumConformers()):
        e = _mmff_energy(mol, cid)
        if e is None:
            e = _uff_energy(mol, cid)
        if e is None:
            continue
        if best_e is None or e < best_e:
            best_e = e
            best_idx = cid

    return _copy_single_conformer(mol, best_idx)


def conformer_to_pdb(mol: "Chem.Mol", conf_id: int = 0) -> str:
    """Return PDB block string for a given conformer."""
    _require_rdkit()
    if mol is None or mol.GetNumConformers() < 1:
        raise ValueError("Mol has no conformers")
    return str(Chem.MolToPDBBlock(mol, confId=int(conf_id)))


def align_conformers(mol: "Chem.Mol") -> None:
    """Align all conformers in a multi-conformer Mol to the first conformer."""
    _require_rdkit()
    if mol is None or mol.GetNumConformers() < 2:
        return
    try:
        AllChem.AlignMolConformers(mol)
    except Exception:
        return


def conformer_energies(mol: "Chem.Mol", *, prefer: str = "MMFF") -> List[float]:
    """Compute per-conformer energies (MMFF preferred, UFF fallback)."""
    _require_rdkit()
    if mol is None or mol.GetNumConformers() < 1:
        raise ValueError("Mol has no conformers")
    prefer_ff = (prefer or "MMFF").upper()
    out: List[float] = []
    for cid in range(mol.GetNumConformers()):
        e: Optional[float]
        if prefer_ff == "MMFF":
            e = _mmff_energy(mol, cid)
            if e is None:
                e = _uff_energy(mol, cid)
        else:
            e = _uff_energy(mol, cid)
            if e is None:
                e = _mmff_energy(mol, cid)
        out.append(float(e) if e is not None else float("nan"))
    return out


def align_molecules_to_reference(mols: List["Chem.Mol"], reference_idx: int = 0) -> None:
    """
    Best-effort alignment of multiple molecules to a reference (in-place).

    Uses MCS-based atom mapping when possible, with a short timeout.
    """
    _require_rdkit()
    if not mols:
        return
    ref_i = int(reference_idx)
    ref_i = max(0, min(ref_i, len(mols) - 1))
    ref = mols[ref_i]
    if ref is None or ref.GetNumConformers() < 1:
        return

    for i, m in enumerate(mols):
        if i == ref_i or m is None or m.GetNumConformers() < 1:
            continue
        try:
            res = rdFMCS.FindMCS([ref, m], timeout=2, ringMatchesRingOnly=True, completeRingsOnly=True)
            if not res or not res.smartsString:
                continue
            patt = Chem.MolFromSmarts(res.smartsString)
            if patt is None:
                continue
            ref_match = ref.GetSubstructMatch(patt)
            mov_match = m.GetSubstructMatch(patt)
            if not ref_match or not mov_match or len(ref_match) != len(mov_match):
                continue
            amap = list(zip(mov_match, ref_match))
            rdMolAlign.AlignMol(m, ref, atomMap=amap)
        except Exception:
            continue


__all__ = [
    "RDKIT_AVAILABLE",
    "generate_conformers",
    "optimize_conformer",
    "get_lowest_energy_conformer",
    "conformer_to_pdb",
    "align_conformers",
    "conformer_energies",
    "align_molecules_to_reference",
]


