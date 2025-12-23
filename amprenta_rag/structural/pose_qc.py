"""Pose QC utilities: clashes, ligand efficiency, PLIP analysis."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from amprenta_rag.database.models import DockingPose
from amprenta_rag.structural.plip_parser import PlipResult, parse_plip_xml
from amprenta_rag.structural.plip_runner import run_plip


@dataclass(frozen=True)
class PoseQualityMetrics:
    num_hbonds: int
    num_hydrophobic: int
    num_salt_bridges: int
    num_pi_stacking: int
    num_pi_cation: int
    num_halogen: int
    num_metal: int
    total_interactions: int
    has_clashes: bool
    ligand_efficiency: Optional[float]


@dataclass(frozen=True)
class PoseInteractionRecord:
    interaction_type: str
    ligand_atom: Optional[str]
    protein_residue: Optional[str]
    distance: Optional[float]
    angle: Optional[float]


def calc_ligand_efficiency(affinity: float, smiles: str) -> Optional[float]:
    """Ligand efficiency = affinity / heavy atom count (kcal/mol per heavy atom).

    Affinity from Vina is typically negative; we return a negative LE value.
    """
    try:
        from rdkit import Chem
    except Exception:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    ha = mol.GetNumHeavyAtoms()
    if ha <= 0:
        return None
    try:
        return float(affinity) / float(ha)
    except Exception:
        return None


def _parse_coords_from_pdb(pdb_path: Path) -> tuple[List[tuple[float, float, float]], List[tuple[float, float, float]]]:
    """Return (protein_coords, ligand_coords) using ATOM vs HETATM split."""
    prot: List[tuple[float, float, float]] = []
    lig: List[tuple[float, float, float]] = []
    for ln in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        # PDB fixed columns: x=30-38, y=38-46, z=46-54 (1-indexed)
        try:
            x = float(ln[30:38].strip())
            y = float(ln[38:46].strip())
            z = float(ln[46:54].strip())
        except Exception:
            parts = ln.split()
            floats: List[float] = []
            for p in parts:
                try:
                    floats.append(float(p))
                except Exception:
                    continue
            if len(floats) >= 3:
                x, y, z = floats[-3], floats[-2], floats[-1]
            else:
                continue
        if ln.startswith("HETATM"):
            lig.append((x, y, z))
        else:
            prot.append((x, y, z))
    return prot, lig


def check_clashes(complex_pdb: str, threshold: float = 2.0) -> bool:
    """Return True if any ligand atom is within threshold Å of any protein atom."""
    p = Path(complex_pdb)
    if not p.exists():
        raise FileNotFoundError(str(p))
    prot, lig = _parse_coords_from_pdb(p)
    if not prot or not lig:
        return False

    thr2 = float(threshold) ** 2
    for lx, ly, lz in lig:
        for px, py, pz in prot:
            dx = lx - px
            dy = ly - py
            dz = lz - pz
            if (dx * dx + dy * dy + dz * dz) < thr2:
                return True
    return False


def _pdbqt_to_pdb(pdbqt_path: Path, pdb_path: Path) -> bool:
    if shutil.which("obabel") is None:
        return False
    cmd = ["obabel", "-i", "pdbqt", str(pdbqt_path), "-o", "pdb", "-O", str(pdb_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=str(pdb_path.parent))
    return proc.returncode == 0 and pdb_path.exists()


def _find_plip_xml(output_dir: Path) -> Optional[Path]:
    xmls = list(output_dir.rglob("*.xml"))
    if not xmls:
        return None
    # Prefer report.xml if present
    for x in xmls:
        if x.name.lower() == "report.xml":
            return x
    return xmls[0]


def analyze_pose(pose: DockingPose, receptor_pdb: str) -> tuple[PoseQualityMetrics, List[PoseInteractionRecord]]:
    """Analyze a docked pose against a receptor PDB using PLIP + basic QC."""
    if not pose.pose_pdbqt_path:
        raise ValueError("pose.pose_pdbqt_path is missing")

    receptor_path = Path(receptor_pdb)
    if not receptor_path.exists():
        raise FileNotFoundError(str(receptor_path))

    pose_pdbqt = Path(pose.pose_pdbqt_path)
    if not pose_pdbqt.exists():
        raise FileNotFoundError(str(pose_pdbqt))

    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        ligand_pdb = td_path / "ligand.pdb"
        if not _pdbqt_to_pdb(pose_pdbqt, ligand_pdb):
            raise RuntimeError("Failed to convert pose PDBQT to PDB (need obabel)")

        complex_pdb = td_path / "complex.pdb"
        complex_pdb.write_text(
            receptor_path.read_text(encoding="utf-8", errors="replace")
            + "\n"
            + ligand_pdb.read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8",
        )

        plip_out = td_path / "plip"
        ok = run_plip(str(complex_pdb), str(plip_out))
        if not ok:
            raise RuntimeError("PLIP execution failed")

        xml_path = _find_plip_xml(plip_out)
        if not xml_path:
            raise RuntimeError("PLIP output XML not found")

        plip: PlipResult = parse_plip_xml(str(xml_path))
        has_clashes = check_clashes(str(complex_pdb), threshold=2.0)

        le: Optional[float] = None
        try:
            affinity = float(pose.binding_affinity) if pose.binding_affinity is not None else None
        except Exception:
            affinity = None

        smiles = getattr(getattr(pose, "compound", None), "smiles", None)
        if affinity is not None and isinstance(smiles, str) and smiles:
            le = calc_ligand_efficiency(affinity, smiles)

        qm = PoseQualityMetrics(
            num_hbonds=plip.num_hbonds,
            num_hydrophobic=plip.num_hydrophobic,
            num_salt_bridges=plip.num_salt_bridges,
            num_pi_stacking=plip.num_pi_stacking,
            num_pi_cation=plip.num_pi_cation,
            num_halogen=plip.num_halogen,
            num_metal=plip.num_metal,
            total_interactions=plip.total_interactions,
            has_clashes=bool(has_clashes),
            ligand_efficiency=le,
        )

        interactions = [
            PoseInteractionRecord(
                interaction_type=i.interaction_type,
                ligand_atom=i.ligand_atom,
                protein_residue=i.protein_residue,
                distance=i.distance,
                angle=i.angle,
            )
            for i in plip.interactions
        ]

        return qm, interactions


__all__ = [
    "PoseQualityMetrics",
    "PoseInteractionRecord",
    "calc_ligand_efficiency",
    "check_clashes",
    "analyze_pose",
]


