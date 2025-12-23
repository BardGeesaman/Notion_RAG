"""Parse PLIP XML reports into interaction records."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
from xml.etree import ElementTree as ET


@dataclass(frozen=True)
class PlipInteraction:
    interaction_type: str
    ligand_atom: Optional[str] = None
    protein_residue: Optional[str] = None
    distance: Optional[float] = None
    angle: Optional[float] = None


@dataclass(frozen=True)
class PlipResult:
    interactions: List[PlipInteraction]

    num_hbonds: int
    num_hydrophobic: int
    num_salt_bridges: int
    num_pi_stacking: int
    num_pi_cation: int
    num_halogen: int
    num_metal: int
    total_interactions: int


_TYPE_TO_TAGS = {
    "hbond": ("hydrogen_bonds", "hydrogen_bond"),
    "hydrophobic": ("hydrophobic_interactions", "hydrophobic_interaction"),
    "saltbridge": ("salt_bridges", "salt_bridge"),
    "pistacking": ("pi_stacks", "pi_stack"),
    "pication": ("pi_cation_interactions", "pi_cation_interaction"),
    "halogen": ("halogen_bonds", "halogen_bond"),
    "metal": ("metal_complexes", "metal_complex"),
}


def _first_text(node: ET.Element, *paths: str) -> Optional[str]:
    for p in paths:
        child = node.find(p)
        if child is not None and child.text:
            t = child.text.strip()
            if t:
                return t
    return None


def _first_float(node: ET.Element, *paths: str) -> Optional[float]:
    t = _first_text(node, *paths)
    if not t:
        return None
    try:
        return float(t)
    except Exception:
        return None


def parse_plip_xml(xml_path: str) -> PlipResult:
    """Parse a PLIP XML file produced by `plip -x`."""
    p = Path(xml_path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    root = ET.fromstring(p.read_text(encoding="utf-8", errors="replace"))

    interactions: List[PlipInteraction] = []
    counts = {k: 0 for k in _TYPE_TO_TAGS.keys()}

    # PLIP structure can have multiple bindingsites; search globally for interaction blocks.
    for itype, (container_tag, item_tag) in _TYPE_TO_TAGS.items():
        for container in root.findall(f".//{container_tag}"):
            for item in container.findall(f"./{item_tag}"):
                ligand_atom = _first_text(item, "ligatom", "ligatom_idx", "lig_atom", "ligand_atom")
                # Build residue identifier if available: CHAIN:RES:NUM
                reschain = _first_text(item, "reschain")
                restype = _first_text(item, "restype")
                resnr = _first_text(item, "resnr")
                protein_residue = None
                if reschain or restype or resnr:
                    protein_residue = ":".join([x for x in [reschain, restype, resnr] if x])
                if not protein_residue:
                    protein_residue = _first_text(item, "residue", "protres")

                dist = _first_float(item, "dist", "distance")
                ang = _first_float(item, "angle", "ang")
                interactions.append(
                    PlipInteraction(
                        interaction_type=itype,
                        ligand_atom=ligand_atom,
                        protein_residue=protein_residue,
                        distance=dist,
                        angle=ang,
                    )
                )
                counts[itype] += 1

    # Map into output fields
    num_hbonds = counts["hbond"]
    num_hydrophobic = counts["hydrophobic"]
    num_salt_bridges = counts["saltbridge"]
    num_pi_stacking = counts["pistacking"]
    num_pi_cation = counts["pication"]
    num_halogen = counts["halogen"]
    num_metal = counts["metal"]
    total = sum(counts.values())

    return PlipResult(
        interactions=interactions,
        num_hbonds=num_hbonds,
        num_hydrophobic=num_hydrophobic,
        num_salt_bridges=num_salt_bridges,
        num_pi_stacking=num_pi_stacking,
        num_pi_cation=num_pi_cation,
        num_halogen=num_halogen,
        num_metal=num_metal,
        total_interactions=total,
    )


__all__ = ["PlipInteraction", "PlipResult", "parse_plip_xml"]


