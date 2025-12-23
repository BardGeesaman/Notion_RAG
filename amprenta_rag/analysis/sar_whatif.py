"""SAR What-If utilities (SMILES validation, 2D rendering, property comparison, scaffold hops)."""

from __future__ import annotations

import base64
import io
from typing import Dict, List

import pandas as pd


def _require_rdkit():
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Draw  # type: ignore
        from rdkit.Chem import rdChemReactions  # type: ignore
    except ImportError as e:
        raise ImportError("RDKit is required for SAR What-If utilities.") from e
    return Chem, Draw, rdChemReactions


def validate_smiles(smiles: str) -> bool:
    """Return True if SMILES parses with RDKit."""
    Chem, _, _ = _require_rdkit()
    if not smiles or not str(smiles).strip():
        return False
    mol = Chem.MolFromSmiles(str(smiles).strip())
    return mol is not None


def render_structure_2d(smiles: str) -> bytes:
    """Render a SMILES to a 2D PNG image (bytes)."""
    Chem, Draw, _ = _require_rdkit()
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        raise ValueError("Invalid SMILES")
    img = Draw.MolToImage(mol, size=(320, 240))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def _img_b64(smiles: str) -> str:
    png = render_structure_2d(smiles)
    return base64.b64encode(png).decode("utf-8")


def compare_properties(smiles_list: List[str]) -> pd.DataFrame:
    """Compare ADMET properties for a list of SMILES.

    Returns DataFrame columns:
      - smiles
      - structure_img (base64 PNG)
      - logp, logs, herg (flattened values if available)
    """
    from amprenta_rag.ml.admet import get_admet_predictor

    cleaned = [str(s).strip() for s in (smiles_list or []) if s and str(s).strip()]
    if not cleaned:
        return pd.DataFrame(columns=["smiles", "structure_img"])

    predictor = get_admet_predictor()
    preds = predictor.predict(cleaned, endpoints=["logp", "logs", "herg"], include_shap=False)

    rows: List[Dict[str, object]] = []
    for p in preds:
        s = str(p.get("smiles") or "").strip()
        row: Dict[str, object] = {"smiles": s}
        try:
            row["structure_img"] = _img_b64(s)
        except Exception:
            row["structure_img"] = None

        # Flatten common endpoint outputs
        for ep in ("logp", "logs", "herg"):
            val = p.get(ep)
            if isinstance(val, dict):
                if "value" in val:
                    row[ep] = val.get("value")
                elif "probability" in val:
                    row[ep] = val.get("probability")
                else:
                    row[ep] = None
            else:
                row[ep] = val

        if p.get("error"):
            row["error"] = p.get("error")
        rows.append(row)

    return pd.DataFrame(rows)


TRANSFORMATIONS: Dict[str, Dict[str, str]] = {
    "benzene_to_pyridine": {
        "label": "Benzene → Pyridine",
        "description": "Replace a phenyl ring with a pyridine ring (one C→N).",
        "exact_smiles": "c1ccccc1",
        "query_smarts": "c1ccccc1",
        "replacement_smiles": "n1ccccc1",
    },
    "cyclohexane_to_piperidine": {
        "label": "Cyclohexane → Piperidine",
        "description": "Replace a cyclohexane ring with piperidine (one C→N).",
        "exact_smiles": "C1CCCCC1",
        "query_smarts": "C1CCCCC1",
        "replacement_smiles": "N1CCCCC1",
    },
}


def scaffold_hop(smiles: str, transformation: str) -> List[str]:
    """Apply a predefined scaffold hop transformation to a SMILES.

    Returns a list of product SMILES (unique).
    """
    Chem, _, _ = _require_rdkit()
    key = (transformation or "").strip()
    if key not in TRANSFORMATIONS:
        raise ValueError("Unknown transformation")

    q = Chem.MolFromSmarts(TRANSFORMATIONS[key]["query_smarts"])
    if q is None:
        raise ValueError("Invalid transformation definition")
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        raise ValueError("Invalid SMILES")

    # MVP reliability: exact-match shortcut for the canonical scaffold definitions.
    exact = TRANSFORMATIONS[key].get("exact_smiles")
    if exact:
        try:
            in_can = Chem.MolToSmiles(mol, canonical=True)
            ex_mol = Chem.MolFromSmiles(str(exact))
            if ex_mol is not None and in_can == Chem.MolToSmiles(ex_mol, canonical=True):
                rep = Chem.MolFromSmiles(TRANSFORMATIONS[key]["replacement_smiles"])
                if rep is not None:
                    return [Chem.MolToSmiles(rep, canonical=True)]
        except Exception:
            pass

    matches = mol.GetSubstructMatches(q)
    if not matches:
        return []

    out: List[str] = []
    for match in matches:
        # MVP: mutate the first atom in the matched substructure (C -> N)
        idx = int(match[0])
        rw = Chem.RWMol(mol)
        atom = rw.GetAtomWithIdx(idx)
        atom.SetAtomicNum(7)
        # Preserve aromaticity if the atom is aromatic in the source
        if atom.GetIsAromatic():
            atom.SetIsAromatic(True)
        try:
            new = rw.GetMol()
            Chem.SanitizeMol(new)
            s = Chem.MolToSmiles(new, canonical=True)
            if s and s not in out:
                out.append(s)
        except Exception:
            continue

    return out


__all__ = [
    "validate_smiles",
    "render_structure_2d",
    "compare_properties",
    "scaffold_hop",
    "TRANSFORMATIONS",
]


