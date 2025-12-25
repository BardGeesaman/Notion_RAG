"""Structural alert filters (PAINS / Brenk / Lilly) using RDKit.

These filters are intended for fast screening of compounds for common assay
interference motifs and chemical liabilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple


try:
    from rdkit import Chem  # type: ignore
    from rdkit.Chem import FilterCatalog  # type: ignore
    from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore

    RDKIT_AVAILABLE = True
except Exception:  # noqa: BLE001
    Chem = None  # type: ignore
    FilterCatalog = None  # type: ignore
    FilterCatalogParams = None  # type: ignore
    RDKIT_AVAILABLE = False


@dataclass
class AlertResult:
    alert_type: str  # "PAINS", "LILLY", "BRENK"
    pattern_name: str
    description: str
    severity: str  # "high" | "medium" | "low"
    matched_smarts: str

    def to_dict(self) -> dict:
        return {
            "alert_type": self.alert_type,
            "pattern_name": self.pattern_name,
            "description": self.description,
            "severity": self.severity,
            "matched_smarts": self.matched_smarts,
        }


def _mol_from_smiles(smiles: str):
    if not RDKIT_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)  # type: ignore[union-attr]
    return mol


class PAINSFilter:
    """RDKit FilterCatalog-based PAINS filters (A/B/C)."""

    def __init__(self) -> None:
        self.catalog = None
        self.pattern_count = 0
        if not RDKIT_AVAILABLE:
            return
        params = FilterCatalogParams()  # type: ignore[operator]
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)  # type: ignore[union-attr]
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)  # type: ignore[union-attr]
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)  # type: ignore[union-attr]
        self.catalog = FilterCatalog.FilterCatalog(params)  # type: ignore[union-attr]
        try:
            self.pattern_count = int(self.catalog.GetNumEntries())  # type: ignore[union-attr]
        except Exception:  # noqa: BLE001
            self.pattern_count = 0

    def check(self, smiles: str) -> List[AlertResult]:
        if not RDKIT_AVAILABLE or self.catalog is None:
            return []
        mol = _mol_from_smiles(smiles)
        if mol is None:
            return []
        matches = self.catalog.GetMatches(mol)  # type: ignore[union-attr]
        out: List[AlertResult] = []
        for m in matches:
            try:
                desc = str(m.GetDescription())
            except Exception:  # noqa: BLE001
                desc = "PAINS_match"
            # RDKit FilterCatalog does not expose SMARTS directly on the entry (binary serialized query).
            out.append(
                AlertResult(
                    alert_type="PAINS",
                    pattern_name=desc,
                    description=f"PAINS substructure match: {desc}",
                    severity="high",
                    matched_smarts="",
                )
            )
        return out


class BRENKFilter:
    """RDKit FilterCatalog-based Brenk unwanted substructures."""

    def __init__(self) -> None:
        self.catalog = None
        self.pattern_count = 0
        if not RDKIT_AVAILABLE:
            return
        params = FilterCatalogParams()  # type: ignore[operator]
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)  # type: ignore[union-attr]
        self.catalog = FilterCatalog.FilterCatalog(params)  # type: ignore[union-attr]
        try:
            self.pattern_count = int(self.catalog.GetNumEntries())  # type: ignore[union-attr]
        except Exception:  # noqa: BLE001
            self.pattern_count = 0

    def check(self, smiles: str) -> List[AlertResult]:
        if not RDKIT_AVAILABLE or self.catalog is None:
            return []
        mol = _mol_from_smiles(smiles)
        if mol is None:
            return []
        matches = self.catalog.GetMatches(mol)  # type: ignore[union-attr]
        out: List[AlertResult] = []
        for m in matches:
            try:
                desc = str(m.GetDescription())
            except Exception:  # noqa: BLE001
                desc = "BRENK_match"
            out.append(
                AlertResult(
                    alert_type="BRENK",
                    pattern_name=desc,
                    description=f"Brenk unwanted substructure match: {desc}",
                    severity="medium",
                    matched_smarts="",
                )
            )
        return out


# Lilly-style structural alerts (MVP subset). Name -> (SMARTS, description)
LILLY_PATTERNS: Dict[str, Tuple[str, str]] = {
    # Reactive
    "acyl_halide": ("[CX3](=O)[F,Cl,Br,I]", "Acyl halide (highly reactive electrophile)"),
    "michael_acceptor": ("[CX3]=[CX3][CX3](=O)", "Michael acceptor (alpha,beta-unsaturated carbonyl)"),
    "aldehyde": ("[CX3H1](=O)[#6]", "Aldehyde (reactive carbonyl)"),
    "isocyanate": ("N=C=O", "Isocyanate (reactive)"),
    "isothiocyanate": ("N=C=S", "Isothiocyanate (reactive)"),
    "aziridine": ("C1NC1", "Aziridine (strained, reactive)"),
    # Toxicophores
    "nitro_aromatic": ("[NX3](=O)=O[a,c]", "Nitro group on aromatic system (toxicophore risk)"),
    "aniline": ("Nc1ccccc1", "Aniline / aromatic amine (toxicity risk)"),
    "hydrazine": ("NN", "Hydrazine motif (toxicity risk)"),
    "alkyl_halide": ("[CX4][Cl,Br,I]", "Alkyl halide (alkylating potential)"),
    "nitroso": ("N=O", "Nitroso group (toxicity risk)"),
    # Metabolic / instability
    "epoxide": ("C1OC1", "Epoxide (reactive metabolite risk)"),
    "peroxide": ("OO", "Peroxide (unstable/reactive)"),
    "thiocarbonyl": ("[CX3](=S)[#6,#7,#8,#16]", "Thiocarbonyl (reactive/metabolic liability)"),
    "quinone": ("O=c1ccc(=O)cc1", "Quinone (redox cycling potential)"),
    "catechol": ("c1ccc(O)c(O)c1", "Catechol (redox/covalent liability potential)"),
    "phenol": ("c1ccc(O)cc1", "Phenol (possible metabolic liability)"),
    "sulfonate_ester": ("S(=O)(=O)O[CX4]", "Sulfonate ester (reactive electrophile)"),
    "azide": ("N=[N+]=[N-]", "Azide (reactive/toxicity risk)"),
    "thiol": ("[SX2H]", "Thiol (reactivity/oxidation risk)"),
}


class LillyFilter:
    """SMARTS-based Lilly alert subset (MVP)."""

    def __init__(self) -> None:
        self._queries: List[Tuple[str, str, str]] = []
        self.pattern_count = len(LILLY_PATTERNS)
        if not RDKIT_AVAILABLE:
            return
        for name, (smarts, desc) in LILLY_PATTERNS.items():
            q = Chem.MolFromSmarts(smarts)  # type: ignore[union-attr]
            if q is None:
                continue
            self._queries.append((name, smarts, desc))

    def check(self, smiles: str) -> List[AlertResult]:
        if not RDKIT_AVAILABLE:
            return []
        mol = _mol_from_smiles(smiles)
        if mol is None:
            return []

        out: List[AlertResult] = []
        for name, smarts, desc in self._queries:
            try:
                q = Chem.MolFromSmarts(smarts)  # type: ignore[union-attr]
                if q is None:
                    continue
                if mol.HasSubstructMatch(q):
                    out.append(
                        AlertResult(
                            alert_type="LILLY",
                            pattern_name=name,
                            description=desc,
                            severity="medium",
                            matched_smarts=smarts,
                        )
                    )
            except Exception:  # noqa: BLE001
                continue
        return out


__all__ = [
    "AlertResult",
    "PAINSFilter",
    "BRENKFilter",
    "LillyFilter",
    "LILLY_PATTERNS",
    "RDKIT_AVAILABLE",
]


