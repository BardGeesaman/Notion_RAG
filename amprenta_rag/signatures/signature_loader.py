"""
Signature loader for multi-omics signatures.

Loads signature definitions from TSV files or Notion with components including:
- Feature name (gene, protein, metabolite, or lipid species)
- Feature type (gene, protein, metabolite, lipid - auto-detected if not provided)
- Direction (↑ / ↓ / neutral / complex)
- Weight (float or None)
"""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from amprenta_rag.signatures.feature_type_inference import infer_feature_type


@dataclass
class SignatureComponent:
    """
    A single component of a multi-omics signature.

    Attributes:
        feature_name: Feature name (gene, protein, metabolite, or lipid species)
        feature_type: Feature type (gene, protein, metabolite, lipid) - auto-inferred if None
        direction: Direction of change (↑, ↓, neutral, complex, or None)
        weight: Optional weight for this component (default: 1.0)
    """

    feature_name: str
    feature_type: Optional[str] = None
    direction: Optional[str] = None
    weight: Optional[float] = None

    def __post_init__(self):
        """Normalize direction, infer feature type, and set defaults."""
        # Infer feature type if not provided
        if not self.feature_type:
            self.feature_type = infer_feature_type(self.feature_name)

        # Normalize direction symbols
        if self.direction:
            direction_lower = self.direction.lower().strip()
            if direction_lower in ["↑", "up", "increased", "increase", "+"]:
                self.direction = "↑"
            elif direction_lower in ["↓", "down", "decreased", "decrease", "-"]:
                self.direction = "↓"
            elif direction_lower in ["neutral", "none", "0", ""]:
                self.direction = "neutral"
            elif direction_lower in ["complex", "mixed"]:
                self.direction = "complex"
            # Keep original if not recognized

        # Default weight to 1.0 if not provided
        if self.weight is None:
            self.weight = 1.0

    @property
    def species(self) -> str:
        """
        Backward compatibility: return feature_name as species.
        
        Allows existing code that uses .species to continue working.
        """
        return self.feature_name


@dataclass
class Signature:
    """
    A complete multi-omics signature definition.

    Attributes:
        name: Signature name/identifier
        components: List of signature components (can include genes, proteins, metabolites, lipids)
        description: Optional description
        modalities: Set of feature types present in this signature (auto-computed)
    """

    name: str
    components: List[SignatureComponent]
    description: Optional[str] = None
    modalities: Optional[List[str]] = None

    def __post_init__(self):
        """Auto-compute modalities from components if not provided."""
        if self.modalities is None:
            self.modalities = list(
                set(
                    comp.feature_type
                    for comp in self.components
                    if comp.feature_type
                )
            )


def load_signature_from_tsv(tsv_path: Path) -> Signature:
    """
    Load a multi-omics signature from a TSV file.

    Expected TSV format (with feature_type column):
        feature_type    feature_name        direction   weight
        gene            TP53                up          1.0
        protein         P04637              down        0.8
        metabolite      Glutamate           up          1.0
        lipid           Cer(d18:1/16:0)     up          1.0

    OR (without feature_type - auto-detected):
        feature_name        direction   weight
        TP53                up          1.0
        Cer(d18:1/16:0)     up          1.0

    Columns:
        - feature_name (required): Feature identifier/name
        - feature_type (optional): gene, protein, metabolite, lipid (auto-detected if missing)
        - direction (optional): ↑, ↓, neutral, complex, up, down, or empty
        - weight (optional): Numeric weight (default: 1.0)
        - Legacy: "species" column is also supported (treated as feature_name)

    Args:
        tsv_path: Path to TSV/CSV file

    Returns:
        Signature object with loaded components

    Raises:
        FileNotFoundError: If TSV file doesn't exist
        ValueError: If TSV format is invalid
    """
    if not tsv_path.is_file():
        raise FileNotFoundError(f"Signature TSV not found: {tsv_path}")

    components: List[SignatureComponent] = []
    signature_name = tsv_path.stem  # Use filename as signature name

    with tsv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        # Try comma if tab doesn't work
        if not reader.fieldnames:
            f.seek(0)
            reader = csv.DictReader(f, delimiter=",")

        if not reader.fieldnames:
            raise ValueError(f"Invalid TSV format: {tsv_path}")

        # Normalize column names (case-insensitive)
        fieldnames_lower = {f.lower(): f for f in reader.fieldnames}

        # Find feature name column (supports multiple names)
        feature_name_col = None
        for key in ["feature_name", "feature", "name", "species", "metabolite", "lipid", "gene", "protein"]:
            if key in fieldnames_lower:
                feature_name_col = fieldnames_lower[key]
                break

        # Find feature type column
        feature_type_col = None
        for key in ["feature_type", "type", "featuretype"]:
            if key in fieldnames_lower:
                feature_type_col = fieldnames_lower[key]
                break

        # Find direction column
        direction_col = None
        for key in ["direction", "dir", "change"]:
            if key in fieldnames_lower:
                direction_col = fieldnames_lower[key]
                break

        # Find weight column
        weight_col = None
        for key in ["weight", "w", "importance"]:
            if key in fieldnames_lower:
                weight_col = fieldnames_lower[key]
                break

        if not feature_name_col:
            raise ValueError(
                f"TSV must have a feature name column (feature_name, name, species, etc.): {tsv_path}"
            )

        for row in reader:
            feature_name = row.get(feature_name_col, "").strip()
            if not feature_name:
                continue  # Skip empty rows

            # Get feature type (if provided)
            feature_type = None
            if feature_type_col:
                feature_type_raw = row.get(feature_type_col, "").strip()
                if feature_type_raw:
                    feature_type = feature_type_raw.lower()
                    # Normalize feature type values
                    if feature_type in ["gene", "g", "genes"]:
                        feature_type = "gene"
                    elif feature_type in ["protein", "p", "proteins", "prot"]:
                        feature_type = "protein"
                    elif feature_type in ["metabolite", "m", "metabolites", "met"]:
                        feature_type = "metabolite"
                    elif feature_type in ["lipid", "l", "lipids"]:
                        feature_type = "lipid"
                    else:
                        # If not recognized, let inference handle it
                        feature_type = None

            # Get direction
            direction = row.get(direction_col, "").strip() if direction_col else None
            if not direction:
                direction = None

            # Get weight
            weight_str = row.get(weight_col, "").strip() if weight_col else None
            weight = None
            if weight_str:
                try:
                    weight = float(weight_str)
                except ValueError:
                    weight = None  # Invalid weight, use default

            components.append(
                SignatureComponent(
                    feature_name=feature_name,
                    feature_type=feature_type,  # Will be auto-inferred if None
                    direction=direction,
                    weight=weight,
                )
            )

    if not components:
        raise ValueError(f"No valid components found in signature: {tsv_path}")

    return Signature(
        name=signature_name,
        components=components,
    )


def load_signatures_from_postgres() -> List[Signature]:
    """
    Load all signatures from Postgres database.
    
    Returns:
        List of Signature objects
    """
    from amprenta_rag.database.base import get_db
    from amprenta_rag.database.models import Signature as SignatureModel
    
    db = next(get_db())
    try:
        db_signatures = db.query(SignatureModel).all()
        signatures = []
        
        for db_sig in db_signatures:
            # Convert components if stored
            components = []
            if hasattr(db_sig, 'components') and db_sig.components:
                for comp_data in db_sig.components:
                    if isinstance(comp_data, dict):
                        components.append(SignatureComponent(
                            feature_name=comp_data.get('feature_name', ''),
                            feature_type=comp_data.get('feature_type'),
                            direction=comp_data.get('direction'),
                            weight=comp_data.get('weight'),
                        ))
            
            sig = Signature(
                name=db_sig.name or str(db_sig.id),
                components=components,
                description=db_sig.description if hasattr(db_sig, 'description') else None,
            )
            # Store the database ID for reference
            sig.id = db_sig.id
            signatures.append(sig)
        
        return signatures
    finally:
        db.close()
