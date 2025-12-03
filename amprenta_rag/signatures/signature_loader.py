"""
Signature loader for lipid signatures.

Loads signature definitions from TSV files or Notion with components including:
- Species name
- Direction (↑ / ↓ / neutral / complex)
- Weight (float or None)
"""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass
class SignatureComponent:
    """
    A single component of a lipid signature.
    
    Attributes:
        species: Lipid species name
        direction: Direction of change (↑, ↓, neutral, complex, or None)
        weight: Optional weight for this component (default: 1.0)
    """
    species: str
    direction: Optional[str] = None
    weight: Optional[float] = None
    
    def __post_init__(self):
        """Normalize direction and set default weight."""
        if self.direction:
            # Normalize direction symbols
            direction_lower = self.direction.lower().strip()
            if direction_lower in ['↑', 'up', 'increased', 'increase', '+']:
                self.direction = '↑'
            elif direction_lower in ['↓', 'down', 'decreased', 'decrease', '-']:
                self.direction = '↓'
            elif direction_lower in ['neutral', 'none', '0', '']:
                self.direction = 'neutral'
            elif direction_lower in ['complex', 'mixed']:
                self.direction = 'complex'
            else:
                # Keep original if not recognized
                pass
        
        # Default weight to 1.0 if not provided
        if self.weight is None:
            self.weight = 1.0


@dataclass
class Signature:
    """
    A complete lipid signature definition.
    
    Attributes:
        name: Signature name/identifier
        components: List of signature components
        description: Optional description
    """
    name: str
    components: List[SignatureComponent]
    description: Optional[str] = None


def load_signature_from_tsv(tsv_path: Path) -> Signature:
    """
    Load a signature from a TSV file.
    
    Expected TSV format:
        species	direction	weight
        Cer(d18:1/16:0)	↑	1.0
        Cer(d18:1/18:0)	↓	0.8
        SM(d18:1/16:0)	↑	1.2
    
    Columns:
        - species (required): Lipid species name
        - direction (optional): ↑, ↓, neutral, complex, or empty
        - weight (optional): Numeric weight (default: 1.0)
    
    Args:
        tsv_path: Path to TSV file
        
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
        
        species_col = None
        direction_col = None
        weight_col = None
        
        for key in ['species', 'metabolite', 'lipid', 'name']:
            if key in fieldnames_lower:
                species_col = fieldnames_lower[key]
                break
        
        for key in ['direction', 'dir', 'change']:
            if key in fieldnames_lower:
                direction_col = fieldnames_lower[key]
                break
        
        for key in ['weight', 'w', 'importance']:
            if key in fieldnames_lower:
                weight_col = fieldnames_lower[key]
                break
        
        if not species_col:
            raise ValueError(f"TSV must have a 'species' column: {tsv_path}")
        
        for row in reader:
            species = row.get(species_col, "").strip()
            if not species:
                continue  # Skip empty rows
            
            direction = row.get(direction_col, "").strip() if direction_col else None
            if not direction:
                direction = None
            
            weight_str = row.get(weight_col, "").strip() if weight_col else None
            weight = None
            if weight_str:
                try:
                    weight = float(weight_str)
                except ValueError:
                    weight = None  # Invalid weight, use default
            
            components.append(SignatureComponent(
                species=species,
                direction=direction,
                weight=weight,
            ))
    
    if not components:
        raise ValueError(f"No valid components found in signature: {tsv_path}")
    
    return Signature(
        name=signature_name,
        components=components,
    )


def load_signature_from_notion(signature_page_id: str) -> Signature:
    """
    Load a signature from a Notion page.
    
    This is a placeholder for future Notion integration.
    For now, signatures should be loaded from TSV files.
    
    Args:
        signature_page_id: Notion page ID for the signature
        
    Returns:
        Signature object
        
    Raises:
        NotImplementedError: Not yet implemented
    """
    raise NotImplementedError(
        "Notion signature loading not yet implemented. "
        "Please use TSV files for now."
    )

