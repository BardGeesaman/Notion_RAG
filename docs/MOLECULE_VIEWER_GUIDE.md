# 3D Molecule Viewer Guide

## Overview

The 3D Molecule Viewer provides interactive visualization of molecular structures and protein conformations using py3Dmol. It generates optimized 3D conformers from SMILES strings, overlays multiple molecules for comparison, and displays protein structures from PDB/AlphaFold databases.

**Use Cases:**
- Visualize compound structures before synthesis or purchase
- Compare binding poses from virtual screening results
- Overlay R-group variants in SAR analysis
- Inspect protein binding sites with ligands
- Generate publication-quality molecular images

---

## Conformer Generation

The module uses RDKit's **ETKDG (Experimental Torsion-angle Knowledge Distance Geometry)** method for 3D structure generation:

### Process

1. **Parse SMILES** - Convert 2D string to RDKit molecule object
2. **Add hydrogens** - Explicit hydrogens for accurate geometry
3. **Embed conformer** - ETKDG generates initial 3D coordinates
4. **Energy minimization** - MMFF94 force field optimization (200 iterations)
5. **Alignment** (optional) - Align multiple conformers to common scaffold

**Parameters:**
- Random seed for reproducibility
- Max iterations: 200 for MMFF
- Convergence threshold: 1e-4 kcal/mol

**Performance:** ~500ms per conformer on typical small molecules (<50 heavy atoms)

---

## Dashboard Pages

### Molecule Viewer Page

Standalone page for uploading and visualizing molecules:

1. **Input Methods:**
   - SMILES string input
   - Upload SDF/MOL file
   - Select from existing compounds

2. **Viewer Controls:**
   - Style selector: stick, sphere, line, cartoon
   - Color scheme: element colors, hydrophobicity, charge
   - Background: white, black, transparent
   - Spin toggle (auto-rotation)

3. **Actions:**
   - Download PDB file
   - Take screenshot (PNG)
   - Copy embed code for reports

### Compound 3D View (Expander)

Available on Compound Detail pages:

- Automatically loads structure from compound SMILES
- Embedded viewer with stick style default
- Click to open in full Molecule Viewer page

### Protein 3D View (Expander)

Available on Protein Structure pages:

- Displays PDB or AlphaFold structure
- Cartoon style for secondary structure
- Ligand overlay if co-crystal structure
- Binding site highlighting

---

## API Reference

### POST /api/viz3d/conformers

Generate 3D conformer from SMILES.

**Request:**
```json
{
  "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "seed": 42
}
```

**Response:**
```json
{
  "pdb": "ATOM    1  C1  UNL     1      ...",
  "success": true
}
```

### POST /api/viz3d/overlay

Generate and align multiple conformers.

**Request:**
```json
{
  "smiles_list": ["CCO", "CCC", "CCCC"],
  "align": true
}
```

**Response:**
```json
{
  "pdbs": ["ATOM 1 ...", "ATOM 1 ...", "ATOM 1 ..."],
  "success": true
}
```

### GET /api/viz3d/protein/{structure_id}

Retrieve protein structure in PDB format.

**Response:**
```json
{
  "pdb": "ATOM    1  N   MET A   1      ...",
  "source": "pdb" | "alphafold",
  "resolution": 2.1
}
```

---

## Styles

**stick** - Ball-and-stick representation (default for small molecules)

**sphere** - Space-filling CPK model (shows van der Waals surface)

**line** - Wire-frame (fast rendering for large molecules)

**cartoon** - Ribbon diagram (for proteins, shows secondary structure)

**Color schemes:**
- **element** - CPK colors (C=gray, O=red, N=blue)
- **hydrophobicity** - Gradient from hydrophilic (blue) to hydrophobic (orange)
- **charge** - Red (negative) to blue (positive)

---

## Troubleshooting

### py3Dmol Not Loading

**Symptoms:** Viewer shows blank or "No 3D view available"

**Causes:**
- JavaScript blocked by browser
- Streamlit running in restricted iframe
- py3Dmol CDN unavailable

**Solution:** Module automatically falls back to 2D structure image (RDKit SVG rendering)

### Invalid SMILES Handling

**Symptoms:** Error "Failed to generate conformer"

**Causes:**
- Malformed SMILES string
- RDKit cannot parse (e.g., organometallics, radicals)
- Kekulization failure

**Solution:** Try with `sanitize=False` or use 2D view only. Check SMILES validity in RDKit first.

### Slow Rendering (Large Molecules)

**Symptoms:** Viewer takes >5 seconds to load

**Causes:**
- Molecule has >100 heavy atoms
- Complex ring systems requiring many ETKDG iterations

**Solution:** Use "line" style instead of "stick/sphere", or skip energy minimization for faster (lower quality) conformers.

---

## Related Documentation

- [SAR Delta Explorer](VOILA_SAR_QUICKSTART.md) - R-group visualization
- [Virtual Screening Guide](VIRTUAL_SCREENING_GUIDE.md) - Docking pose visualization
- [Compound Management](USER_GUIDE.md) - Structure search and registration

