# Compound-Target Network Guide

## Overview

The Compound-Target Network provides interactive visualization of drug discovery relationships using Cytoscape.js. It displays compounds and their biological targets as a graph, with edges weighted by bioactivity (IC50, EC50, Ki, Kd values). This enables visual exploration of polypharmacology, target selectivity, and structure-activity relationships across multiple targets.

**Use Cases:**
- Identify promiscuous compounds (many targets) vs. selective compounds (few targets)
- Discover target polypharmacology patterns for drug repurposing
- Visualize SAR series and their target coverage
- Find compounds with similar target profiles
- Prioritize selective hits for medicinal chemistry optimization

---

## Building Networks

### Initial Network Creation

1. **Start from Compound:**
   - Navigate to Compound Detail page
   - Click "View in Network" button
   - Network loads compound + connected targets

2. **Start from Target:**
   - Navigate to Target Detail page (if available)
   - Click "View in Network" button
   - Network loads target + active compounds

3. **Custom Query:**
   - Compound-Target Network page
   - Select compounds by IC50 range or structural similarity
   - Select targets by protein family or pathway
   - Click "Build Network"

### Expanding Networks

**Expand Compound Node:**
- Right-click compound → "Expand targets"
- Adds all targets with bioactivity for that compound

**Expand Target Node:**
- Right-click target → "Expand compounds"
- Adds all active compounds for that target

**Filters:**
- IC50 range slider (10 nM to 10 μM)
- Activity type checkboxes (IC50, EC50, Ki, Kd)
- Max nodes limit (default: 100, prevents layout explosion)

---

## Visualization

### Node Styling

**Compounds:**
- Shape: Ellipse
- Color: Light blue
- Size: Fixed (30px) or scaled by activity count
- Label: Compound name or SMILES

**Targets:**
- Shape: Diamond
- Color: Light green
- Size: Fixed (30px) or scaled by compound count
- Label: Gene symbol or protein name

### Edge Styling

**Width:** Scaled by pIC50 (thicker = more potent)
- pIC50 = -log10(IC50 in M)
- pIC50 > 7 (IC50 < 100 nM): thick edge
- pIC50 5-7: medium edge
- pIC50 < 5: thin edge

**Color:** By activity type
- IC50: Blue
- EC50: Green
- Ki: Orange
- Kd: Purple

**Hover:** Displays bioactivity value, units, and assay details

---

## Details Panel

Click any node or edge to view details in side panel:

**Compound Node:**
- SMILES structure (2D image)
- Molecular properties (MW, LogP, HBA/HBD)
- Target count and activity range
- "View Compound Details" button

**Target Node:**
- Gene symbol and protein name
- UniProt ID and pathway links
- Compound count and potency distribution
- "View Target Details" button (if available)

**Edge:**
- Activity value with units
- Assay type and conditions
- Publication reference (if available)
- ChEMBL assay ID link

---

## API Reference

### POST /api/network/compound-target

Build compound-target network from seed entities.

**Request:**
```json
{
  "compound_ids": ["uuid1", "uuid2"],
  "target_ids": ["uuid3"],
  "ic50_max": 1000,
  "activity_types": ["IC50", "Ki"],
  "max_nodes": 100
}
```

**Response:**
```json
{
  "nodes": [
    {"id": "uuid1", "type": "compound", "label": "Aspirin", "data": {...}},
    {"id": "uuid3", "type": "target", "label": "COX2", "data": {...}}
  ],
  "edges": [
    {"source": "uuid1", "target": "uuid3", "pic50": 7.2, "activity_type": "IC50"}
  ]
}
```

### GET /api/network/expand/compound/{compound_id}

Get all targets for a compound.

**Response:**
```json
{
  "targets": [{"id": "uuid", "name": "COX2", "pic50": 7.2}]
}
```

### GET /api/network/expand/target/{target_id}

Get all active compounds for a target.

**Response:**
```json
{
  "compounds": [{"id": "uuid", "smiles": "CC(=O)O...", "pic50": 7.2}]
}
```

---

## Integration

**Compound Detail Page:**
- "View in Network" button below structure image
- Opens Compound-Target Network with compound pre-loaded

**Virtual Screening Results:**
- "Network View" tab
- Visualize docked compounds and their predicted targets

**Signature Analysis:**
- "Target Network" for signature-associated targets
- Identify polypharmacology patterns in signature hits

---

## Related Documentation

- [Graph Explorer](USER_GUIDE.md) - General knowledge graph navigation
- [Virtual Screening Guide](VIRTUAL_SCREENING_GUIDE.md) - Docking workflow
- [Bioactivity Data](API_REFERENCE.md) - ChEMBL/PubChem integration

