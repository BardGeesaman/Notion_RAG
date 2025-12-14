#!/usr/bin/env python3
"""
Seed HTS test data for development/testing.

Creates:
- 1 HTSCampaign (TEST-HTS-001)
- 96 HTSResult records (A01-H12)
- Test compounds if they don't exist

Idempotent: skips if campaign already exists.
"""

import random
from datetime import datetime, timezone
from uuid import uuid4

from amprenta_rag.database.session import db_session
from amprenta_rag.models.chemistry import HTSCampaign, HTSResult, Compound
from amprenta_rag.chemistry.normalization import normalize_smiles, compute_molecular_descriptors


def generate_well_positions():
    """Generate all 96 well positions (A01-H12)."""
    wells = []
    for row in range(8):  # A-H
        row_letter = chr(ord('A') + row)
        for col in range(1, 13):  # 1-12
            wells.append(f"{row_letter}{col:02d}")
    return wells


def get_or_create_compound(db, compound_id: str, smiles: str) -> Compound:
    """Get existing compound or create new one."""
    # Check by compound_id first
    existing = db.query(Compound).filter(Compound.compound_id == compound_id).first()
    if existing:
        return existing
    
    # Normalize SMILES and compute descriptors
    canonical, inchi_key, formula = normalize_smiles(smiles)
    
    # Check by inchi_key (structure-based lookup)
    if inchi_key:
        existing_by_structure = db.query(Compound).filter(Compound.inchi_key == inchi_key).first()
        if existing_by_structure:
            return existing_by_structure
    
    # Only create if not found by either compound_id or inchi_key
    descriptors = compute_molecular_descriptors(canonical or smiles)
    
    compound = Compound(
        id=uuid4(),
        compound_id=compound_id,
        smiles=canonical or smiles,
        inchi_key=inchi_key,
        canonical_smiles=canonical,
        molecular_formula=formula,
        molecular_weight=descriptors.get("molecular_weight"),
        logp=descriptors.get("logp"),
        hbd_count=descriptors.get("hbd_count"),
        hba_count=descriptors.get("hba_count"),
        rotatable_bonds=descriptors.get("rotatable_bonds"),
        aromatic_rings=descriptors.get("aromatic_rings"),
    )
    db.add(compound)
    return compound


def seed_hts_data():
    """Seed HTS test data."""
    try:
        with db_session() as db:
            # Check if campaign already exists (idempotent)
            existing_campaign = db.query(HTSCampaign).filter(
                HTSCampaign.campaign_id == "TEST-HTS-001"
            ).first()
            
            if existing_campaign:
                print("Campaign TEST-HTS-001 already exists. Skipping seed.")
                return
            
            print("Creating HTS campaign TEST-HTS-001...")
            
            # Create campaign
            campaign = HTSCampaign(
                id=uuid4(),
                campaign_id="TEST-HTS-001",
                campaign_name="Test Kinase Inhibitor Screen",
                assay_type="Biochemical",
                target="CDK2",
                total_wells=96,
                hit_count=0,  # Will be updated after results
                run_date=datetime.now(timezone.utc),
            )
            db.add(campaign)
            db.flush()  # Get campaign.id
            
            # Generate well positions
            well_positions = generate_well_positions()
            
            # Create results with random values
            results = []
            normalized_values = []
            
            for i, well_pos in enumerate(well_positions):
                # Generate random normalized value (0-100, higher for potential hits)
                normalized_value = random.uniform(0, 100)
                normalized_values.append((i, normalized_value))
                
                # Correlated raw value (slightly different)
                raw_value = normalized_value * random.uniform(0.9, 1.1)
                
                # Random z-score (-2 to +3)
                z_score = random.uniform(-2, 3)
                
                # Create compound ID for this well
                compound_id = f"TEST-CMP-{i+1:03d}"
                smiles = "C" * (i % 10 + 4)  # Simple SMILES: CCCC, CCCCC, etc.
                
                # Get or create compound
                compound = get_or_create_compound(db, compound_id, smiles)
                db.flush()  # Ensure compound is persisted
                
                result = HTSResult(
                    id=uuid4(),
                    result_id=f"TEST-HTS-001-{well_pos}",
                    campaign_id=campaign.id,
                    compound_id=compound.id,
                    well_position=well_pos,
                    raw_value=round(raw_value, 2),
                    normalized_value=round(normalized_value, 2),
                    z_score=round(z_score, 2),
                    hit_flag=False,  # Will set top ~10 as hits
                    hit_category=None,
                )
                results.append(result)
                db.add(result)
            
            db.flush()
            
            # Mark top ~10 wells as hits (highest normalized_value)
            normalized_values.sort(key=lambda x: x[1], reverse=True)
            hit_count = min(10, len(normalized_values))
            
            hit_indices = {normalized_values[i][0] for i in range(hit_count)}
            
            for i, result in enumerate(results):
                if i in hit_indices:
                    result.hit_flag = True
                    result.hit_category = "Primary Hit"
            
            # Update campaign hit_count
            campaign.hit_count = hit_count
            
            db.commit()
            
            print(f"✓ Created campaign: {campaign.campaign_id}")
            print(f"✓ Created {len(results)} HTS results")
            print(f"✓ Marked {hit_count} wells as hits")
            print(f"✓ Created/verified {len(set(r.compound_id for r in results))} compounds")
    except Exception as e:
        print(f"Error seeding HTS data: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    seed_hts_data()

