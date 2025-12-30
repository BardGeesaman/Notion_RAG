#!/usr/bin/env python3
"""
Seed biophysical assay demo data (idempotent).

Creates:
- SPR experiments with realistic kinetic parameters and sensorgrams
- MST experiments with dose-response curves and affinity data
- DSC experiments with thermograms and thermal stability parameters
- Demo compounds and targets for linking

Usage:
  python scripts/seed_biophysical_data.py
  python scripts/seed_biophysical_data.py --reset
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple
from uuid import uuid4, UUID

import numpy as np

# Ensure repo root is on sys.path when running as `python scripts/seed_biophysical_data.py`
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import User
from amprenta_rag.database.models_biophysical import (
    SPRExperiment,
    SPRSensorgram,
    MSTExperiment,
    MSTDoseResponse,
    DSCExperiment,
    DSCScan,
    SPRExperimentType,
)
from amprenta_rag.models.chemistry import Compound
from amprenta_rag.database.models import ProteinStructure


# Demo compounds for biophysical studies
DEMO_COMPOUNDS = [
    {"smiles": "CC(=O)Nc1ccc(O)cc1", "name": "Acetaminophen", "mw": 151.16},
    {"smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "name": "Ibuprofen", "mw": 206.28},
    {"smiles": "CC1=CC=C(C=C1)C(=O)C2=CC=CC=C2", "name": "Deoxybenzoin", "mw": 196.24},
    {"smiles": "CC(C)(C)NCC(c1ccc(O)c(CO)c1)O", "name": "Salbutamol", "mw": 239.31},
    {"smiles": "CN1CCN(CC1)c2ccc3nc4ccccc4nc3c2", "name": "Acridine Orange", "mw": 265.36},
]

# Demo targets for biophysical studies
DEMO_TARGETS = [
    {"name": "Human Serum Albumin", "sequence": "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL"},
    {"name": "Lysozyme", "sequence": "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"},
    {"name": "Carbonic Anhydrase II", "sequence": "MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNWRPAQPLKNRQIKASFK"},
    {"name": "Thrombin", "sequence": "TFGSGEADCGLRPLFEKKSLEDKTERELLESYIDGRIVEGSDAEIGMSPWQVMLFRKSPQELLCGASLISDRWVLTAAHCLLYPPWDKNFTENDLLVRIGKHSRTRYERNIEKISMLEKIYIHPRYNWRENLDRDIALMKLKKPVAFSDYIHPVCLPDRETAASLLQAGYKGRVTGWGNLKETWTANVGKGQPSVLQVVNLPIVERPVCKDSTRIRITDNMFCAGYKPDEGKRGDACEGDSGGPFVMKSPFNNRWYQMGIVSWGEGCDRDGKYGFYTHVFRLKKWIQKVIDQFGE"},
    {"name": "Trypsin", "sequence": "IVNGEEAVPGSWPWQVSLQDKTGFHFCGGSLINENWVVTAAHCGVTTSDVVVAGEFDQGSSSEKIQKLKIAKVFKNSKYNSLTINNDITLLKLSTAASFSQTVSAVCLPSASDDFAAGTTCVTTGWGLTRYTNANTPDRLQQASLPLLSNTNCKKYWGTKIKDAMICAGASGVSSCMGDSGGPLVCKKNGAWTLVGIVSWGSSTCSTSTPGVYARVTALVNWVQQTLAAN"},
]

# SPR experimental parameters
SPR_EXPERIMENTS = [
    {
        "name": "Compound-HSA Binding Study",
        "instrument": "Biacore T200",
        "chip_type": "CM5",
        "ka": 2.5e5,
        "kd_rate": 1.2e-3,
        "rmax": 150.0,
        "temperature": 25.0,
        "buffer": "HBS-EP+ (10 mM HEPES, 150 mM NaCl, 3 mM EDTA, 0.05% P20, pH 7.4)",
    },
    {
        "name": "Small Molecule-Lysozyme Interaction",
        "instrument": "Biacore 8K",
        "chip_type": "CM5",
        "ka": 1.8e4,
        "kd_rate": 3.5e-3,
        "rmax": 85.0,
        "temperature": 25.0,
        "buffer": "PBS-T (10 mM phosphate, 137 mM NaCl, 2.7 mM KCl, 0.05% Tween-20, pH 7.4)",
    },
    {
        "name": "Drug-Carbonic Anhydrase Binding",
        "instrument": "Biacore S200",
        "chip_type": "CM5",
        "ka": 5.2e5,
        "kd_rate": 8.7e-4,
        "rmax": 200.0,
        "temperature": 37.0,
        "buffer": "HBS-EP+ with 5% DMSO",
    },
]

# MST experimental parameters
MST_EXPERIMENTS = [
    {
        "name": "Thrombin Inhibitor Screening",
        "instrument": "Monolith NT.115",
        "kd_value": 2.5e-8,
        "hill_coefficient": 1.2,
        "temperature": 25.0,
        "buffer": "20 mM Tris-HCl, 150 mM NaCl, 0.05% Tween-20, pH 7.5",
    },
    {
        "name": "Trypsin Competitive Assay",
        "instrument": "Monolith NT.Automated",
        "kd_value": 1.8e-7,
        "hill_coefficient": 0.9,
        "temperature": 25.0,
        "buffer": "PBS with 0.1% BSA and 0.01% Triton X-100",
    },
    {
        "name": "HSA Drug Binding",
        "instrument": "Monolith NT.115",
        "kd_value": 4.2e-6,
        "hill_coefficient": 1.1,
        "temperature": 37.0,
        "buffer": "Physiological buffer (140 mM NaCl, 5 mM KCl, 1 mM MgCl2, pH 7.4)",
    },
]

# DSC experimental parameters
DSC_EXPERIMENTS = [
    {
        "name": "Lysozyme Thermal Stability",
        "instrument": "MicroCal VP-DSC",
        "tm_value": 72.5,
        "delta_h": -125.8,
        "cooperativity": 1.85,
        "buffer": "20 mM sodium phosphate, pH 7.0",
    },
    {
        "name": "HSA Denaturation Study",
        "instrument": "MicroCal Auto-DSC",
        "tm_value": 68.2,
        "delta_h": -98.3,
        "cooperativity": 1.45,
        "buffer": "PBS, pH 7.4",
    },
    {
        "name": "Carbonic Anhydrase Stability",
        "instrument": "TA Nano DSC",
        "tm_value": 64.8,
        "delta_h": -110.5,
        "cooperativity": 1.65,
        "buffer": "50 mM Tris-HCl, 100 mM NaCl, pH 8.0",
    },
]


def clear_biophysical_data(db) -> None:
    """Clear existing biophysical demo data."""
    print("🗑️  Clearing existing biophysical data...")
    
    # Delete in order to respect foreign key constraints
    db.query(SPRSensorgram).filter(SPRSensorgram.experiment_id.in_(
        db.query(SPRExperiment.id).filter(SPRExperiment.experiment_name.like("Demo %"))
    )).delete(synchronize_session=False)
    
    db.query(MSTDoseResponse).filter(MSTDoseResponse.experiment_id.in_(
        db.query(MSTExperiment.id).filter(MSTExperiment.experiment_name.like("Demo %"))
    )).delete(synchronize_session=False)
    
    db.query(DSCScan).filter(DSCScan.experiment_id.in_(
        db.query(DSCExperiment.id).filter(DSCExperiment.experiment_name.like("Demo %"))
    )).delete(synchronize_session=False)
    
    db.query(SPRExperiment).filter(SPRExperiment.experiment_name.like("Demo %")).delete(synchronize_session=False)
    db.query(MSTExperiment).filter(MSTExperiment.experiment_name.like("Demo %")).delete(synchronize_session=False)
    db.query(DSCExperiment).filter(DSCExperiment.experiment_name.like("Demo %")).delete(synchronize_session=False)
    
    # Note: Keeping existing compounds and targets to avoid schema field issues
    db.commit()
    print("✅ Cleared existing data")


def get_or_create_demo_user(db) -> UUID:
    """Get or create demo user for experiments."""
    user = db.query(User).filter(User.username == "demo_user").first()
    
    if not user:
        user = User(
            id=uuid4(),
            username="demo_user",
            email="demo@example.com",
            password_hash="fake_hash_for_demo"
        )
        db.add(user)
        db.commit()
        print("👤 Created demo user")
    
    user_id = user.id
    db.expunge(user)
    return user_id


def get_or_create_demo_compounds(db) -> List[UUID]:
    """Get or create demo compounds for biophysical studies."""
    compound_ids = []
    
    for comp_data in DEMO_COMPOUNDS:
        # Check if compound already exists by SMILES
        compound = db.query(Compound).filter(Compound.smiles == comp_data["smiles"]).first()
        
        if not compound:
            compound = Compound(
                id=uuid4(),
                compound_id=f"DEMO_{comp_data['name'].upper().replace(' ', '_')}_{uuid4().hex[:8]}",
                smiles=comp_data["smiles"],
                molecular_weight=comp_data["mw"],
                created_at=datetime.now(timezone.utc)
            )
            db.add(compound)
            db.commit()
            print(f"🧪 Created compound: {comp_data['name']}")
        
        compound_id = compound.id
        db.expunge(compound)
        compound_ids.append(compound_id)
    
    return compound_ids


def get_or_create_demo_targets(db) -> List[UUID]:
    """Get or create demo protein targets for biophysical studies."""
    target_ids = []
    
    # For simplicity, we'll just use any existing protein structures
    # In practice, these would be linked through proper protein/target models
    existing_targets = db.query(ProteinStructure).limit(len(DEMO_TARGETS)).all()
    
    if len(existing_targets) >= len(DEMO_TARGETS):
        for target in existing_targets[:len(DEMO_TARGETS)]:
            db.expunge(target)
            target_ids.append(target.id)
            print(f"🎯 Using existing target: {target.id}")
    else:
        # Create minimal demo targets if none exist
        for i, target_data in enumerate(DEMO_TARGETS[:3]):  # Limit to 3 for demo
            target = ProteinStructure(
                id=uuid4(),
                source="demo",
                created_at=datetime.now(timezone.utc)
            )
            db.add(target)
            db.commit()
            target_id = target.id
            db.expunge(target)
            target_ids.append(target_id)
            print(f"🎯 Created demo target: {target_id}")
    
    return target_ids


def seed_spr_experiments(db, compound_ids: List[UUID], target_ids: List[UUID], user_id: UUID) -> List[UUID]:
    """Create SPR experiments with sensorgrams."""
    spr_ids = []
    
    for i, spr_data in enumerate(SPR_EXPERIMENTS):
        compound_id = compound_ids[i % len(compound_ids)]
        target_id = target_ids[i % len(target_ids)]
        
        # Create SPR experiment
        experiment = SPRExperiment(
            id=uuid4(),
            experiment_name=f"Demo {spr_data['name']}",
            compound_id=compound_id,
            target_id=target_id,
            target_name=f"Demo Target {i + 1}",
            experiment_type=SPRExperimentType.KINETICS,
            instrument_model=spr_data["instrument"],
            chip_type=spr_data["chip_type"],
            temperature_celsius=spr_data["temperature"],
            buffer_composition=spr_data["buffer"],
            ka=spr_data["ka"],
            kd=spr_data["kd_rate"],
            kd_equilibrium=spr_data["kd_rate"] / spr_data["ka"],
            rmax=spr_data["rmax"],
            chi_squared=np.random.uniform(0.5, 2.5),
            processing_status="completed",
            response_units="RU",
            concentration_units="M",
            replicate_number=1,
            replicate_group_id=uuid4(),
            created_by_id=user_id,
            created_at=datetime.now(timezone.utc)
        )
        db.add(experiment)
        db.commit()
        
        experiment_id = experiment.id
        db.expunge(experiment)
        spr_ids.append(experiment_id)
        
        # Create sensorgrams for different concentrations
        concentrations = [0.1e-9, 1e-9, 10e-9, 100e-9, 1e-6]  # nM to M
        
        for cycle, conc in enumerate(concentrations, 1):
            # Generate realistic sensorgram data
            time_points = np.linspace(0, 600, 1000)  # 10 minutes
            response = simulate_spr_sensorgram(
                time_points, conc, spr_data["ka"], spr_data["kd_rate"], spr_data["rmax"]
            )
            
            sensorgram = SPRSensorgram(
                id=uuid4(),
                experiment_id=experiment_id,
                cycle_phase="association",
                time_seconds=[float(x) for x in time_points.tolist()],
                response_values=[float(x) for x in response.tolist()],
                created_at=datetime.now(timezone.utc)
            )
            db.add(sensorgram)
        
        db.commit()
        print(f"📈 Created SPR experiment: {spr_data['name']} with {len(concentrations)} sensorgrams")
    
    return spr_ids


def seed_mst_experiments(db, compound_ids: List[UUID], target_ids: List[UUID], user_id: UUID) -> List[UUID]:
    """Create MST experiments with dose-response curves."""
    mst_ids = []
    
    for i, mst_data in enumerate(MST_EXPERIMENTS):
        compound_id = compound_ids[i % len(compound_ids)]
        target_id = target_ids[i % len(target_ids)]
        
        # Create MST experiment
        experiment = MSTExperiment(
            id=uuid4(),
            experiment_name=f"Demo {mst_data['name']}",
            compound_id=compound_id,
            target_id=target_id,
            target_name=f"Demo Target {i + 1}",
            instrument_model=mst_data["instrument"],
            temperature_celsius=mst_data["temperature"],
            buffer_composition=mst_data["buffer"],
            kd_value=mst_data["kd_value"],
            kd_error=mst_data["kd_value"] * 0.15,  # 15% error
            hill_coefficient=mst_data["hill_coefficient"],
            binding_amplitude=np.random.uniform(15, 35),
            signal_to_noise=np.random.uniform(8, 20),
            aggregation_detected=False,
            processing_status="completed",
            response_units="‰",
            concentration_units="M",
            replicate_number=1,
            replicate_group_id=uuid4(),
            created_by_id=user_id,
            created_at=datetime.now(timezone.utc)
        )
        db.add(experiment)
        db.commit()
        
        experiment_id = experiment.id
        db.expunge(experiment)
        mst_ids.append(experiment_id)
        
        # Create dose-response points
        concentrations = np.logspace(-12, -6, 12)  # 1 pM to 1 μM
        
        for conc in concentrations:
            # Calculate Fnorm using Hill equation
            fnorm = calculate_hill_response(conc, mst_data["kd_value"], mst_data["hill_coefficient"])
            # Add realistic noise
            fnorm_noise = np.random.normal(0, 0.5)
            fnorm_error = np.random.uniform(0.3, 0.8)
            
            dose_point = MSTDoseResponse(
                id=uuid4(),
                experiment_id=experiment_id,
                ligand_concentration=float(conc),
                thermophoresis_response=float(fnorm + fnorm_noise),
                thermophoresis_error=float(fnorm_error),
                created_at=datetime.now(timezone.utc)
            )
            db.add(dose_point)
        
        db.commit()
        print(f"🌡️  Created MST experiment: {mst_data['name']} with {len(concentrations)} dose points")
    
    return mst_ids


def seed_dsc_experiments(db, compound_ids: List[UUID], target_ids: List[UUID], user_id: UUID) -> List[UUID]:
    """Create DSC experiments with thermograms."""
    dsc_ids = []
    
    for i, dsc_data in enumerate(DSC_EXPERIMENTS):
        compound_id = compound_ids[i % len(compound_ids)] if i < len(compound_ids) else None
        target_id = target_ids[i % len(target_ids)]
        
        # Create DSC experiment
        experiment = DSCExperiment(
            id=uuid4(),
            experiment_name=f"Demo {dsc_data['name']}",
            compound_id=compound_id,
            target_id=target_id,
            target_name=f"Demo Target {i + 1}",
            instrument_model=dsc_data["instrument"],
            buffer_composition=dsc_data["buffer"],
            tm_value=dsc_data["tm_value"],
            tm_error=np.random.uniform(0.2, 0.8),
            delta_h=dsc_data["delta_h"],
            delta_h_error=abs(dsc_data["delta_h"]) * 0.1,
            cooperativity=dsc_data["cooperativity"],
            processing_status="completed",
            response_units="kcal/mol/°C",
            concentration_units="mg/mL",
            replicate_number=1,
            replicate_group_id=uuid4(),
            created_by_id=user_id,
            created_at=datetime.now(timezone.utc)
        )
        db.add(experiment)
        db.commit()
        
        experiment_id = experiment.id
        db.expunge(experiment)
        dsc_ids.append(experiment_id)
        
        # Create thermogram scan
        temperatures = np.linspace(20, 90, 500)  # 20-90°C
        heat_capacity = simulate_dsc_thermogram(
            temperatures, dsc_data["tm_value"], dsc_data["delta_h"], dsc_data["cooperativity"]
        )
        
        scan = DSCScan(
            id=uuid4(),
            experiment_id=experiment_id,
            scan_type="sample",
            scan_number=1,
            temperature_celsius=[float(x) for x in temperatures.tolist()],
            heat_capacity=[float(x) for x in heat_capacity.tolist()],
            baseline_subtracted=True,
            created_at=datetime.now(timezone.utc)
        )
        db.add(scan)
        db.commit()
        
        print(f"🔥 Created DSC experiment: {dsc_data['name']} with thermogram")
    
    return dsc_ids


def simulate_spr_sensorgram(time: np.ndarray, conc: float, ka: float, kd: float, rmax: float) -> np.ndarray:
    """Simulate SPR sensorgram response."""
    response = np.zeros_like(time)
    association_end = 300  # seconds
    
    for i, t in enumerate(time):
        if t <= association_end:
            # Association phase
            kobs = conc * ka + kd
            response[i] = rmax * conc * ka / (conc * ka + kd) * (1 - np.exp(-kobs * t))
        else:
            # Dissociation phase
            req = rmax * conc * ka / (conc * ka + kd) * (1 - np.exp(-kobs * association_end))
            response[i] = req * np.exp(-kd * (t - association_end))
    
    # Add realistic noise
    noise = np.random.normal(0, 0.5, len(response))
    return response + noise


def calculate_hill_response(conc: float, kd: float, n: float) -> float:
    """Calculate Hill equation response."""
    return 1 / (1 + (kd / conc) ** n)


def simulate_dsc_thermogram(temperature: np.ndarray, tm: float, delta_h: float, cooperativity: float) -> np.ndarray:
    """Simulate DSC thermogram."""
    R = 1.987e-3  # Gas constant in kcal/mol/K
    T_kelvin = temperature + 273.15
    tm_kelvin = tm + 273.15
    
    # Van't Hoff equation with cooperativity
    K = np.exp(-delta_h * cooperativity / R * (1/T_kelvin - 1/tm_kelvin))
    
    # Fraction unfolded
    fu = K / (1 + K)
    
    # Heat capacity (derivative of enthalpy)
    cp_baseline = 1.5  # kcal/mol/°C
    cp_transition = (delta_h * cooperativity)**2 / (R * T_kelvin**2) * fu * (1 - fu)
    
    # Add baseline noise
    noise = np.random.normal(0, 0.05, len(temperature))
    return cp_baseline + cp_transition + noise


def main():
    """Main seeding function."""
    parser = argparse.ArgumentParser(description="Seed biophysical assay demo data")
    parser.add_argument("--reset", action="store_true", help="Clear existing data before seeding")
    args = parser.parse_args()
    
    print("🧬 Seeding biophysical assay demo data...")
    
    with db_session() as db:
        if args.reset:
            clear_biophysical_data(db)
        
        # Get or create demo entities
        user_id = get_or_create_demo_user(db)
        compound_ids = get_or_create_demo_compounds(db)
        target_ids = get_or_create_demo_targets(db)
        
        # Seed experiments
        spr_ids = seed_spr_experiments(db, compound_ids, target_ids, user_id)
        mst_ids = seed_mst_experiments(db, compound_ids, target_ids, user_id)
        dsc_ids = seed_dsc_experiments(db, compound_ids, target_ids, user_id)
        
        print(f"\n✅ Successfully created:")
        print(f"   📈 {len(spr_ids)} SPR experiments")
        print(f"   🌡️  {len(mst_ids)} MST experiments") 
        print(f"   🔥 {len(dsc_ids)} DSC experiments")
        print(f"   🧪 {len(compound_ids)} demo compounds")
        print(f"   🎯 {len(target_ids)} demo targets")
        print(f"\n🎉 Biophysical demo data seeding complete!")


if __name__ == "__main__":
    main()
