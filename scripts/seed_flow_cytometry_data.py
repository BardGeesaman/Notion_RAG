#!/usr/bin/env python3
"""
Seed flow cytometry demo data (idempotent).

Creates:
- A synthetic FCS dataset with 10K events and 12 parameters
- Sample gates (lymphocytes, CD4+, CD8+, etc.)
- Population statistics for each gate
- Realistic flow cytometry parameter ranges and distributions

Usage:
  python scripts/seed_flow_cytometry_data.py
  python scripts/seed_flow_cytometry_data.py --reset
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple
from uuid import uuid4

import numpy as np
import pandas as pd

# Ensure repo root is on sys.path when running as `python scripts/seed_flow_cytometry_data.py`
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, User
from amprenta_rag.database.models_flow_cytometry import (
    FlowCytometryDataset,
    FlowCytometryParameter,
    FlowCytometryGate,
    FlowCytometryPopulation,
    GateType,
)
from amprenta_rag.flow_cytometry.fcs_parser import save_events_parquet
from amprenta_rag.flow_cytometry.gating import (
    apply_polygon_gate,
    apply_rectangle_gate, 
    apply_quadrant_gate,
    compute_population_stats,
)

# Dataset constants
DATASET_NAME = "Demo Flow Cytometry Data"
DATASET_DESCRIPTION = "Synthetic flow cytometry data for demonstration and testing"
SAMPLE_ID = "DEMO_PBMC_001"
CYTOMETER_MODEL = "BD FACSCanto II"
N_EVENTS = 10000
N_PARAMETERS = 12

# Parameter definitions (typical flow cytometry panel)
PARAMETERS = [
    {"name": "FSC-A", "display": "Forward Scatter Area", "range_min": 0, "range_max": 262144},
    {"name": "FSC-H", "display": "Forward Scatter Height", "range_min": 0, "range_max": 262144},
    {"name": "SSC-A", "display": "Side Scatter Area", "range_min": 0, "range_max": 262144},
    {"name": "SSC-H", "display": "Side Scatter Height", "range_min": 0, "range_max": 262144},
    {"name": "FITC-A", "display": "CD3 FITC", "range_min": 0, "range_max": 262144},
    {"name": "PE-A", "display": "CD4 PE", "range_min": 0, "range_max": 262144},
    {"name": "PerCP-A", "display": "CD8 PerCP", "range_min": 0, "range_max": 262144},
    {"name": "APC-A", "display": "CD19 APC", "range_min": 0, "range_max": 262144},
    {"name": "PE-Cy7-A", "display": "CD14 PE-Cy7", "range_min": 0, "range_max": 262144},
    {"name": "FITC-H", "display": "CD3 FITC Height", "range_min": 0, "range_max": 262144},
    {"name": "PE-H", "display": "CD4 PE Height", "range_min": 0, "range_max": 262144},
    {"name": "Time", "display": "Time", "range_min": 0, "range_max": 1000},
]

# Gate definitions
SAMPLE_GATES = [
    {
        "name": "Lymphocytes",
        "type": GateType.POLYGON,
        "x_param": "FSC-A",
        "y_param": "SSC-A",
        "definition": {
            "vertices": [
                [20000, 15000], [80000, 15000], [120000, 60000], 
                [100000, 120000], [40000, 100000], [15000, 40000]
            ]
        }
    },
    {
        "name": "CD3+ T Cells",
        "type": GateType.RECTANGLE,
        "x_param": "FITC-A",
        "y_param": "SSC-A",
        "definition": {
            "x_min": 1000, "x_max": 100000,
            "y_min": 10000, "y_max": 150000
        }
    },
    {
        "name": "CD4+ T Helper",
        "type": GateType.QUADRANT,
        "x_param": "PE-A",
        "y_param": "PerCP-A",
        "definition": {
            "x_thresh": 1000,
            "y_thresh": 1000
        }
    },
    {
        "name": "B Cells",
        "type": GateType.RECTANGLE,
        "x_param": "APC-A",
        "y_param": "SSC-A",
        "definition": {
            "x_min": 2000, "x_max": 50000,
            "y_min": 15000, "y_max": 80000
        }
    }
]

# Stable prefixes for idempotent operation
DATASET_PREFIX = "FLOW-DEMO-"
GATE_PREFIX = "GATE-DEMO-"


def _now_utc() -> datetime:
    """Return current UTC datetime."""
    return datetime.now(timezone.utc)


def generate_synthetic_events(n_events: int) -> pd.DataFrame:
    """Generate synthetic flow cytometry events with realistic distributions."""
    print(f"Generating {n_events:,} synthetic events...")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    events = {}
    
    # Forward scatter (cell size) - log-normal distribution
    fsc_mean, fsc_std = 50000, 30000
    events["FSC-A"] = np.random.lognormal(np.log(fsc_mean), 0.8, n_events)
    events["FSC-H"] = events["FSC-A"] * np.random.uniform(0.8, 1.2, n_events)
    
    # Side scatter (cell granularity) - correlated with FSC but more variable
    ssc_base = events["FSC-A"] * np.random.uniform(0.3, 1.5, n_events)
    events["SSC-A"] = np.maximum(1000, ssc_base + np.random.normal(0, 15000, n_events))
    events["SSC-H"] = events["SSC-A"] * np.random.uniform(0.7, 1.3, n_events)
    
    # Fluorescence channels - mixture of negative and positive populations
    
    # CD3 (T cell marker) - ~70% of lymphocytes positive
    cd3_neg = np.random.lognormal(np.log(500), 0.5, int(n_events * 0.3))
    cd3_pos = np.random.lognormal(np.log(15000), 0.6, int(n_events * 0.7))
    events["FITC-A"] = np.concatenate([cd3_neg, cd3_pos])
    np.random.shuffle(events["FITC-A"])
    events["FITC-H"] = events["FITC-A"] * np.random.uniform(0.8, 1.2, n_events)
    
    # CD4 (T helper cells) - ~40% of T cells positive
    cd4_neg = np.random.lognormal(np.log(400), 0.4, int(n_events * 0.6))
    cd4_pos = np.random.lognormal(np.log(8000), 0.5, int(n_events * 0.4))
    events["PE-A"] = np.concatenate([cd4_neg, cd4_pos])
    np.random.shuffle(events["PE-A"])
    events["PE-H"] = events["PE-A"] * np.random.uniform(0.8, 1.2, n_events)
    
    # CD8 (Cytotoxic T cells) - ~30% of T cells positive
    cd8_neg = np.random.lognormal(np.log(600), 0.4, int(n_events * 0.7))
    cd8_pos = np.random.lognormal(np.log(12000), 0.5, int(n_events * 0.3))
    events["PerCP-A"] = np.concatenate([cd8_neg, cd8_pos])
    np.random.shuffle(events["PerCP-A"])
    
    # CD19 (B cell marker) - ~15% positive
    cd19_neg = np.random.lognormal(np.log(300), 0.3, int(n_events * 0.85))
    cd19_pos = np.random.lognormal(np.log(25000), 0.4, int(n_events * 0.15))
    events["APC-A"] = np.concatenate([cd19_neg, cd19_pos])
    np.random.shuffle(events["APC-A"])
    
    # CD14 (Monocyte marker) - ~10% positive
    cd14_neg = np.random.lognormal(np.log(400), 0.3, int(n_events * 0.9))
    cd14_pos = np.random.lognormal(np.log(30000), 0.4, int(n_events * 0.1))
    events["PE-Cy7-A"] = np.concatenate([cd14_neg, cd14_pos])
    np.random.shuffle(events["PE-Cy7-A"])
    
    # Time parameter (acquisition time)
    events["Time"] = np.linspace(0, 1000, n_events) + np.random.normal(0, 5, n_events)
    
    # Clip all values to parameter ranges
    for param in PARAMETERS:
        param_name = param["name"]
        if param_name in events:
            events[param_name] = np.clip(
                events[param_name], 
                param["range_min"], 
                param["range_max"]
            ).astype(np.float32)
    
    df = pd.DataFrame(events)
    print(f"✅ Generated events with shape: {df.shape}")
    return df


def create_flow_dataset(reset: bool = False) -> UUID:
    """Create or retrieve the demo flow cytometry dataset."""
    
    with db_session() as db:
        # Check if dataset already exists
        existing_dataset = db.query(Dataset).filter(
            Dataset.name.like(f"{DATASET_PREFIX}%")
        ).first()
        
        if existing_dataset and not reset:
            print(f"✅ Found existing dataset: {existing_dataset.name}")
            flow_dataset = db.query(FlowCytometryDataset).filter(
                FlowCytometryDataset.dataset_id == existing_dataset.id
            ).first()
            return flow_dataset.id if flow_dataset else None
        
        if existing_dataset and reset:
            print(f"🗑️  Removing existing dataset: {existing_dataset.name}")
            # Remove related flow cytometry data
            flow_dataset = db.query(FlowCytometryDataset).filter(
                FlowCytometryDataset.dataset_id == existing_dataset.id
            ).first()
            if flow_dataset:
                # Remove populations
                db.query(FlowCytometryPopulation).filter(
                    FlowCytometryPopulation.flow_dataset_id == flow_dataset.id
                ).delete()
                # Remove gates  
                db.query(FlowCytometryGate).filter(
                    FlowCytometryGate.flow_dataset_id == flow_dataset.id
                ).delete()
                # Remove parameters
                db.query(FlowCytometryParameter).filter(
                    FlowCytometryParameter.flow_dataset_id == flow_dataset.id
                ).delete()
                # Remove flow dataset
                db.delete(flow_dataset)
            
            # Remove main dataset
            db.delete(existing_dataset)
            db.commit()
        
        # Get or create a user for the dataset
        user = db.query(User).first()
        if not user:
            user = User(
                id=uuid4(),
                username="demo_user",
                email="demo@example.com",
                password_hash="demo_hash"
            )
            db.add(user)
            db.commit()
        
        # Create new dataset
        dataset_id = uuid4()
        dataset = Dataset(
            id=dataset_id,
            name=f"{DATASET_PREFIX}{datetime.now().strftime('%Y%m%d_%H%M%S')}",
            description=DATASET_DESCRIPTION,
            omics_type="flow_cytometry",
            created_by_id=user.id,  # Use FK field, not relationship field
            created_at=_now_utc()
        )
        db.add(dataset)
        db.commit()
        
        # Generate synthetic events
        events_df = generate_synthetic_events(N_EVENTS)
        
        # Save events to parquet
        parquet_dir = Path(tempfile.gettempdir()) / "flow_cytometry_demo"
        parquet_dir.mkdir(exist_ok=True)
        parquet_path = parquet_dir / f"events_{dataset_id}.parquet"
        
        save_events_parquet(events_df.values, list(events_df.columns), str(parquet_path))
        print(f"✅ Saved events to: {parquet_path}")
        
        # Create FlowCytometryDataset
        flow_dataset_id = uuid4()
        flow_dataset = FlowCytometryDataset(
            id=flow_dataset_id,
            dataset_id=dataset.id,
            events_parquet_path=str(parquet_path),
            file_size_bytes=parquet_path.stat().st_size,
            n_events=len(events_df),
            n_parameters=len(events_df.columns),
            cytometer_model=CYTOMETER_MODEL,
            sample_id=SAMPLE_ID,
            processing_status="completed",
            ingested_at=_now_utc(),
            processed_at=_now_utc()
        )
        db.add(flow_dataset)
        db.commit()
        
        # Create parameters
        for i, param in enumerate(PARAMETERS):
            parameter = FlowCytometryParameter(
                id=uuid4(),
                flow_dataset_id=flow_dataset.id,
                parameter_name=param["name"],
                parameter_index=i,
                min_value=float(param["range_min"]),
                max_value=float(param["range_max"])
            )
            db.add(parameter)
        
        db.commit()
        print(f"✅ Created dataset: {dataset.name}")
        print(f"✅ Created flow dataset with {N_EVENTS:,} events and {N_PARAMETERS} parameters")
        
        return flow_dataset.id


def create_sample_gates(flow_dataset_id: UUID, reset: bool = False) -> List[FlowCytometryGate]:
    """Create sample gates for the demo dataset."""
    
    with db_session() as db:
        # Get the flow dataset
        flow_dataset = db.query(FlowCytometryDataset).filter(
            FlowCytometryDataset.id == flow_dataset_id
        ).first()
        if not flow_dataset:
            raise ValueError(f"FlowCytometryDataset {flow_dataset_id} not found")
        
        # Check for existing gates
        existing_gates = db.query(FlowCytometryGate).filter(
            FlowCytometryGate.flow_dataset_id == flow_dataset.id,
            FlowCytometryGate.gate_name.like(f"{GATE_PREFIX}%")
        ).all()
        
        if existing_gates and not reset:
            print(f"✅ Found {len(existing_gates)} existing gates")
            return existing_gates
        
        if existing_gates and reset:
            print(f"🗑️  Removing {len(existing_gates)} existing gates")
            for gate in existing_gates:
                # Remove populations for this gate
                db.query(FlowCytometryPopulation).filter(
                    FlowCytometryPopulation.gate_id == gate.id
                ).delete()
                db.delete(gate)
            db.commit()
        
        # Load events data for population calculations
        events_df = pd.read_parquet(flow_dataset.events_parquet_path)
        total_events = len(events_df)
        
        # Get parameter lookup by name
        parameters = db.query(FlowCytometryParameter).filter(
            FlowCytometryParameter.flow_dataset_id == flow_dataset.id
        ).all()
        param_lookup = {param.parameter_name: param.id for param in parameters}
        
        created_gates = []
        
        for gate_def in SAMPLE_GATES:
            gate_id = uuid4()
            
            # Get parameter IDs
            x_param_id = param_lookup.get(gate_def["x_param"])
            y_param_id = param_lookup.get(gate_def["y_param"]) if gate_def.get("y_param") else None
            
            if not x_param_id:
                print(f"⚠️  Parameter {gate_def['x_param']} not found, skipping gate {gate_def['name']}")
                continue
            
            gate = FlowCytometryGate(
                id=gate_id,
                flow_dataset_id=flow_dataset.id,
                gate_name=f"{GATE_PREFIX}{gate_def['name']}",
                gate_type=gate_def["type"],
                x_parameter_id=x_param_id,
                y_parameter_id=y_param_id,
                gate_definition=gate_def["definition"],
                is_active=True,
                created_at=_now_utc()
            )
            db.add(gate)
            created_gates.append(gate)
        
        db.commit()
        print(f"✅ Created {len(created_gates)} sample gates")
        
        # Compute population statistics for each gate
        for gate in created_gates:
            compute_gate_population(gate, events_df, total_events)
        
        return created_gates


def compute_gate_population(gate: FlowCytometryGate, events_df: pd.DataFrame, total_events: int) -> None:
    """Compute population statistics for a gate."""
    
    with db_session() as db:
        # Get parameter names from the gate's parameter relationships
        x_param = db.query(FlowCytometryParameter).filter(
            FlowCytometryParameter.id == gate.x_parameter_id
        ).first()
        y_param = db.query(FlowCytometryParameter).filter(
            FlowCytometryParameter.id == gate.y_parameter_id
        ).first() if gate.y_parameter_id else None
        
        if not x_param:
            print(f"⚠️  X parameter not found for gate {gate.gate_name}")
            return
            
        # Apply gate to get boolean mask
        x_data = events_df[x_param.parameter_name].values
        y_data = events_df[y_param.parameter_name].values if y_param else None
        
        if gate.gate_type == GateType.POLYGON:
            if y_data is None:
                print(f"⚠️  Polygon gate {gate.gate_name} requires Y parameter")
                return
            vertices = gate.gate_definition["vertices"]
            mask = apply_polygon_gate(
                np.column_stack([x_data, y_data]),
                0, 1, vertices
            )
        elif gate.gate_type == GateType.RECTANGLE:
            if y_data is None:
                print(f"⚠️  Rectangle gate {gate.gate_name} requires Y parameter")
                return
            bounds = gate.gate_definition
            mask = apply_rectangle_gate(
                np.column_stack([x_data, y_data]),
                0, 1, bounds
            )
        elif gate.gate_type == GateType.QUADRANT:
            if y_data is None:
                print(f"⚠️  Quadrant gate {gate.gate_name} requires Y parameter")
                return
            # For quadrant gates, take Q1 (upper right)
            x_thresh = gate.gate_definition["x_thresh"]
            y_thresh = gate.gate_definition["y_thresh"]
            quadrants = apply_quadrant_gate(
                np.column_stack([x_data, y_data]),
                0, 1, x_thresh, y_thresh
            )
            mask = quadrants["Q1"]  # CD4+CD8- for typical flow panel
        else:
            print(f"⚠️  Unknown gate type: {gate.gate_type}")
            return
        
        # Compute population statistics
        param_names = list(events_df.columns)
        stats = compute_population_stats(
            events_df.values, mask, param_names, 
            parent_mask=None,  # No parent gate for demo gates
            total_events=total_events
        )
        
        # Create population record
        population = FlowCytometryPopulation(
            id=uuid4(),
            flow_dataset_id=gate.flow_dataset_id,
            gate_id=gate.id,
            population_name=f"{gate.gate_name} Population",
            event_count=stats.n_events,
            parent_event_count=total_events,
            percentage_of_parent=stats.pct_of_parent,
            percentage_of_total=stats.pct_of_total,
            parameter_statistics={
                param_name: {
                    "median": stats.median_values.get(param_name) if stats.median_values is not None else None,
                    "mean": stats.mean_values.get(param_name) if stats.mean_values is not None else None,
                    "cv": stats.cv_values.get(param_name) if stats.cv_values is not None else None
                }
                for param_name in param_names
            } if any([stats.median_values, stats.mean_values, stats.cv_values]) else None,
            analyzed_at=_now_utc()
        )
        db.add(population)
        db.commit()
        
        print(f"✅ Created population for {gate.gate_name}: {stats.n_events:,} events ({stats.pct_of_total:.1f}% of total)")


def main(reset: bool = False) -> None:
    """Main seeding function."""
    print("🧬 Seeding Flow Cytometry Demo Data")
    print("=" * 50)
    
    if reset:
        print("🔄 Reset mode: Will recreate all data")
    
    # Create dataset and events
    flow_dataset_id = create_flow_dataset(reset=reset)
    
    # Create sample gates and populations (pass ID to avoid DetachedInstanceError)
    gates = create_sample_gates(flow_dataset_id, reset=reset)
    
    # Get final dataset info for display
    with db_session() as db:
        flow_dataset = db.query(FlowCytometryDataset).filter(
            FlowCytometryDataset.id == flow_dataset_id
        ).first()
        
        dataset_id = flow_dataset.id
        n_events = flow_dataset.n_events
        n_parameters = flow_dataset.n_parameters
    
    print("\n" + "=" * 50)
    print("✅ Flow Cytometry Demo Data Seeding Complete!")
    print(f"📊 Dataset: {dataset_id}")
    print(f"📈 Events: {n_events:,}")
    print(f"🎯 Gates: {len(gates)}")
    print(f"🔬 Parameters: {n_parameters}")
    print("\nYou can now:")
    print("1. View data in the Flow Cytometry dashboard")
    print("2. Test API endpoints with the seeded data")
    print("3. Run integration tests")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Seed flow cytometry demo data")
    parser.add_argument(
        "--reset", 
        action="store_true", 
        help="Reset existing data (delete and recreate)"
    )
    args = parser.parse_args()
    
    try:
        main(reset=args.reset)
    except Exception as e:
        print(f"❌ Error seeding data: {e}")
        sys.exit(1)
