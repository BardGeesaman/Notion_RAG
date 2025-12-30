"""
Flow cytometry ingestion and processing service.

This service orchestrates the complete flow cytometry workflow:
- FCS file ingestion and metadata extraction
- Background processing with event data transformation
- Gate application and population analysis
- Database persistence of results

Follows the pattern established by single_cell/ingest_service.py with threading-based
background processing for MVP deployment.
"""

from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional
from uuid import UUID

import numpy as np

from amprenta_rag.database.models import Dataset
from amprenta_rag.database.models_flow_cytometry import (
    FlowCytometryDataset,
    FlowCytometryGate,
    FlowCytometryParameter,
    FlowCytometryPopulation,
)
from amprenta_rag.database.session import db_session
from amprenta_rag.flow_cytometry.fcs_parser import (
    extract_metadata,
    load_fcs,
    load_events_parquet,
    save_events_parquet,
    validate_fcs,
)
from amprenta_rag.flow_cytometry.gating import (
    apply_polygon_gate,
    apply_quadrant_gate,
    apply_rectangle_gate,
    compute_population_stats,
)
from amprenta_rag.flow_cytometry.transforms import (
    apply_compensation,
    auto_logicle_params,
    logicle_transform,
    subsample_events,
)

logger = logging.getLogger(__name__)


class GateCreate:
    """Simplified gate creation schema for MVP."""
    
    def __init__(
        self,
        gate_name: str,
        gate_type: str,
        gate_definition: Dict,
        x_parameter_id: UUID,
        y_parameter_id: Optional[UUID] = None,
        parent_gate_id: Optional[UUID] = None,
        boolean_operator: Optional[str] = None,
        operand_gate_ids: Optional[List[UUID]] = None
    ):
        self.gate_name = gate_name
        self.gate_type = gate_type
        self.gate_definition = gate_definition
        self.x_parameter_id = x_parameter_id
        self.y_parameter_id = y_parameter_id
        self.parent_gate_id = parent_gate_id
        self.boolean_operator = boolean_operator
        self.operand_gate_ids = operand_gate_ids or []


def ingest_fcs(fcs_path: str, dataset_id: UUID | None = None) -> FlowCytometryDataset:
    """
    Ingest FCS file and create FlowCytometryDataset with background processing.
    
    Args:
        fcs_path: Path to FCS file
        dataset_id: Optional existing Dataset ID, creates new one if None
        
    Returns:
        FlowCytometryDataset instance
        
    Raises:
        FileNotFoundError: If FCS file doesn't exist
        ValueError: If FCS file is invalid
    """
    fcs_file = Path(fcs_path)
    if not fcs_file.exists():
        raise FileNotFoundError(f"FCS file not found: {fcs_path}")
    
    # Validate FCS file before processing
    issues = validate_fcs(fcs_path)
    if issues:
        raise ValueError(f"Invalid FCS file: {'; '.join(issues)}")
    
    with db_session() as db:
        ds_id = dataset_id
        if ds_id is None:
            # Create minimal Dataset row for FK relationship
            ds = Dataset(
                name=fcs_file.stem,
                omics_type="flow_cytometry",
                description=f"Imported from FCS file: {fcs_file.name}"
            )
            db.add(ds)
            db.commit()
            db.refresh(ds)
            ds_id = ds.id
        
        # Generate Parquet path for event storage
        parquet_path = f"flow_data/{ds_id}/{fcs_file.stem}_events.parquet"
        
        # Create FlowCytometryDataset record
        flow_dataset = FlowCytometryDataset(
            dataset_id=ds_id,
            events_parquet_path=parquet_path,
            file_size_bytes=fcs_file.stat().st_size,
            processing_status="pending",
            processing_log=None,
            ingested_at=datetime.now(timezone.utc)
        )
        db.add(flow_dataset)
        db.commit()
        db.refresh(flow_dataset)
        
        flow_dataset_id = flow_dataset.id
    
    # Start background processing
    thread = threading.Thread(
        target=_process_fcs_async,
        args=(flow_dataset_id, fcs_path),
        daemon=True
    )
    thread.start()
    
    # Return the created dataset
    with db_session() as db:
        result = db.query(FlowCytometryDataset).filter_by(id=flow_dataset_id).first()
        if not result:
            raise RuntimeError("Failed to load created FlowCytometryDataset")
        return result


def _process_fcs_async(flow_dataset_id: UUID, fcs_path: str) -> None:
    """
    Background processing of FCS file: parse, transform, save events.
    
    Args:
        flow_dataset_id: FlowCytometryDataset ID
        fcs_path: Path to original FCS file
    """
    try:
        with db_session() as db:
            flow_dataset = db.query(FlowCytometryDataset).filter_by(id=flow_dataset_id).first()
            if not flow_dataset:
                logger.error(f"FlowCytometryDataset {flow_dataset_id} not found")
                return
            
            flow_dataset.processing_status = "running"
            flow_dataset.processing_log = "Loading FCS file..."
            db.add(flow_dataset)
            db.commit()
        
        logger.info(f"Processing FCS file: {fcs_path}")
        
        # Load and parse FCS file
        meta, data = load_fcs(fcs_path)
        flow_metadata = extract_metadata(meta)
        
        logger.info(f"Loaded {data.shape[0]} events with {data.shape[1]} parameters")
        
        with db_session() as db:
            flow_dataset = db.query(FlowCytometryDataset).filter_by(id=flow_dataset_id).first()
            if not flow_dataset:
                return
            
            # Update dataset metadata
            flow_dataset.n_events = data.shape[0]
            flow_dataset.n_parameters = data.shape[1]
            flow_dataset.acquisition_date = flow_metadata.acquisition_date
            flow_dataset.cytometer_model = flow_metadata.cytometer_model
            flow_dataset.cytometer_serial = flow_metadata.cytometer_serial
            flow_dataset.acquisition_software = flow_metadata.acquisition_software
            flow_dataset.sample_id = flow_metadata.sample_id
            flow_dataset.sample_volume_ul = flow_metadata.sample_volume_ul
            flow_dataset.dilution_factor = flow_metadata.dilution_factor
            flow_dataset.staining_protocol = flow_metadata.staining_protocol
            flow_dataset.processing_log = "Creating parameter records..."
            db.add(flow_dataset)
            
            # Clear existing parameters (idempotent reprocessing)
            db.query(FlowCytometryParameter).filter_by(flow_dataset_id=flow_dataset_id).delete()
            
            # Create parameter records
            for i, param_name in enumerate(flow_metadata.parameter_names):
                param = FlowCytometryParameter(
                    flow_dataset_id=flow_dataset_id,
                    parameter_index=i,
                    parameter_name=param_name,
                    parameter_short_name=flow_metadata.parameter_short_names[i] if i < len(flow_metadata.parameter_short_names) else None,
                    fluorophore=flow_metadata.fluorophores[i] if i < len(flow_metadata.fluorophores) else None,
                    data_range=flow_metadata.data_ranges[i] if i < len(flow_metadata.data_ranges) else None,
                    min_value=float(np.min(data[:, i])),
                    max_value=float(np.max(data[:, i]))
                )
                db.add(param)
            
            db.commit()
        
        # Apply transformations if compensation matrix available
        transformed_data = data.copy()
        
        if flow_metadata.spillover_matrix is not None:
            logger.info("Applying compensation matrix")
            try:
                transformed_data = apply_compensation(transformed_data, flow_metadata.spillover_matrix)
            except Exception as e:
                logger.warning(f"Compensation failed: {e}")
        
        # Apply logicle transformation to fluorescence parameters
        # Assume parameters with names containing common fluorophores need transformation
        fluor_keywords = ['fitc', 'pe', 'apc', 'alexa', 'cy', 'pacific', 'brilliant']
        
        for i, param_name in enumerate(flow_metadata.parameter_names):
            param_lower = param_name.lower()
            if any(keyword in param_lower for keyword in fluor_keywords):
                logger.debug(f"Applying logicle transform to {param_name}")
                try:
                    # Auto-estimate logicle parameters
                    T, M, W, A = auto_logicle_params(transformed_data[:, i])
                    transformed_data[:, i] = logicle_transform(transformed_data[:, i], T, M, W, A)
                except Exception as e:
                    logger.warning(f"Logicle transform failed for {param_name}: {e}")
        
        # Subsample if dataset is very large (>100K events)
        if transformed_data.shape[0] > 100000:
            logger.info(f"Subsampling from {transformed_data.shape[0]} to 100000 events")
            transformed_data = subsample_events(transformed_data, max_events=100000)
        
        # Save events to Parquet
        parquet_path = Path("flow_data") / str(flow_dataset_id) / "events.parquet"
        parquet_path.parent.mkdir(parents=True, exist_ok=True)
        
        save_events_parquet(transformed_data, flow_metadata.parameter_names, str(parquet_path))
        
        # Update final status
        with db_session() as db:
            flow_dataset = db.query(FlowCytometryDataset).filter_by(id=flow_dataset_id).first()
            if flow_dataset:
                flow_dataset.events_parquet_path = str(parquet_path)
                flow_dataset.n_events = transformed_data.shape[0]  # Update if subsampled
                flow_dataset.processing_status = "completed"
                flow_dataset.processing_log = "Processing completed successfully"
                flow_dataset.processed_at = datetime.now(timezone.utc)
                db.add(flow_dataset)
                db.commit()
        
        logger.info(f"FCS processing completed for dataset {flow_dataset_id}")
        
    except Exception as e:
        logger.exception(f"FCS processing failed for dataset {flow_dataset_id}")
        
        with db_session() as db:
            flow_dataset = db.query(FlowCytometryDataset).filter_by(id=flow_dataset_id).first()
            if flow_dataset:
                flow_dataset.processing_status = "failed"
                flow_dataset.processing_log = str(e)[:4000]  # Truncate long error messages
                flow_dataset.processed_at = datetime.now(timezone.utc)
                db.add(flow_dataset)
                db.commit()


def apply_gate_to_dataset(flow_dataset_id: UUID, gate_definition: GateCreate) -> FlowCytometryGate:
    """
    Apply a gate to a flow cytometry dataset and compute population statistics.
    
    Args:
        flow_dataset_id: FlowCytometryDataset ID
        gate_definition: Gate configuration
        
    Returns:
        Created FlowCytometryGate instance
        
    Raises:
        ValueError: If dataset not found or gate parameters invalid
    """
    with db_session() as db:
        # Validate dataset exists and is processed
        flow_dataset = db.query(FlowCytometryDataset).filter_by(id=flow_dataset_id).first()
        if not flow_dataset:
            raise ValueError(f"FlowCytometryDataset {flow_dataset_id} not found")
        
        if flow_dataset.processing_status != "completed":
            raise ValueError(f"Dataset processing not completed (status: {flow_dataset.processing_status})")
        
        # Validate parameters exist
        x_param = db.query(FlowCytometryParameter).filter_by(id=gate_definition.x_parameter_id).first()
        if not x_param:
            raise ValueError(f"X parameter {gate_definition.x_parameter_id} not found")
        
        y_param = None
        if gate_definition.y_parameter_id:
            y_param = db.query(FlowCytometryParameter).filter_by(id=gate_definition.y_parameter_id).first()
            if not y_param:
                raise ValueError(f"Y parameter {gate_definition.y_parameter_id} not found")
        
        # Create gate record
        gate = FlowCytometryGate(
            flow_dataset_id=flow_dataset_id,
            gate_name=gate_definition.gate_name,
            gate_type=gate_definition.gate_type,
            gate_definition=gate_definition.gate_definition,
            x_parameter_id=gate_definition.x_parameter_id,
            y_parameter_id=gate_definition.y_parameter_id,
            parent_gate_id=gate_definition.parent_gate_id,
            boolean_operator=gate_definition.boolean_operator,
            operand_gate_ids=gate_definition.operand_gate_ids,
            is_active=True,
            created_at=datetime.now(timezone.utc)
        )
        db.add(gate)
        db.commit()
        db.refresh(gate)
        
        gate_id = gate.id
    
    # Apply gate and compute population in separate session
    _apply_gate_and_compute_population(gate_id)
    
    # Return the created gate
    with db_session() as db:
        result = db.query(FlowCytometryGate).filter_by(id=gate_id).first()
        if not result:
            raise RuntimeError("Failed to load created gate")
        return result


def _apply_gate_and_compute_population(gate_id: UUID) -> None:
    """
    Apply gate to events and compute population statistics.
    
    Args:
        gate_id: FlowCytometryGate ID
    """
    with db_session() as db:
        gate = db.query(FlowCytometryGate).filter_by(id=gate_id).first()
        if not gate:
            logger.error(f"Gate {gate_id} not found")
            return
        
        flow_dataset = gate.flow_dataset
        if not flow_dataset or not flow_dataset.events_parquet_path:
            logger.error(f"No event data found for dataset {gate.flow_dataset_id}")
            return
        
        # Load event data
        try:
            events, param_names = load_events_parquet(flow_dataset.events_parquet_path)
        except Exception as e:
            logger.error(f"Failed to load events from {flow_dataset.events_parquet_path}: {e}")
            return
        
        # Get parameter indices
        x_param = gate.x_parameter
        y_param = gate.y_parameter
        
        x_idx = x_param.parameter_index
        y_idx = y_param.parameter_index if y_param else None
        
        # Apply gate based on type
        try:
            if gate.gate_type == "polygon":
                vertices = gate.gate_definition["vertices"]
                mask = apply_polygon_gate(events, x_idx, y_idx, vertices)
                
            elif gate.gate_type == "rectangle":
                bounds = gate.gate_definition
                mask = apply_rectangle_gate(events, x_idx, y_idx, bounds)
                
            elif gate.gate_type == "quadrant":
                x_thresh = gate.gate_definition["x_threshold"]
                y_thresh = gate.gate_definition["y_threshold"]
                quadrant = gate.gate_definition.get("quadrant", "Q2")  # Default to Q2
                
                quadrants = apply_quadrant_gate(events, x_idx, y_idx, x_thresh, y_thresh)
                mask = quadrants[quadrant]
                
            else:
                logger.error(f"Unsupported gate type: {gate.gate_type}")
                return
                
        except Exception as e:
            logger.error(f"Gate application failed: {e}")
            return
        
        # Get parent mask if hierarchical gating
        parent_mask = None
        if gate.parent_gate_id:
            parent_population = db.query(FlowCytometryPopulation).filter_by(
                flow_dataset_id=gate.flow_dataset_id,
                gate_id=gate.parent_gate_id
            ).first()
            if parent_population:
                # For simplicity, assume parent mask covers all events for now
                # In full implementation, would need to reconstruct parent mask
                parent_mask = np.ones(len(events), dtype=bool)
        
        # Compute population statistics
        stats = compute_population_stats(
            events=events,
            mask=mask,
            param_names=param_names,
            parent_mask=parent_mask,
            total_events=len(events)
        )
        
        # Create or update population record
        existing_pop = db.query(FlowCytometryPopulation).filter_by(
            flow_dataset_id=gate.flow_dataset_id,
            gate_id=gate_id
        ).first()
        
        if existing_pop:
            # Update existing
            population = existing_pop
        else:
            # Create new
            population = FlowCytometryPopulation(
                flow_dataset_id=gate.flow_dataset_id,
                gate_id=gate_id
            )
            db.add(population)
        
        # Update population statistics
        population.event_count = stats.n_events
        population.parent_event_count = int(np.sum(parent_mask)) if parent_mask is not None else len(events)
        population.percentage_of_parent = stats.pct_of_parent
        population.percentage_of_total = stats.pct_of_total
        population.parameter_statistics = {
            "median": stats.median_values,
            "mean": stats.mean_values,
            "cv": stats.cv_values
        }
        population.analysis_software = "amprenta_rag"
        population.analysis_version = "1.0"
        population.analyzed_at = datetime.now(timezone.utc)
        
        db.commit()
        
        logger.info(f"Population computed for gate {gate_id}: {stats.n_events} events ({stats.pct_of_total:.1f}% of total)")


def recompute_populations(flow_dataset_id: UUID) -> None:
    """
    Recompute all populations for a flow cytometry dataset.
    
    Args:
        flow_dataset_id: FlowCytometryDataset ID
    """
    with db_session() as db:
        gates = db.query(FlowCytometryGate).filter_by(
            flow_dataset_id=flow_dataset_id,
            is_active=True
        ).all()
        
        gate_ids = [gate.id for gate in gates]
    
    # Recompute each gate in sequence
    for gate_id in gate_ids:
        try:
            _apply_gate_and_compute_population(gate_id)
        except Exception as e:
            logger.error(f"Failed to recompute population for gate {gate_id}: {e}")


__all__ = [
    "GateCreate",
    "ingest_fcs",
    "apply_gate_to_dataset", 
    "recompute_populations",
]
