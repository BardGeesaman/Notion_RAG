"""
Biophysical assay ingestion and processing service.

This service orchestrates the complete biophysical workflow for SPR, MST, and DSC:
- File ingestion and metadata extraction
- Background processing with analysis pipelines
- Entity linking to compounds and targets
- Database persistence of results

Follows the pattern established by flow_cytometry/ingest_service.py with threading-based
background processing for MVP deployment.
"""

from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Union
from uuid import UUID, uuid4

import numpy as np

from amprenta_rag.database.models_biophysical import (
    DSCExperiment,
    DSCScan,
    MSTDoseResponse,
    MSTExperiment,
    SPRExperiment,
    SPRExperimentType,
    SPRSensorgram,
)
from amprenta_rag.database.session import db_session
from amprenta_rag.biophysical.spr_parser import parse_biacore_csv, parse_biacore_sensorgram_txt
from amprenta_rag.biophysical.mst_parser import parse_nanotemper_xlsx, parse_mst_csv
from amprenta_rag.biophysical.dsc_parser import parse_microcal_csv, parse_ta_instruments_txt
from amprenta_rag.biophysical.spr_analysis import fit_1_to_1_langmuir, global_fit
from amprenta_rag.biophysical.mst_analysis import fit_dose_response, quality_check
from amprenta_rag.biophysical.dsc_analysis import fit_two_state_unfolding, detect_peaks

logger = logging.getLogger(__name__)


def ingest_spr_file(
    file_path: str,
    experiment_id: Optional[UUID] = None,
    compound_id: Optional[UUID] = None,
    target_id: Optional[UUID] = None,
    target_name: Optional[str] = None,
    user_id: Optional[UUID] = None,
    db: Optional = None
) -> SPRExperiment:
    """
    Ingest SPR file: parse metadata, create experiment record, spawn background processing.
    
    Args:
        file_path: Path to SPR data file (CSV or TXT)
        experiment_id: Optional existing experiment ID
        compound_id: Optional compound to link to
        target_id: Optional target protein to link to
        target_name: Optional target name if target_id not provided
        user_id: User performing the ingestion
        db: Optional database session
        
    Returns:
        SPRExperiment object with processing_status='pending'
        
    Raises:
        ValueError: If file is invalid or required parameters missing
        FileNotFoundError: If file doesn't exist
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"SPR file not found: {file_path}")
    
    logger.info(f"Starting SPR ingestion: {file_path}")
    
    # Determine file format and parse metadata
    try:
        if file_path.suffix.lower() == '.csv':
            spr_data = parse_biacore_csv(str(file_path))
        elif file_path.suffix.lower() == '.txt':
            spr_data = parse_biacore_sensorgram_txt(str(file_path))
        else:
            raise ValueError(f"Unsupported SPR file format: {file_path.suffix}")
    except Exception as e:
        logger.error(f"Failed to parse SPR file {file_path}: {e}")
        raise ValueError(f"Invalid SPR file format: {e}") from e
    
    # Determine experiment type based on data
    if len(spr_data.sensorgrams) > 1:
        experiment_type = SPRExperimentType.KINETICS
    else:
        experiment_type = SPRExperimentType.AFFINITY
    
    # Create database session if not provided
    close_db = db is None
    if db is None:
        db = db_session()
    
    try:
        # Create SPR experiment record
        spr_experiment = SPRExperiment(
            id=experiment_id or uuid4(),
            experiment_name=f"SPR_{file_path.stem}",
            compound_id=compound_id,
            target_id=target_id,
            target_name=target_name or spr_data.ligand_name,
            experiment_type=experiment_type,
            ligand_name=spr_data.ligand_name,
            temperature_celsius=spr_data.temperature,
            buffer_composition=spr_data.buffer,
            chip_type=spr_data.chip_type,
            flow_rate=spr_data.flow_rate,
            response_units="RU",
            concentration_units="M",
            replicate_number=1,
            replicate_group_id=experiment_id,
            processing_status="pending",
            created_by_id=user_id
        )
        
        db.add(spr_experiment)
        db.commit()
        
        # Get the experiment ID for background processing
        experiment_uuid = spr_experiment.id
        
        # Detach from session before returning
        db.expunge(spr_experiment)
        
        # Start background processing
        processing_thread = threading.Thread(
            target=_process_spr_async,
            args=(experiment_uuid,),
            daemon=True
        )
        processing_thread.start()
        
        logger.info(f"SPR experiment created: {experiment_uuid}, processing started")
        return spr_experiment
        
    except Exception as e:
        logger.error(f"Failed to create SPR experiment: {e}")
        db.rollback()
        raise
    finally:
        if close_db:
            db.close()


def ingest_mst_file(
    file_path: str,
    experiment_id: Optional[UUID] = None,
    compound_id: Optional[UUID] = None,
    target_id: Optional[UUID] = None,
    target_name: Optional[str] = None,
    user_id: Optional[UUID] = None,
    db: Optional = None
) -> MSTExperiment:
    """
    Ingest MST file: parse metadata, create experiment record, spawn background processing.
    
    Args:
        file_path: Path to MST data file (XLSX or CSV)
        experiment_id: Optional existing experiment ID
        compound_id: Optional compound to link to
        target_id: Optional target protein to link to
        target_name: Optional target name if target_id not provided
        user_id: User performing the ingestion
        db: Optional database session
        
    Returns:
        MSTExperiment object with processing_status='pending'
        
    Raises:
        ValueError: If file is invalid or required parameters missing
        FileNotFoundError: If file doesn't exist
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"MST file not found: {file_path}")
    
    logger.info(f"Starting MST ingestion: {file_path}")
    
    # Determine file format and parse metadata
    try:
        if file_path.suffix.lower() in ['.xlsx', '.xls']:
            mst_data = parse_nanotemper_xlsx(str(file_path))
        elif file_path.suffix.lower() == '.csv':
            mst_data = parse_mst_csv(str(file_path))
        else:
            raise ValueError(f"Unsupported MST file format: {file_path.suffix}")
    except Exception as e:
        logger.error(f"Failed to parse MST file {file_path}: {e}")
        raise ValueError(f"Invalid MST file format: {e}") from e
    
    # Create database session if not provided
    close_db = db is None
    if db is None:
        db = db_session()
    
    try:
        # Create MST experiment record
        mst_experiment = MSTExperiment(
            id=experiment_id or uuid4(),
            experiment_name=f"MST_{file_path.stem}",
            compound_id=compound_id,
            target_id=target_id,
            target_name=target_name or "Unknown Target",
            temperature_celsius=mst_data.temperature,
            buffer_composition=mst_data.buffer,
            capillary_type=mst_data.capillary_type,
            mst_power=mst_data.mst_power,
            excitation_power=mst_data.led_power,
            response_units="‰",
            concentration_units="M",
            replicate_number=1,
            replicate_group_id=experiment_id,
            processing_status="pending",
            created_by_id=user_id
        )
        
        db.add(mst_experiment)
        db.commit()
        
        # Get the experiment ID for background processing
        experiment_uuid = mst_experiment.id
        
        # Detach from session before returning
        db.expunge(mst_experiment)
        
        # Start background processing
        processing_thread = threading.Thread(
            target=_process_mst_async,
            args=(experiment_uuid,),
            daemon=True
        )
        processing_thread.start()
        
        logger.info(f"MST experiment created: {experiment_uuid}, processing started")
        return mst_experiment
        
    except Exception as e:
        logger.error(f"Failed to create MST experiment: {e}")
        db.rollback()
        raise
    finally:
        if close_db:
            db.close()


def ingest_dsc_file(
    file_path: str,
    experiment_id: Optional[UUID] = None,
    compound_id: Optional[UUID] = None,
    protein_id: Optional[UUID] = None,
    protein_name: Optional[str] = None,
    user_id: Optional[UUID] = None,
    db: Optional = None
) -> DSCExperiment:
    """
    Ingest DSC file: parse metadata, create experiment record, spawn background processing.
    
    Args:
        file_path: Path to DSC data file (CSV or TXT)
        experiment_id: Optional existing experiment ID
        compound_id: Optional compound to link to (for ligand binding studies)
        protein_id: Optional protein to link to (same as target_id)
        protein_name: Optional protein name if protein_id not provided
        user_id: User performing the ingestion
        db: Optional database session
        
    Returns:
        DSCExperiment object with processing_status='pending'
        
    Raises:
        ValueError: If file is invalid or required parameters missing
        FileNotFoundError: If file doesn't exist
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"DSC file not found: {file_path}")
    
    logger.info(f"Starting DSC ingestion: {file_path}")
    
    # Determine file format and parse metadata
    try:
        if file_path.suffix.lower() == '.csv':
            dsc_data = parse_microcal_csv(str(file_path))
        elif file_path.suffix.lower() == '.txt':
            dsc_data = parse_ta_instruments_txt(str(file_path))
        else:
            raise ValueError(f"Unsupported DSC file format: {file_path.suffix}")
    except Exception as e:
        logger.error(f"Failed to parse DSC file {file_path}: {e}")
        raise ValueError(f"Invalid DSC file format: {e}") from e
    
    # Create database session if not provided
    close_db = db is None
    if db is None:
        db = db_session()
    
    try:
        # Create DSC experiment record
        dsc_experiment = DSCExperiment(
            id=experiment_id or uuid4(),
            experiment_name=f"DSC_{file_path.stem}",
            compound_id=compound_id,
            target_id=protein_id,  # DSC uses target_id for protein
            target_name=protein_name or dsc_data.protein_name,
            scan_rate=dsc_data.scan_rate,
            protein_concentration=dsc_data.protein_concentration,
            buffer_composition=dsc_data.buffer,
            response_units="kcal/mol/°C",
            concentration_units="mg/mL",
            replicate_number=1,
            replicate_group_id=experiment_id,
            processing_status="pending",
            created_by_id=user_id
        )
        
        db.add(dsc_experiment)
        db.commit()
        
        # Get the experiment ID for background processing
        experiment_uuid = dsc_experiment.id
        
        # Detach from session before returning
        db.expunge(dsc_experiment)
        
        # Start background processing
        processing_thread = threading.Thread(
            target=_process_dsc_async,
            args=(experiment_uuid,),
            daemon=True
        )
        processing_thread.start()
        
        logger.info(f"DSC experiment created: {experiment_uuid}, processing started")
        return dsc_experiment
        
    except Exception as e:
        logger.error(f"Failed to create DSC experiment: {e}")
        db.rollback()
        raise
    finally:
        if close_db:
            db.close()


def _process_spr_async(spr_id: UUID) -> None:
    """
    Background task for SPR processing.
    
    1. Parse file (spr_parser)
    2. Run kinetic analysis (spr_analysis)
    3. Store sensorgrams in DB
    4. Update processing_status='completed' or 'failed'
    
    Args:
        spr_id: UUID of the SPR experiment to process
    """
    logger.info(f"Starting background SPR processing: {spr_id}")
    
    with db_session() as db:
        try:
            # Re-query experiment
            spr_experiment = db.query(SPRExperiment).filter(SPRExperiment.id == spr_id).first()
            if not spr_experiment:
                logger.error(f"SPR experiment not found: {spr_id}")
                return
            
            # Update status to running
            spr_experiment.processing_status = "running"
            db.commit()
            
            # Parse file again for full data
            file_path = spr_experiment.raw_data_path
            if file_path.endswith('.csv'):
                spr_data = parse_biacore_csv(file_path)
            else:
                spr_data = parse_biacore_sensorgram_txt(file_path)
            
            # Process each sensorgram
            kinetic_fits = []
            for i, sensorgram in enumerate(spr_data.sensorgrams):
                try:
                    # Run kinetic analysis
                    fit = fit_1_to_1_langmuir(sensorgram)
                    kinetic_fits.append(fit)
                    
                    # Store sensorgram in database
                    spr_sensorgram = SPRSensorgram(
                        experiment_id=spr_id,
                        cycle_phase="association",  # Simplified for now
                        time_seconds=sensorgram.time.tolist(),
                        response_values=sensorgram.response.tolist(),
                        flow_cell=sensorgram.flow_cell,
                        reference_subtracted=sensorgram.reference_subtracted
                    )
                    db.add(spr_sensorgram)
                    
                except Exception as e:
                    logger.warning(f"Failed to process sensorgram {i}: {e}")
                    continue
            
            # Update experiment with analysis results
            if kinetic_fits:
                # Use first fit or global fit if multiple concentrations
                if len(kinetic_fits) > 1:
                    # Perform global fitting
                    try:
                        global_fit_result = global_fit(spr_data.sensorgrams)
                        spr_experiment.ka = global_fit_result.ka
                        spr_experiment.kd = global_fit_result.kd
                        spr_experiment.kd_equilibrium = global_fit_result.kd_affinity
                    except Exception as e:
                        logger.warning(f"Global fitting failed, using first fit: {e}")
                        best_fit = kinetic_fits[0]
                        spr_experiment.ka = best_fit.ka
                        spr_experiment.kd = best_fit.kd
                        spr_experiment.kd_equilibrium = best_fit.kd_affinity
                else:
                    best_fit = kinetic_fits[0]
                    spr_experiment.ka = best_fit.ka
                    spr_experiment.kd = best_fit.kd
                    spr_experiment.kd_equilibrium = best_fit.kd_affinity
                
                # Update other parameters
                spr_experiment.rmax = np.mean([fit.rmax for fit in kinetic_fits])
                spr_experiment.chi_squared = np.mean([fit.chi_squared for fit in kinetic_fits])
            
            # Mark as completed
            spr_experiment.processing_status = "completed"
            spr_experiment.processed_at = datetime.now(timezone.utc)
            
            db.commit()
            logger.info(f"SPR processing completed: {spr_id}")
            
        except Exception as e:
            logger.error(f"SPR processing failed for {spr_id}: {e}")
            try:
                # Mark as failed
                spr_experiment = db.query(SPRExperiment).filter(SPRExperiment.id == spr_id).first()
                if spr_experiment:
                    spr_experiment.processing_status = "failed"
                    spr_experiment.error_message = str(e)
                    db.commit()
            except Exception as commit_error:
                logger.error(f"Failed to update error status: {commit_error}")


def _process_mst_async(mst_id: UUID) -> None:
    """
    Background task for MST processing.
    
    1. Parse file (mst_parser)
    2. Run affinity analysis (mst_analysis)
    3. Store dose-response data in DB
    4. Update processing_status='completed' or 'failed'
    
    Args:
        mst_id: UUID of the MST experiment to process
    """
    logger.info(f"Starting background MST processing: {mst_id}")
    
    with db_session() as db:
        try:
            # Re-query experiment
            mst_experiment = db.query(MSTExperiment).filter(MSTExperiment.id == mst_id).first()
            if not mst_experiment:
                logger.error(f"MST experiment not found: {mst_id}")
                return
            
            # Update status to running
            mst_experiment.processing_status = "running"
            db.commit()
            
            # Parse file again for full data
            file_path = mst_experiment.raw_data_path
            if file_path.endswith('.xlsx') or file_path.endswith('.xls'):
                mst_data = parse_nanotemper_xlsx(file_path)
            else:
                mst_data = parse_mst_csv(file_path)
            
            # Process dose-response data
            concentrations = np.array([point.concentration for point in mst_data.dose_points])
            fnorm_values = np.array([point.fnorm for point in mst_data.dose_points])
            errors = np.array([point.fnorm_error for point in mst_data.dose_points if point.fnorm_error])
            
            # Run affinity analysis
            try:
                if len(errors) == len(concentrations):
                    affinity_fit = fit_dose_response(concentrations, fnorm_values, errors)
                else:
                    affinity_fit = fit_dose_response(concentrations, fnorm_values)
                
                # Quality check
                quality_metrics = quality_check(affinity_fit, concentrations, fnorm_values)
                
                # Store dose-response points in database
                for point in mst_data.dose_points:
                    mst_dose_point = MSTDoseResponse(
                        experiment_id=mst_id,
                        ligand_concentration=point.concentration,
                        thermophoresis_response=point.fnorm,
                        thermophoresis_error=point.fnorm_error,
                        initial_fluorescence=point.cold_fluorescence,
                        thermophoresis_signal=point.hot_fluorescence
                    )
                    db.add(mst_dose_point)
                
                # Update experiment with analysis results
                mst_experiment.kd_value = affinity_fit.kd
                mst_experiment.kd_error = affinity_fit.kd_error
                mst_experiment.hill_coefficient = affinity_fit.hill_coefficient
                mst_experiment.binding_amplitude = affinity_fit.amplitude
                mst_experiment.r_squared = affinity_fit.r_squared
                mst_experiment.signal_to_noise = quality_metrics.signal_to_noise
                
            except Exception as e:
                logger.warning(f"MST analysis failed: {e}")
                # Still store raw data points
                for point in mst_data.dose_points:
                    mst_dose_point = MSTDoseResponse(
                        experiment_id=mst_id,
                        ligand_concentration=point.concentration,
                        thermophoresis_response=point.fnorm,
                        thermophoresis_error=point.fnorm_error,
                        initial_fluorescence=point.cold_fluorescence,
                        thermophoresis_signal=point.hot_fluorescence
                    )
                    db.add(mst_dose_point)
            
            # Mark as completed
            mst_experiment.processing_status = "completed"
            mst_experiment.processed_at = datetime.now(timezone.utc)
            
            db.commit()
            logger.info(f"MST processing completed: {mst_id}")
            
        except Exception as e:
            logger.error(f"MST processing failed for {mst_id}: {e}")
            try:
                # Mark as failed
                mst_experiment = db.query(MSTExperiment).filter(MSTExperiment.id == mst_id).first()
                if mst_experiment:
                    mst_experiment.processing_status = "failed"
                    mst_experiment.error_message = str(e)
                    db.commit()
            except Exception as commit_error:
                logger.error(f"Failed to update error status: {commit_error}")


def _process_dsc_async(dsc_id: UUID) -> None:
    """
    Background task for DSC processing.
    
    1. Parse file (dsc_parser)
    2. Run thermal analysis (dsc_analysis)
    3. Store thermogram data in DB
    4. Update processing_status='completed' or 'failed'
    
    Args:
        dsc_id: UUID of the DSC experiment to process
    """
    logger.info(f"Starting background DSC processing: {dsc_id}")
    
    with db_session() as db:
        try:
            # Re-query experiment
            dsc_experiment = db.query(DSCExperiment).filter(DSCExperiment.id == dsc_id).first()
            if not dsc_experiment:
                logger.error(f"DSC experiment not found: {dsc_id}")
                return
            
            # Update status to running
            dsc_experiment.processing_status = "running"
            db.commit()
            
            # Parse file again for full data
            file_path = dsc_experiment.raw_data_path
            if file_path.endswith('.csv'):
                dsc_data = parse_microcal_csv(file_path)
            else:
                dsc_data = parse_ta_instruments_txt(file_path)
            
            # Process each scan
            thermal_fits = []
            for i, scan in enumerate(dsc_data.scans):
                try:
                    # Run thermal analysis
                    thermal_fit = fit_two_state_unfolding(scan.temperature, scan.heat_capacity)
                    thermal_fits.append(thermal_fit)
                    
                    # Detect peaks
                    peaks = detect_peaks(scan.temperature, scan.heat_capacity)
                    
                    # Store scan in database
                    dsc_scan = DSCScan(
                        experiment_id=dsc_id,
                        scan_number=scan.scan_number,
                        temperature_celsius=scan.temperature.tolist(),
                        heat_capacity=scan.heat_capacity.tolist(),
                        scan_type=scan.scan_type or "sample"
                    )
                    db.add(dsc_scan)
                    
                except Exception as e:
                    logger.warning(f"Failed to process DSC scan {i}: {e}")
                    # Still store raw scan data
                    dsc_scan = DSCScan(
                        experiment_id=dsc_id,
                        scan_number=scan.scan_number,
                        temperature_celsius=scan.temperature.tolist(),
                        heat_capacity=scan.heat_capacity.tolist(),
                        scan_type=scan.scan_type or "sample"
                    )
                    db.add(dsc_scan)
                    continue
            
            # Update experiment with analysis results
            if thermal_fits:
                # Use best fit (highest cooperativity or lowest error)
                best_fit = max(thermal_fits, key=lambda f: f.cooperativity)
                
                dsc_experiment.tm_value = best_fit.tm
                dsc_experiment.tm_error = best_fit.tm_error
                dsc_experiment.delta_h = best_fit.delta_h
                dsc_experiment.delta_cp = best_fit.delta_cp
                dsc_experiment.cooperativity = best_fit.cooperativity
            
            # Mark as completed
            dsc_experiment.processing_status = "completed"
            dsc_experiment.processed_at = datetime.now(timezone.utc)
            
            db.commit()
            logger.info(f"DSC processing completed: {dsc_id}")
            
        except Exception as e:
            logger.error(f"DSC processing failed for {dsc_id}: {e}")
            try:
                # Mark as failed
                dsc_experiment = db.query(DSCExperiment).filter(DSCExperiment.id == dsc_id).first()
                if dsc_experiment:
                    dsc_experiment.processing_status = "failed"
                    dsc_experiment.error_message = str(e)
                    db.commit()
            except Exception as commit_error:
                logger.error(f"Failed to update error status: {commit_error}")


def link_to_compound(
    experiment_id: UUID,
    compound_id: UUID,
    assay_type: str
) -> None:
    """
    Link biophysical experiment to a compound.
    
    Args:
        experiment_id: UUID of the experiment
        compound_id: UUID of the compound to link
        assay_type: Type of assay ("spr", "mst", "dsc")
        
    Raises:
        ValueError: If assay_type is invalid or experiment not found
    """
    if assay_type not in ["spr", "mst", "dsc"]:
        raise ValueError(f"Invalid assay type: {assay_type}")
    
    logger.info(f"Linking {assay_type.upper()} experiment {experiment_id} to compound {compound_id}")
    
    with db_session() as db:
        try:
            if assay_type == "spr":
                experiment = db.query(SPRExperiment).filter(SPRExperiment.id == experiment_id).first()
            elif assay_type == "mst":
                experiment = db.query(MSTExperiment).filter(MSTExperiment.id == experiment_id).first()
            else:  # dsc
                experiment = db.query(DSCExperiment).filter(DSCExperiment.id == experiment_id).first()
            
            if not experiment:
                raise ValueError(f"{assay_type.upper()} experiment not found: {experiment_id}")
            
            experiment.compound_id = compound_id
            db.commit()
            
            logger.info(f"Successfully linked {assay_type.upper()} experiment to compound")
            
        except Exception as e:
            logger.error(f"Failed to link experiment to compound: {e}")
            db.rollback()
            raise


def link_to_target(
    experiment_id: UUID,
    target_id: UUID,
    assay_type: str
) -> None:
    """
    Link biophysical experiment to a target protein.
    
    Args:
        experiment_id: UUID of the experiment
        target_id: UUID of the target protein to link
        assay_type: Type of assay ("spr", "mst", "dsc")
        
    Raises:
        ValueError: If assay_type is invalid or experiment not found
    """
    if assay_type not in ["spr", "mst", "dsc"]:
        raise ValueError(f"Invalid assay type: {assay_type}")
    
    logger.info(f"Linking {assay_type.upper()} experiment {experiment_id} to target {target_id}")
    
    with db_session() as db:
        try:
            if assay_type == "spr":
                experiment = db.query(SPRExperiment).filter(SPRExperiment.id == experiment_id).first()
            elif assay_type == "mst":
                experiment = db.query(MSTExperiment).filter(MSTExperiment.id == experiment_id).first()
            else:  # dsc
                experiment = db.query(DSCExperiment).filter(DSCExperiment.id == experiment_id).first()
            
            if not experiment:
                raise ValueError(f"{assay_type.upper()} experiment not found: {experiment_id}")
            
            experiment.target_id = target_id
            db.commit()
            
            logger.info(f"Successfully linked {assay_type.upper()} experiment to target")
            
        except Exception as e:
            logger.error(f"Failed to link experiment to target: {e}")
            db.rollback()
            raise


def get_compound_biophysical_profile(compound_id: UUID) -> Dict[str, List]:
    """
    Get all biophysical data for a compound across assay types.
    
    Args:
        compound_id: UUID of the compound
        
    Returns:
        Dictionary with keys 'spr', 'mst', 'dsc' containing lists of experiments
    """
    logger.info(f"Retrieving biophysical profile for compound: {compound_id}")
    
    with db_session() as db:
        profile = {
            'spr': [],
            'mst': [],
            'dsc': []
        }
        
        try:
            # Get SPR experiments
            spr_experiments = db.query(SPRExperiment).filter(
                SPRExperiment.compound_id == compound_id
            ).all()
            
            for exp in spr_experiments:
                db.expunge(exp)
                profile['spr'].append(exp)
            
            # Get MST experiments
            mst_experiments = db.query(MSTExperiment).filter(
                MSTExperiment.compound_id == compound_id
            ).all()
            
            for exp in mst_experiments:
                db.expunge(exp)
                profile['mst'].append(exp)
            
            # Get DSC experiments
            dsc_experiments = db.query(DSCExperiment).filter(
                DSCExperiment.compound_id == compound_id
            ).all()
            
            for exp in dsc_experiments:
                db.expunge(exp)
                profile['dsc'].append(exp)
            
            logger.info(f"Found {len(profile['spr'])} SPR, {len(profile['mst'])} MST, "
                       f"{len(profile['dsc'])} DSC experiments for compound")
            
            return profile
            
        except Exception as e:
            logger.error(f"Failed to retrieve compound profile: {e}")
            return profile


def get_processing_status(experiment_id: UUID, assay_type: str) -> str:
    """
    Get current processing status of an experiment.
    
    Args:
        experiment_id: UUID of the experiment
        assay_type: Type of assay ("spr", "mst", "dsc")
        
    Returns:
        Processing status string ("pending", "running", "completed", "failed")
        
    Raises:
        ValueError: If assay_type is invalid or experiment not found
    """
    if assay_type not in ["spr", "mst", "dsc"]:
        raise ValueError(f"Invalid assay type: {assay_type}")
    
    with db_session() as db:
        try:
            if assay_type == "spr":
                experiment = db.query(SPRExperiment).filter(SPRExperiment.id == experiment_id).first()
            elif assay_type == "mst":
                experiment = db.query(MSTExperiment).filter(MSTExperiment.id == experiment_id).first()
            else:  # dsc
                experiment = db.query(DSCExperiment).filter(DSCExperiment.id == experiment_id).first()
            
            if not experiment:
                raise ValueError(f"{assay_type.upper()} experiment not found: {experiment_id}")
            
            return experiment.processing_status
            
        except Exception as e:
            logger.error(f"Failed to get processing status: {e}")
            raise


def reprocess_experiment(experiment_id: UUID, assay_type: str) -> None:
    """
    Re-run analysis on existing experiment data.
    
    Args:
        experiment_id: UUID of the experiment to reprocess
        assay_type: Type of assay ("spr", "mst", "dsc")
        
    Raises:
        ValueError: If assay_type is invalid or experiment not found
    """
    if assay_type not in ["spr", "mst", "dsc"]:
        raise ValueError(f"Invalid assay type: {assay_type}")
    
    logger.info(f"Reprocessing {assay_type.upper()} experiment: {experiment_id}")
    
    with db_session() as db:
        try:
            # Check experiment exists and reset status
            if assay_type == "spr":
                experiment = db.query(SPRExperiment).filter(SPRExperiment.id == experiment_id).first()
            elif assay_type == "mst":
                experiment = db.query(MSTExperiment).filter(MSTExperiment.id == experiment_id).first()
            else:  # dsc
                experiment = db.query(DSCExperiment).filter(DSCExperiment.id == experiment_id).first()
            
            if not experiment:
                raise ValueError(f"{assay_type.upper()} experiment not found: {experiment_id}")
            
            # Reset processing status
            experiment.processing_status = "pending"
            experiment.error_message = None
            db.commit()
            
            # Start background processing
            if assay_type == "spr":
                processing_func = _process_spr_async
            elif assay_type == "mst":
                processing_func = _process_mst_async
            else:  # dsc
                processing_func = _process_dsc_async
            
            processing_thread = threading.Thread(
                target=processing_func,
                args=(experiment_id,),
                daemon=True
            )
            processing_thread.start()
            
            logger.info(f"Reprocessing started for {assay_type.upper()} experiment")
            
        except Exception as e:
            logger.error(f"Failed to start reprocessing: {e}")
            db.rollback()
            raise


# Update module exports
__all__ = [
    "ingest_spr_file",
    "ingest_mst_file", 
    "ingest_dsc_file",
    "link_to_compound",
    "link_to_target",
    "get_compound_biophysical_profile",
    "get_processing_status",
    "reprocess_experiment",
]
