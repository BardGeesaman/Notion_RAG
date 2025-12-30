"""Quality control metrics for HTS imaging assays."""

from __future__ import annotations

import logging
import numpy as np
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple
from uuid import UUID

from sqlalchemy.orm import Session
from sqlalchemy import and_, func

from amprenta_rag.models.chemistry import HTSPlate, HTSWell
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.imaging.aggregation import PlateAggregator
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class QCMetrics:
    """Quality control metrics for HTS imaging assays."""
    
    def __init__(self):
        """Initialize QC metrics calculator."""
        self.aggregator = PlateAggregator()
    
    def calculate_zprime(
        self,
        plate_id: UUID,
        feature_name: str,
        positive_control_positions: List[str],
        negative_control_positions: List[str],
        aggregation_method: str = 'mean',
        db: Session = None
    ) -> Dict[str, Any]:
        """
        Calculate Z' factor for assay quality assessment.
        
        Z' = 1 - (3 * (σp + σn)) / |μp - μn|
        where σp, σn are standard deviations and μp, μn are means
        of positive and negative controls.
        
        Z' > 0.5: Excellent assay
        0.5 > Z' > 0: Marginal assay  
        Z' < 0: Poor assay
        
        Args:
            plate_id: Plate UUID
            feature_name: Name of feature to analyze
            positive_control_positions: List of positive control well positions
            negative_control_positions: List of negative control well positions
            aggregation_method: Method to aggregate cell features ('mean', 'median')
            db: Database session
        
        Returns:
            Dictionary with Z' factor and control statistics
        """
        try:
            logger.info(f"Calculating Z' factor for feature {feature_name} in plate {plate_id}")
            
            # Get aggregated plate features
            df = self.aggregator.aggregate_plate_features(
                plate_id, [feature_name], [aggregation_method], db
            )
            
            if df.empty:
                raise ValueError(f"No feature data available for plate {plate_id}")
            
            feature_col = f'{feature_name}_{aggregation_method}'
            if feature_col not in df.columns:
                raise ValueError(f"Feature {feature_col} not found in aggregated data")
            
            # Get positive control values
            positive_data = df[df['well_position'].isin(positive_control_positions)][feature_col].dropna()
            if len(positive_data) == 0:
                raise ValueError(f"No data found for positive control positions: {positive_control_positions}")
            
            # Get negative control values
            negative_data = df[df['well_position'].isin(negative_control_positions)][feature_col].dropna()
            if len(negative_data) == 0:
                raise ValueError(f"No data found for negative control positions: {negative_control_positions}")
            
            # Calculate statistics
            pos_mean = float(positive_data.mean())
            pos_std = float(positive_data.std())
            neg_mean = float(negative_data.mean())
            neg_std = float(negative_data.std())
            
            # Calculate Z' factor
            mean_diff = abs(pos_mean - neg_mean)
            if mean_diff == 0:
                zprime = float('-inf')
                logger.warning("Zero difference between control means - Z' factor is undefined")
            else:
                zprime = 1 - (3 * (pos_std + neg_std)) / mean_diff
            
            # Interpret Z' factor
            if zprime > 0.5:
                quality = "Excellent"
            elif zprime > 0:
                quality = "Marginal"
            else:
                quality = "Poor"
            
            result = {
                "plate_id": str(plate_id),
                "feature_name": feature_name,
                "zprime_factor": float(zprime),
                "assay_quality": quality,
                "positive_controls": {
                    "positions": positive_control_positions,
                    "count": len(positive_data),
                    "mean": pos_mean,
                    "std": pos_std,
                    "cv": (pos_std / pos_mean * 100) if pos_mean != 0 else float('inf')
                },
                "negative_controls": {
                    "positions": negative_control_positions,
                    "count": len(negative_data),
                    "mean": neg_mean,
                    "std": neg_std,
                    "cv": (neg_std / neg_mean * 100) if neg_mean != 0 else float('inf')
                },
                "signal_window": mean_diff,
                "signal_to_background": pos_mean / neg_mean if neg_mean != 0 else float('inf')
            }
            
            logger.info(f"Z' factor calculated: {zprime:.3f} ({quality})")
            return result
            
        except Exception as e:
            logger.error(f"Failed to calculate Z' factor: {str(e)}")
            raise
    
    def calculate_zscore(
        self,
        plate_id: UUID,
        feature_name: str,
        aggregation_method: str = 'mean',
        control_positions: Optional[List[str]] = None,
        db: Session = None
    ) -> Dict[str, Any]:
        """
        Calculate Z-scores for all wells in a plate.
        
        Z-score = (value - μ) / σ
        where μ and σ are calculated from control wells or entire plate.
        
        Args:
            plate_id: Plate UUID
            feature_name: Name of feature to analyze
            aggregation_method: Method to aggregate cell features
            control_positions: List of control well positions (None to use all wells)
            db: Database session
        
        Returns:
            Dictionary with Z-scores and statistics
        """
        try:
            logger.info(f"Calculating Z-scores for feature {feature_name} in plate {plate_id}")
            
            # Get aggregated plate features
            df = self.aggregator.aggregate_plate_features(
                plate_id, [feature_name], [aggregation_method], db
            )
            
            if df.empty:
                raise ValueError(f"No feature data available for plate {plate_id}")
            
            feature_col = f'{feature_name}_{aggregation_method}'
            if feature_col not in df.columns:
                raise ValueError(f"Feature {feature_col} not found in aggregated data")
            
            # Calculate reference statistics (from controls or all wells)
            if control_positions:
                reference_data = df[df['well_position'].isin(control_positions)][feature_col].dropna()
                reference_type = "controls"
            else:
                reference_data = df[feature_col].dropna()
                reference_type = "all_wells"
            
            if len(reference_data) == 0:
                raise ValueError("No reference data available for Z-score calculation")
            
            ref_mean = float(reference_data.mean())
            ref_std = float(reference_data.std())
            
            if ref_std == 0:
                logger.warning("Zero standard deviation in reference data - Z-scores will be undefined")
                ref_std = 1.0  # Avoid division by zero
            
            # Calculate Z-scores for all wells
            df['zscore'] = (df[feature_col] - ref_mean) / ref_std
            
            # Identify outliers (|Z-score| > 2 or 3)
            outliers_2sigma = df[abs(df['zscore']) > 2]['well_position'].tolist()
            outliers_3sigma = df[abs(df['zscore']) > 3]['well_position'].tolist()
            
            # Create well-to-zscore mapping
            zscore_data = {}
            for _, row in df.iterrows():
                if pd.notna(row['zscore']):
                    zscore_data[row['well_position']] = float(row['zscore'])
            
            result = {
                "plate_id": str(plate_id),
                "feature_name": feature_name,
                "reference_type": reference_type,
                "reference_statistics": {
                    "mean": ref_mean,
                    "std": ref_std,
                    "count": len(reference_data)
                },
                "zscore_data": zscore_data,
                "outliers_2sigma": outliers_2sigma,
                "outliers_3sigma": outliers_3sigma,
                "zscore_statistics": {
                    "mean": float(df['zscore'].mean()),
                    "std": float(df['zscore'].std()),
                    "min": float(df['zscore'].min()),
                    "max": float(df['zscore'].max())
                }
            }
            
            logger.info(f"Calculated Z-scores for {len(zscore_data)} wells")
            return result
            
        except Exception as e:
            logger.error(f"Failed to calculate Z-scores: {str(e)}")
            raise
    
    def calculate_cv(
        self,
        well_id: UUID,
        feature_name: str,
        db: Session
    ) -> Dict[str, Any]:
        """
        Calculate coefficient of variation (CV) for a well.
        
        CV = (standard deviation / mean) * 100%
        
        Args:
            well_id: Well UUID
            feature_name: Name of feature to analyze
            db: Database session
        
        Returns:
            Dictionary with CV and statistics
        """
        try:
            logger.debug(f"Calculating CV for feature {feature_name} in well {well_id}")
            
            # Get all cell features for this well
            features = db.query(CellFeature).join(CellSegmentation).join(MicroscopyImage).filter(
                MicroscopyImage.well_id == well_id
            ).all()
            
            if not features:
                raise ValueError(f"No features found for well {well_id}")
            
            # Extract feature values
            values = []
            for feature in features:
                if hasattr(feature, feature_name):
                    value = getattr(feature, feature_name)
                    if value is not None:
                        values.append(float(value))
            
            if len(values) == 0:
                raise ValueError(f"No valid {feature_name} values found for well {well_id}")
            
            # Calculate statistics
            values = np.array(values)
            mean_val = float(np.mean(values))
            std_val = float(np.std(values))
            cv = (std_val / mean_val * 100) if mean_val != 0 else float('inf')
            
            result = {
                "well_id": str(well_id),
                "feature_name": feature_name,
                "cell_count": len(values),
                "mean": mean_val,
                "std": std_val,
                "cv_percent": cv,
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "median": float(np.median(values))
            }
            
            logger.debug(f"CV calculated: {cv:.2f}% for {len(values)} cells")
            return result
            
        except Exception as e:
            logger.error(f"Failed to calculate CV: {str(e)}")
            raise
    
    def calculate_plate_uniformity(
        self,
        plate_id: UUID,
        feature_name: str,
        aggregation_method: str = 'mean',
        db: Session = None
    ) -> Dict[str, Any]:
        """
        Calculate plate uniformity metrics.
        
        Args:
            plate_id: Plate UUID
            feature_name: Name of feature to analyze
            aggregation_method: Method to aggregate cell features
            db: Database session
        
        Returns:
            Dictionary with uniformity metrics
        """
        try:
            logger.info(f"Calculating plate uniformity for feature {feature_name}")
            
            # Get aggregated plate features
            df = self.aggregator.aggregate_plate_features(
                plate_id, [feature_name], [aggregation_method], db
            )
            
            if df.empty:
                raise ValueError(f"No feature data available for plate {plate_id}")
            
            feature_col = f'{feature_name}_{aggregation_method}'
            if feature_col not in df.columns:
                raise ValueError(f"Feature {feature_col} not found in aggregated data")
            
            values = df[feature_col].dropna()
            
            if len(values) == 0:
                raise ValueError("No valid feature values for uniformity calculation")
            
            # Calculate uniformity metrics
            mean_val = float(values.mean())
            std_val = float(values.std())
            cv = (std_val / mean_val * 100) if mean_val != 0 else float('inf')
            
            # Calculate edge effects (compare edge vs center wells)
            edge_positions, center_positions = self._classify_edge_center_wells(df['well_position'].tolist())
            
            edge_values = df[df['well_position'].isin(edge_positions)][feature_col].dropna()
            center_values = df[df['well_position'].isin(center_positions)][feature_col].dropna()
            
            edge_effect = None
            if len(edge_values) > 0 and len(center_values) > 0:
                edge_mean = float(edge_values.mean())
                center_mean = float(center_values.mean())
                edge_effect = {
                    "edge_mean": edge_mean,
                    "center_mean": center_mean,
                    "difference": edge_mean - center_mean,
                    "percent_difference": ((edge_mean - center_mean) / center_mean * 100) if center_mean != 0 else float('inf')
                }
            
            result = {
                "plate_id": str(plate_id),
                "feature_name": feature_name,
                "well_count": len(values),
                "uniformity_metrics": {
                    "mean": mean_val,
                    "std": std_val,
                    "cv_percent": cv,
                    "min": float(values.min()),
                    "max": float(values.max()),
                    "range": float(values.max() - values.min()),
                    "q25": float(values.quantile(0.25)),
                    "q75": float(values.quantile(0.75)),
                    "iqr": float(values.quantile(0.75) - values.quantile(0.25))
                },
                "edge_effect": edge_effect
            }
            
            logger.info(f"Plate uniformity CV: {cv:.2f}%")
            return result
            
        except Exception as e:
            logger.error(f"Failed to calculate plate uniformity: {str(e)}")
            raise
    
    def generate_qc_report(
        self,
        plate_id: UUID,
        feature_name: str,
        positive_control_positions: Optional[List[str]] = None,
        negative_control_positions: Optional[List[str]] = None,
        db: Session = None
    ) -> Dict[str, Any]:
        """
        Generate comprehensive QC report for a plate.
        
        Args:
            plate_id: Plate UUID
            feature_name: Name of feature to analyze
            positive_control_positions: Positive control well positions
            negative_control_positions: Negative control well positions
            db: Database session
        
        Returns:
            Comprehensive QC report dictionary
        """
        try:
            logger.info(f"Generating QC report for plate {plate_id}, feature {feature_name}")
            
            report = {
                "plate_id": str(plate_id),
                "feature_name": feature_name,
                "generated_at": pd.Timestamp.now().isoformat()
            }
            
            # Calculate uniformity metrics
            try:
                uniformity = self.calculate_plate_uniformity(plate_id, feature_name, db=db)
                report["uniformity"] = uniformity["uniformity_metrics"]
                report["edge_effect"] = uniformity["edge_effect"]
            except Exception as e:
                logger.warning(f"Failed to calculate uniformity: {str(e)}")
                report["uniformity"] = None
                report["edge_effect"] = None
            
            # Calculate Z' factor if controls are provided
            if positive_control_positions and negative_control_positions:
                try:
                    zprime = self.calculate_zprime(
                        plate_id, feature_name, positive_control_positions, negative_control_positions, db=db
                    )
                    report["zprime"] = zprime
                except Exception as e:
                    logger.warning(f"Failed to calculate Z' factor: {str(e)}")
                    report["zprime"] = None
            else:
                report["zprime"] = None
            
            # Calculate Z-scores
            try:
                zscore_result = self.calculate_zscore(plate_id, feature_name, db=db)
                report["zscore_statistics"] = zscore_result["zscore_statistics"]
                report["outliers_2sigma"] = zscore_result["outliers_2sigma"]
                report["outliers_3sigma"] = zscore_result["outliers_3sigma"]
            except Exception as e:
                logger.warning(f"Failed to calculate Z-scores: {str(e)}")
                report["zscore_statistics"] = None
                report["outliers_2sigma"] = []
                report["outliers_3sigma"] = []
            
            # Overall QC assessment
            qc_flags = []
            if report["uniformity"] and report["uniformity"]["cv_percent"] > 20:
                qc_flags.append("High CV (>20%)")
            
            if report["zprime"] and report["zprime"]["zprime_factor"] < 0:
                qc_flags.append("Poor Z' factor (<0)")
            
            if len(report.get("outliers_3sigma", [])) > len(report.get("outliers_2sigma", [])) * 0.1:
                qc_flags.append("High outlier rate")
            
            report["qc_flags"] = qc_flags
            report["qc_status"] = "PASS" if len(qc_flags) == 0 else "WARNING"
            
            logger.info(f"QC report generated with {len(qc_flags)} flags")
            return report
            
        except Exception as e:
            logger.error(f"Failed to generate QC report: {str(e)}")
            raise
    
    def _classify_edge_center_wells(self, well_positions: List[str]) -> Tuple[List[str], List[str]]:
        """
        Classify wells as edge or center wells.
        
        Args:
            well_positions: List of well position strings
        
        Returns:
            Tuple of (edge_wells, center_wells)
        """
        edge_wells = []
        center_wells = []
        
        # Parse well positions to get row and column ranges
        rows = set()
        cols = set()
        
        for pos in well_positions:
            try:
                # Extract row letter and column number
                row_match = ''.join([c for c in pos if c.isalpha()])
                col_match = ''.join([c for c in pos if c.isdigit()])
                
                if row_match and col_match:
                    rows.add(row_match)
                    cols.add(int(col_match))
            except:
                continue
        
        if not rows or not cols:
            return well_positions, []  # Return all as edge if parsing fails
        
        min_row = min(rows)
        max_row = max(rows)
        min_col = min(cols)
        max_col = max(cols)
        
        # Classify wells
        for pos in well_positions:
            try:
                row_match = ''.join([c for c in pos if c.isalpha()])
                col_match = int(''.join([c for c in pos if c.isdigit()]))
                
                # Edge wells are on the perimeter
                if (row_match == min_row or row_match == max_row or 
                    col_match == min_col or col_match == max_col):
                    edge_wells.append(pos)
                else:
                    center_wells.append(pos)
            except:
                edge_wells.append(pos)  # Default to edge if parsing fails
        
        return edge_wells, center_wells
