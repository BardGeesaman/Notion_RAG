"""Feature aggregation for HTS plates and wells."""

from __future__ import annotations

import logging
import numpy as np
import pandas as pd
from typing import List, Dict, Any, Optional, Union
from uuid import UUID

from sqlalchemy.orm import Session
from sqlalchemy import and_, func

from amprenta_rag.models.chemistry import HTSPlate, HTSWell
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.imaging.feature_extraction import CellMorphologyFeatures, WellAggregatedFeatures
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class PlateAggregator:
    """Aggregate cell features at plate and well levels."""
    
    def __init__(self):
        """Initialize plate aggregator."""
        pass
    
    def aggregate_plate_features(
        self,
        plate_id: UUID,
        feature_names: Optional[List[str]] = None,
        aggregation_methods: Optional[List[str]] = None,
        db: Session = None
    ) -> pd.DataFrame:
        """
        Aggregate cell features across an entire plate.
        
        Args:
            plate_id: Plate UUID
            feature_names: List of features to aggregate (None for all)
            aggregation_methods: List of methods ['mean', 'median', 'std', 'count']
            db: Database session
        
        Returns:
            DataFrame with aggregated features per well
        """
        try:
            if aggregation_methods is None:
                aggregation_methods = ['mean', 'median', 'std', 'count']
            
            logger.info(f"Aggregating plate features for plate {plate_id}")
            
            # Get all wells for the plate
            wells = db.query(HTSWell).filter(HTSWell.plate_id == plate_id).all()
            
            if not wells:
                logger.warning(f"No wells found for plate {plate_id}")
                return pd.DataFrame()
            
            # Collect aggregated features for each well
            well_data = []
            
            for well in wells:
                well_features = self.aggregate_well_features(
                    well.id, feature_names, aggregation_methods, db
                )
                
                if well_features:
                    # Add well metadata
                    well_features.update({
                        'well_id': str(well.id),
                        'well_position': well.position,
                        'compound_id': str(well.compound_id) if well.compound_id else None,
                        'concentration': well.concentration
                    })
                    well_data.append(well_features)
            
            if not well_data:
                logger.warning(f"No feature data found for plate {plate_id}")
                return pd.DataFrame()
            
            # Convert to DataFrame
            df = pd.DataFrame(well_data)
            
            # Sort by well position
            df = df.sort_values('well_position').reset_index(drop=True)
            
            logger.info(f"Aggregated features for {len(df)} wells in plate {plate_id}")
            return df
            
        except Exception as e:
            logger.error(f"Failed to aggregate plate features: {str(e)}")
            raise
    
    def aggregate_well_features(
        self,
        well_id: UUID,
        feature_names: Optional[List[str]] = None,
        aggregation_methods: Optional[List[str]] = None,
        db: Session = None
    ) -> Dict[str, Any]:
        """
        Aggregate cell features for a single well.
        
        Args:
            well_id: Well UUID
            feature_names: List of features to aggregate (None for all)
            aggregation_methods: List of methods ['mean', 'median', 'std', 'count']
            db: Database session
        
        Returns:
            Dictionary with aggregated features
        """
        try:
            if aggregation_methods is None:
                aggregation_methods = ['mean', 'median', 'std', 'count']
            
            # Get all cell features for this well
            features = db.query(CellFeature).join(CellSegmentation).join(MicroscopyImage).filter(
                MicroscopyImage.well_id == well_id
            ).all()
            
            if not features:
                logger.debug(f"No features found for well {well_id}")
                return {}
            
            # Convert to DataFrame for easier aggregation
            feature_data = []
            for feature in features:
                feature_dict = {
                    'cell_id': feature.cell_id,
                    'area': feature.area,
                    'perimeter': feature.perimeter,
                    'circularity': feature.circularity,
                    'eccentricity': feature.eccentricity,
                    'solidity': feature.solidity,
                    'centroid_x': feature.centroid_x,
                    'centroid_y': feature.centroid_y
                }
                
                # Add intensity features if available
                if feature.intensity_features:
                    for channel, intensities in feature.intensity_features.items():
                        if isinstance(intensities, dict):
                            for stat, value in intensities.items():
                                feature_dict[f'{channel}_{stat}'] = value
                
                # Add custom features if available
                if feature.custom_features:
                    for key, value in feature.custom_features.items():
                        if isinstance(value, (int, float)):
                            feature_dict[f'custom_{key}'] = value
                
                feature_data.append(feature_dict)
            
            df = pd.DataFrame(feature_data)
            
            # Filter features if specified
            if feature_names:
                available_features = [col for col in feature_names if col in df.columns]
                if available_features:
                    df = df[available_features]
                else:
                    logger.warning(f"None of the requested features {feature_names} found in data")
                    return {}
            
            # Calculate aggregations
            aggregated = {}
            
            for column in df.select_dtypes(include=[np.number]).columns:
                if column == 'cell_id':
                    continue
                
                series = df[column].dropna()
                if len(series) == 0:
                    continue
                
                for method in aggregation_methods:
                    if method == 'mean':
                        aggregated[f'{column}_mean'] = float(series.mean())
                    elif method == 'median':
                        aggregated[f'{column}_median'] = float(series.median())
                    elif method == 'std':
                        aggregated[f'{column}_std'] = float(series.std())
                    elif method == 'count':
                        aggregated[f'{column}_count'] = int(len(series))
                    elif method == 'min':
                        aggregated[f'{column}_min'] = float(series.min())
                    elif method == 'max':
                        aggregated[f'{column}_max'] = float(series.max())
                    elif method == 'q25':
                        aggregated[f'{column}_q25'] = float(series.quantile(0.25))
                    elif method == 'q75':
                        aggregated[f'{column}_q75'] = float(series.quantile(0.75))
            
            # Add cell count
            aggregated['total_cell_count'] = len(features)
            
            logger.debug(f"Aggregated {len(aggregated)} features for well {well_id}")
            return aggregated
            
        except Exception as e:
            logger.error(f"Failed to aggregate well features: {str(e)}")
            raise
    
    def get_plate_heatmap_data(
        self,
        plate_id: UUID,
        feature_name: str,
        aggregation_method: str = 'mean',
        db: Session = None
    ) -> Dict[str, float]:
        """
        Get plate heatmap data for a specific feature.
        
        Args:
            plate_id: Plate UUID
            feature_name: Name of feature to visualize
            aggregation_method: Aggregation method ('mean', 'median', etc.)
            db: Database session
        
        Returns:
            Dictionary mapping well positions to feature values
        """
        try:
            logger.info(f"Generating heatmap data for feature {feature_name} in plate {plate_id}")
            
            # Get aggregated plate features
            df = self.aggregate_plate_features(
                plate_id, [feature_name], [aggregation_method], db
            )
            
            if df.empty:
                logger.warning(f"No data available for feature {feature_name}")
                return {}
            
            # Create feature column name
            feature_col = f'{feature_name}_{aggregation_method}'
            
            if feature_col not in df.columns:
                logger.warning(f"Feature column {feature_col} not found in aggregated data")
                return {}
            
            # Create heatmap dictionary
            heatmap_data = {}
            for _, row in df.iterrows():
                well_position = row['well_position']
                feature_value = row[feature_col]
                
                if pd.notna(feature_value):
                    heatmap_data[well_position] = float(feature_value)
            
            logger.info(f"Generated heatmap data for {len(heatmap_data)} wells")
            return heatmap_data
            
        except Exception as e:
            logger.error(f"Failed to generate heatmap data: {str(e)}")
            raise
    
    def get_plate_feature_summary(
        self,
        plate_id: UUID,
        db: Session
    ) -> Dict[str, Any]:
        """
        Get summary statistics for all features in a plate.
        
        Args:
            plate_id: Plate UUID
            db: Database session
        
        Returns:
            Dictionary with feature summary statistics
        """
        try:
            # Get all features aggregated
            df = self.aggregate_plate_features(plate_id, db=db)
            
            if df.empty:
                return {
                    "plate_id": str(plate_id),
                    "well_count": 0,
                    "feature_count": 0,
                    "features": {}
                }
            
            # Calculate summary for each numeric feature
            feature_summaries = {}
            numeric_columns = df.select_dtypes(include=[np.number]).columns
            
            for column in numeric_columns:
                if column in ['well_id', 'concentration']:
                    continue
                
                series = df[column].dropna()
                if len(series) > 0:
                    feature_summaries[column] = {
                        'count': int(len(series)),
                        'mean': float(series.mean()),
                        'std': float(series.std()),
                        'min': float(series.min()),
                        'max': float(series.max()),
                        'q25': float(series.quantile(0.25)),
                        'q50': float(series.quantile(0.50)),
                        'q75': float(series.quantile(0.75))
                    }
            
            return {
                "plate_id": str(plate_id),
                "well_count": len(df),
                "feature_count": len(feature_summaries),
                "features": feature_summaries
            }
            
        except Exception as e:
            logger.error(f"Failed to get plate feature summary: {str(e)}")
            raise
    
    def compare_wells(
        self,
        well_ids: List[UUID],
        feature_names: Optional[List[str]] = None,
        db: Session = None
    ) -> pd.DataFrame:
        """
        Compare features across multiple wells.
        
        Args:
            well_ids: List of well UUIDs to compare
            feature_names: List of features to include (None for all)
            db: Database session
        
        Returns:
            DataFrame with features for each well
        """
        try:
            logger.info(f"Comparing features across {len(well_ids)} wells")
            
            well_data = []
            for well_id in well_ids:
                well_features = self.aggregate_well_features(well_id, feature_names, db=db)
                if well_features:
                    well_features['well_id'] = str(well_id)
                    
                    # Get well info
                    well = db.query(HTSWell).filter(HTSWell.id == well_id).first()
                    if well:
                        well_features['well_position'] = well.position
                        well_features['compound_id'] = str(well.compound_id) if well.compound_id else None
                        well_features['concentration'] = well.concentration
                    
                    well_data.append(well_features)
            
            if not well_data:
                logger.warning("No feature data found for any wells")
                return pd.DataFrame()
            
            df = pd.DataFrame(well_data)
            
            # Sort by well position if available
            if 'well_position' in df.columns:
                df = df.sort_values('well_position').reset_index(drop=True)
            
            logger.info(f"Compared features for {len(df)} wells")
            return df
            
        except Exception as e:
            logger.error(f"Failed to compare wells: {str(e)}")
            raise
    
    def get_feature_correlation_matrix(
        self,
        plate_id: UUID,
        feature_names: Optional[List[str]] = None,
        db: Session = None
    ) -> pd.DataFrame:
        """
        Calculate correlation matrix between features.
        
        Args:
            plate_id: Plate UUID
            feature_names: List of features to include (None for all numeric)
            db: Database session
        
        Returns:
            Correlation matrix as DataFrame
        """
        try:
            # Get plate features
            df = self.aggregate_plate_features(plate_id, feature_names, db=db)
            
            if df.empty:
                return pd.DataFrame()
            
            # Select only numeric columns
            numeric_df = df.select_dtypes(include=[np.number])
            
            # Remove ID columns
            id_columns = [col for col in numeric_df.columns if 'id' in col.lower()]
            numeric_df = numeric_df.drop(columns=id_columns, errors='ignore')
            
            if numeric_df.empty:
                return pd.DataFrame()
            
            # Calculate correlation matrix
            correlation_matrix = numeric_df.corr()
            
            logger.info(f"Calculated correlation matrix for {len(correlation_matrix)} features")
            return correlation_matrix
            
        except Exception as e:
            logger.error(f"Failed to calculate correlation matrix: {str(e)}")
            raise
