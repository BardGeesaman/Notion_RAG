"""Feature extraction from segmented cells using scikit-image."""

from __future__ import annotations

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict

logger = logging.getLogger(__name__)


@dataclass
class CellMorphologyFeatures:
    """Morphological features of a single cell."""
    
    # Basic shape features
    area: float
    perimeter: float
    major_axis_length: float
    minor_axis_length: float
    eccentricity: float
    solidity: float
    extent: float
    
    # Position features
    centroid_x: float
    centroid_y: float
    
    # Additional shape descriptors
    circularity: float  # 4π * area / perimeter²
    aspect_ratio: float  # major_axis_length / minor_axis_length
    
    def to_dict(self) -> Dict[str, float]:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


@dataclass
class CellIntensityFeatures:
    """Intensity features of a single cell across channels."""
    
    channel_features: Dict[str, Dict[str, float]]  # channel -> feature_name -> value
    
    def to_dict(self) -> Dict[str, Dict[str, float]]:
        """Convert to dictionary for JSON serialization."""
        return self.channel_features


@dataclass
class WellAggregatedFeatures:
    """Aggregated features at the well level."""
    
    cell_count: int
    morphology_stats: Dict[str, Dict[str, float]]  # feature -> stat -> value
    intensity_stats: Dict[str, Dict[str, Dict[str, float]]]  # channel -> feature -> stat -> value
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "cell_count": self.cell_count,
            "morphology_stats": self.morphology_stats,
            "intensity_stats": self.intensity_stats
        }


class FeatureExtractor:
    """Extract morphological and intensity features from segmented cells."""
    
    def __init__(self):
        """Initialize feature extractor."""
        # Lazy import to avoid overhead
        self._skimage = None
        self._regionprops = None
        
        logger.info("Initialized FeatureExtractor")
    
    @property
    def skimage(self):
        """Lazy import of scikit-image."""
        if self._skimage is None:
            try:
                import skimage
                self._skimage = skimage
            except ImportError as e:
                raise ImportError(
                    "scikit-image is not installed. Install with: pip install scikit-image>=0.21"
                ) from e
        return self._skimage
    
    @property
    def regionprops(self):
        """Lazy import of regionprops."""
        if self._regionprops is None:
            from skimage.measure import regionprops
            self._regionprops = regionprops
        return self._regionprops
    
    def extract_morphology_features(self, masks: np.ndarray) -> List[CellMorphologyFeatures]:
        """
        Extract morphological features from segmentation masks.
        
        Args:
            masks: Segmentation masks with labeled cells (H, W)
        
        Returns:
            List of morphological features for each cell
        """
        if len(masks.shape) != 2:
            raise ValueError("Masks must be 2D array")
        
        # Get region properties
        regions = self.regionprops(masks)
        
        features = []
        for region in regions:
            # Calculate basic features
            area = region.area
            perimeter = region.perimeter
            major_axis = region.major_axis_length
            minor_axis = region.minor_axis_length
            eccentricity = region.eccentricity
            solidity = region.solidity
            extent = region.extent
            centroid_y, centroid_x = region.centroid
            
            # Calculate derived features
            circularity = 4 * np.pi * area / (perimeter ** 2) if perimeter > 0 else 0
            aspect_ratio = major_axis / minor_axis if minor_axis > 0 else 1
            
            cell_features = CellMorphologyFeatures(
                area=area,
                perimeter=perimeter,
                major_axis_length=major_axis,
                minor_axis_length=minor_axis,
                eccentricity=eccentricity,
                solidity=solidity,
                extent=extent,
                centroid_x=centroid_x,
                centroid_y=centroid_y,
                circularity=circularity,
                aspect_ratio=aspect_ratio
            )
            features.append(cell_features)
        
        logger.debug(f"Extracted morphology features for {len(features)} cells")
        return features
    
    def extract_intensity_features(self, masks: np.ndarray, 
                                 images: Dict[str, np.ndarray]) -> List[CellIntensityFeatures]:
        """
        Extract intensity features from multi-channel images.
        
        Args:
            masks: Segmentation masks with labeled cells (H, W)
            images: Dictionary of channel_name -> image arrays (H, W)
        
        Returns:
            List of intensity features for each cell
        """
        if len(masks.shape) != 2:
            raise ValueError("Masks must be 2D array")
        
        # Validate image shapes
        mask_shape = masks.shape
        for channel, image in images.items():
            if image.shape[:2] != mask_shape:
                raise ValueError(f"Image shape {image.shape} doesn't match mask shape {mask_shape}")
        
        # Get unique cell labels (excluding background 0)
        cell_labels = np.unique(masks)
        cell_labels = cell_labels[cell_labels > 0]
        
        features = []
        
        for label in cell_labels:
            # Get mask for this cell
            cell_mask = masks == label
            
            # Extract features for each channel
            channel_features = {}
            
            for channel_name, image in images.items():
                # Get pixel values for this cell
                cell_pixels = image[cell_mask]
                
                if len(cell_pixels) > 0:
                    # Calculate intensity statistics
                    intensity_stats = {
                        "mean": float(np.mean(cell_pixels)),
                        "median": float(np.median(cell_pixels)),
                        "std": float(np.std(cell_pixels)),
                        "min": float(np.min(cell_pixels)),
                        "max": float(np.max(cell_pixels)),
                        "sum": float(np.sum(cell_pixels)),
                        "pixel_count": len(cell_pixels)
                    }
                else:
                    # Handle edge case of empty cell
                    intensity_stats = {
                        "mean": 0.0, "median": 0.0, "std": 0.0,
                        "min": 0.0, "max": 0.0, "sum": 0.0,
                        "pixel_count": 0
                    }
                
                channel_features[channel_name] = intensity_stats
            
            cell_intensity = CellIntensityFeatures(channel_features=channel_features)
            features.append(cell_intensity)
        
        logger.debug(f"Extracted intensity features for {len(features)} cells across {len(images)} channels")
        return features
    
    def extract_all_features(self, masks: np.ndarray, 
                           images: Optional[Dict[str, np.ndarray]] = None) -> Tuple[List[CellMorphologyFeatures], List[CellIntensityFeatures]]:
        """
        Extract both morphological and intensity features.
        
        Args:
            masks: Segmentation masks with labeled cells
            images: Optional dictionary of channel images
        
        Returns:
            Tuple of (morphology_features, intensity_features)
        """
        # Extract morphological features
        morphology_features = self.extract_morphology_features(masks)
        
        # Extract intensity features if images provided
        if images is not None:
            intensity_features = self.extract_intensity_features(masks, images)
        else:
            intensity_features = []
        
        return morphology_features, intensity_features
    
    def aggregate_to_well(self, morphology_features: List[CellMorphologyFeatures],
                         intensity_features: List[CellIntensityFeatures]) -> WellAggregatedFeatures:
        """
        Aggregate cell-level features to well-level statistics.
        
        Args:
            morphology_features: List of morphological features for cells
            intensity_features: List of intensity features for cells
        
        Returns:
            Aggregated features at well level
        """
        cell_count = len(morphology_features)
        
        # Aggregate morphological features
        morphology_stats = self._aggregate_morphology(morphology_features)
        
        # Aggregate intensity features
        intensity_stats = self._aggregate_intensity(intensity_features)
        
        return WellAggregatedFeatures(
            cell_count=cell_count,
            morphology_stats=morphology_stats,
            intensity_stats=intensity_stats
        )
    
    def _aggregate_morphology(self, features: List[CellMorphologyFeatures]) -> Dict[str, Dict[str, float]]:
        """Aggregate morphological features across cells."""
        if not features:
            return {}
        
        # Get all feature names from first cell
        feature_names = list(features[0].to_dict().keys())
        
        aggregated = {}
        
        for feature_name in feature_names:
            # Collect values for this feature across all cells
            values = [getattr(cell, feature_name) for cell in features]
            
            # Calculate statistics
            stats = {
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "q25": float(np.percentile(values, 25)),
                "q75": float(np.percentile(values, 75))
            }
            
            aggregated[feature_name] = stats
        
        return aggregated
    
    def _aggregate_intensity(self, features: List[CellIntensityFeatures]) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Aggregate intensity features across cells."""
        if not features:
            return {}
        
        # Get all channels and intensity feature names
        first_cell = features[0]
        channels = list(first_cell.channel_features.keys())
        
        aggregated = {}
        
        for channel in channels:
            channel_aggregated = {}
            
            # Get intensity feature names for this channel
            intensity_feature_names = list(first_cell.channel_features[channel].keys())
            
            for intensity_feature in intensity_feature_names:
                # Collect values for this intensity feature across all cells
                values = []
                for cell in features:
                    if channel in cell.channel_features and intensity_feature in cell.channel_features[channel]:
                        values.append(cell.channel_features[channel][intensity_feature])
                
                if values:
                    # Calculate statistics
                    stats = {
                        "mean": float(np.mean(values)),
                        "median": float(np.median(values)),
                        "std": float(np.std(values)),
                        "min": float(np.min(values)),
                        "max": float(np.max(values)),
                        "q25": float(np.percentile(values, 25)),
                        "q75": float(np.percentile(values, 75))
                    }
                    
                    channel_aggregated[intensity_feature] = stats
            
            aggregated[channel] = channel_aggregated
        
        return aggregated
    
    def calculate_well_metrics(self, aggregated_features: WellAggregatedFeatures) -> Dict[str, float]:
        """
        Calculate high-level well metrics from aggregated features.
        
        Args:
            aggregated_features: Aggregated features for the well
        
        Returns:
            Dictionary of well-level metrics
        """
        metrics = {
            "cell_count": aggregated_features.cell_count
        }
        
        # Add morphological metrics
        morphology = aggregated_features.morphology_stats
        if "area" in morphology:
            metrics["mean_cell_area"] = morphology["area"]["mean"]
            metrics["cell_area_cv"] = morphology["area"]["std"] / morphology["area"]["mean"] if morphology["area"]["mean"] > 0 else 0
        
        if "circularity" in morphology:
            metrics["mean_circularity"] = morphology["circularity"]["mean"]
        
        if "eccentricity" in morphology:
            metrics["mean_eccentricity"] = morphology["eccentricity"]["mean"]
        
        # Add intensity metrics (for first channel if available)
        intensity = aggregated_features.intensity_stats
        if intensity:
            first_channel = list(intensity.keys())[0]
            if "mean" in intensity[first_channel]:
                metrics[f"{first_channel}_mean_intensity"] = intensity[first_channel]["mean"]["mean"]
                metrics[f"{first_channel}_intensity_cv"] = (
                    intensity[first_channel]["mean"]["std"] / intensity[first_channel]["mean"]["mean"]
                    if intensity[first_channel]["mean"]["mean"] > 0 else 0
                )
        
        return metrics
