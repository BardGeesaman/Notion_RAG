"""CellPose integration service with GPU memory management and tiling support."""

from __future__ import annotations

import logging
import numpy as np
from typing import Optional, Tuple, Dict, Any, List
from pathlib import Path

from .tiling import TileManager, Tile

logger = logging.getLogger(__name__)

# Constants
MAX_IMAGE_SIZE = 2048  # Auto-tile images larger than this
DEFAULT_DIAMETER = 30  # Default cell diameter in pixels
GPU_MEMORY_LIMIT_MB = 8192  # Conservative GPU memory limit


class CellPoseService:
    """
    CellPose wrapper service with automatic tiling and GPU memory management.
    """
    
    def __init__(self, model_type: str = "cyto", gpu: bool = True, 
                 tile_size: int = 1024, overlap: int = 128):
        """
        Initialize CellPose service.
        
        Args:
            model_type: CellPose model type ('cyto', 'nuclei', 'cyto2')
            gpu: Whether to use GPU acceleration
            tile_size: Size of tiles for large images
            overlap: Overlap between tiles in pixels
        """
        self.model_type = model_type
        self.use_gpu = gpu
        self.tile_manager = TileManager(tile_size=tile_size, overlap=overlap)
        
        # Lazy import to avoid overhead if not used
        self._cellpose = None
        self._model = None
        
        logger.info(f"Initialized CellPose service: model={model_type}, gpu={gpu}")
    
    @property
    def cellpose(self):
        """Lazy import of cellpose module."""
        if self._cellpose is None:
            try:
                import cellpose
                self._cellpose = cellpose
            except ImportError as e:
                raise ImportError(
                    "CellPose is not installed. Install with: pip install cellpose>=3.0"
                ) from e
        return self._cellpose
    
    @property
    def model(self):
        """Lazy initialization of CellPose model."""
        if self._model is None:
            try:
                from cellpose import models
                self._model = models.Cellpose(
                    model_type=self.model_type,
                    gpu=self.use_gpu
                )
                logger.info(f"Loaded CellPose model: {self.model_type}")
            except Exception as e:
                logger.error(f"Failed to load CellPose model: {e}")
                raise
        return self._model
    
    def segment(self, image: np.ndarray, diameter: Optional[float] = None,
                channels: Optional[List[int]] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Segment cells in image with automatic tiling for large images.
        
        Args:
            image: Input image (H, W) or (H, W, C)
            diameter: Expected cell diameter in pixels
            channels: Channel configuration [cytoplasm, nucleus] or [0, 0] for grayscale
        
        Returns:
            Tuple of (masks, flows) where masks are labeled segmentation masks
        """
        if diameter is None:
            diameter = DEFAULT_DIAMETER
        
        if channels is None:
            # Default to grayscale
            channels = [0, 0]
        
        # Validate input
        if len(image.shape) < 2:
            raise ValueError("Image must be at least 2D")
        
        height, width = image.shape[:2]
        
        # Check if tiling is needed
        if max(height, width) > MAX_IMAGE_SIZE:
            logger.info(f"Large image ({height}x{width}), using tiling approach")
            return self._segment_tiled(image, diameter, channels)
        else:
            # Check GPU memory requirements
            estimated_memory = self.estimate_memory_mb(image.shape)
            if estimated_memory > GPU_MEMORY_LIMIT_MB and self.use_gpu:
                logger.warning(f"Estimated memory ({estimated_memory}MB) exceeds limit, using CPU fallback")
                return self._segment_cpu_fallback(image, diameter, channels)
            
            logger.info(f"Segmenting image ({height}x{width}) directly")
            return self._segment_direct(image, diameter, channels)
    
    def _segment_direct(self, image: np.ndarray, diameter: float,
                       channels: List[int]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Segment image directly without tiling.
        
        Args:
            image: Input image
            diameter: Expected cell diameter
            channels: Channel configuration
        
        Returns:
            Tuple of (masks, flows)
        """
        try:
            masks, flows, styles, diams = self.model.eval(
                image,
                diameter=diameter,
                channels=channels
            )
            return masks, flows
        except Exception as e:
            if "out of memory" in str(e).lower() and self.use_gpu:
                logger.warning("GPU out of memory, falling back to CPU")
                return self._segment_cpu_fallback(image, diameter, channels)
            else:
                raise
    
    def _segment_cpu_fallback(self, image: np.ndarray, diameter: float,
                            channels: List[int]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fallback to CPU segmentation when GPU runs out of memory.
        
        Args:
            image: Input image
            diameter: Expected cell diameter
            channels: Channel configuration
        
        Returns:
            Tuple of (masks, flows)
        """
        logger.info("Using CPU fallback for segmentation")
        
        # Create CPU-only model
        from cellpose import models
        cpu_model = models.Cellpose(
            model_type=self.model_type,
            gpu=False
        )
        
        try:
            masks, flows, styles, diams = cpu_model.eval(
                image,
                diameter=diameter,
                channels=channels
            )
            return masks, flows
        except Exception as e:
            logger.error(f"CPU fallback also failed: {e}")
            raise
    
    def _segment_tiled(self, image: np.ndarray, diameter: float,
                      channels: List[int]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Segment large image using tiling approach.
        
        Args:
            image: Input image
            diameter: Expected cell diameter
            channels: Channel configuration
        
        Returns:
            Tuple of (masks, flows)
        """
        height, width = image.shape[:2]
        
        # Split image into tiles
        tiles = self.tile_manager.split_image(image)
        logger.info(f"Split image into {len(tiles)} tiles")
        
        # Segment each tile
        mask_tiles = []
        flow_tiles = []
        
        for i, tile in enumerate(tiles):
            try:
                logger.debug(f"Segmenting tile {i+1}/{len(tiles)}")
                tile_masks, tile_flows = self._segment_direct(tile.image, diameter, channels)
                mask_tiles.append((tile_masks, tile))
                flow_tiles.append((tile_flows, tile))
            except Exception as e:
                logger.error(f"Failed to segment tile {i}: {e}")
                # Create empty mask for failed tile
                empty_mask = np.zeros(tile.image.shape[:2], dtype=np.int32)
                empty_flow = np.zeros(tile.image.shape[:2] + (2,), dtype=np.float32)
                mask_tiles.append((empty_mask, tile))
                flow_tiles.append((empty_flow, tile))
        
        # Stitch masks back together
        stitched_masks = self.tile_manager.stitch_masks(mask_tiles, (height, width))
        
        # For flows, we'll create a simple average (flows are less critical for most applications)
        stitched_flows = self._stitch_flows(flow_tiles, (height, width))
        
        logger.info(f"Stitched {len(tiles)} tiles into final mask")
        
        return stitched_masks, stitched_flows
    
    def _stitch_flows(self, flow_tiles: List[Tuple[np.ndarray, Tile]], 
                     original_shape: Tuple[int, int]) -> np.ndarray:
        """
        Stitch flow fields from tiles (simplified approach).
        
        Args:
            flow_tiles: List of (flow, tile) pairs
            original_shape: Shape of original image
        
        Returns:
            Stitched flow field
        """
        height, width = original_shape
        
        # Determine flow shape (usually H, W, 2)
        if flow_tiles and len(flow_tiles[0][0].shape) == 3:
            flow_channels = flow_tiles[0][0].shape[2]
            stitched_flows = np.zeros((height, width, flow_channels), dtype=np.float32)
        else:
            stitched_flows = np.zeros((height, width), dtype=np.float32)
        
        # Simple averaging approach for overlapping regions
        weight_map = np.zeros((height, width), dtype=np.float32)
        
        for flows, tile in flow_tiles:
            y_start, y_end = tile.y_start, tile.y_end
            x_start, x_end = tile.x_start, tile.x_end
            
            # Adjust to image bounds
            y_start = max(0, y_start)
            y_end = min(height, y_end)
            x_start = max(0, x_start)
            x_end = min(width, x_end)
            
            # Calculate corresponding region in tile
            tile_y_start = y_start - tile.y_start
            tile_y_end = tile_y_start + (y_end - y_start)
            tile_x_start = x_start - tile.x_start
            tile_x_end = tile_x_start + (x_end - x_start)
            
            # Extract region from tile flows
            if len(flows.shape) == 3:
                tile_flows_region = flows[tile_y_start:tile_y_end, tile_x_start:tile_x_end, :]
                stitched_flows[y_start:y_end, x_start:x_end, :] += tile_flows_region
            else:
                tile_flows_region = flows[tile_y_start:tile_y_end, tile_x_start:tile_x_end]
                stitched_flows[y_start:y_end, x_start:x_end] += tile_flows_region
            
            # Update weight map
            weight_map[y_start:y_end, x_start:x_end] += 1.0
        
        # Normalize by weights to get average
        weight_map[weight_map == 0] = 1  # Avoid division by zero
        if len(stitched_flows.shape) == 3:
            stitched_flows = stitched_flows / weight_map[:, :, np.newaxis]
        else:
            stitched_flows = stitched_flows / weight_map
        
        return stitched_flows
    
    def estimate_memory_mb(self, image_shape: Tuple[int, ...]) -> float:
        """
        Estimate GPU memory requirements for segmentation.
        
        Args:
            image_shape: Shape of input image
        
        Returns:
            Estimated memory in MB
        """
        # Basic estimation based on image size
        # CellPose typically needs ~4x the image size in memory for intermediate calculations
        
        if len(image_shape) < 2:
            return 0.0
        
        height, width = image_shape[:2]
        channels = image_shape[2] if len(image_shape) > 2 else 1
        
        # Estimate memory usage
        # - Input image: height * width * channels * 4 bytes (float32)
        # - Intermediate tensors: ~3x input size
        # - Model parameters: ~50MB (rough estimate)
        
        input_size_mb = (height * width * channels * 4) / (1024 * 1024)
        intermediate_mb = input_size_mb * 3
        model_mb = 50
        
        total_mb = input_size_mb + intermediate_mb + model_mb
        
        return total_mb
    
    def get_model_info(self) -> Dict[str, Any]:
        """
        Get information about the loaded model.
        
        Returns:
            Dictionary with model information
        """
        info = {
            "model_type": self.model_type,
            "use_gpu": self.use_gpu,
            "tile_size": self.tile_manager.tile_size,
            "overlap": self.tile_manager.overlap,
            "max_image_size": MAX_IMAGE_SIZE,
        }
        
        if self._model is not None:
            info["model_loaded"] = True
            # Add any additional model-specific info
        else:
            info["model_loaded"] = False
        
        return info
    
    def count_cells(self, masks: np.ndarray) -> int:
        """
        Count number of cells in segmentation masks.
        
        Args:
            masks: Segmentation masks with labeled cells
        
        Returns:
            Number of unique cells (excluding background)
        """
        unique_labels = np.unique(masks)
        # Exclude background (label 0)
        cell_labels = unique_labels[unique_labels > 0]
        return len(cell_labels)
    
    def validate_image(self, image: np.ndarray) -> bool:
        """
        Validate input image for segmentation.
        
        Args:
            image: Input image to validate
        
        Returns:
            True if image is valid for segmentation
        """
        if not isinstance(image, np.ndarray):
            return False
        
        if len(image.shape) < 2:
            return False
        
        if image.shape[0] < 10 or image.shape[1] < 10:
            return False
        
        # Check for reasonable data type
        if image.dtype not in [np.uint8, np.uint16, np.float32, np.float64]:
            return False
        
        return True
