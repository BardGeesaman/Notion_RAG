"""Image tiling utilities for processing large microscopy images."""

from __future__ import annotations

import numpy as np
from typing import List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class Tile:
    """Represents a tile of a larger image."""
    
    image: np.ndarray
    x_start: int
    y_start: int
    x_end: int
    y_end: int
    tile_id: int
    
    @property
    def width(self) -> int:
        """Width of the tile."""
        return self.x_end - self.x_start
    
    @property
    def height(self) -> int:
        """Height of the tile."""
        return self.y_end - self.y_start
    
    @property
    def shape(self) -> Tuple[int, int]:
        """Shape of the tile (height, width)."""
        return (self.height, self.width)


class TileManager:
    """Manages tiling of large images for processing."""
    
    def __init__(self, tile_size: int = 1024, overlap: int = 128):
        """
        Initialize TileManager.
        
        Args:
            tile_size: Size of each tile (square)
            overlap: Overlap between adjacent tiles in pixels
        """
        self.tile_size = tile_size
        self.overlap = overlap
        self.effective_tile_size = tile_size - overlap
    
    def split_image(self, image: np.ndarray) -> List[Tile]:
        """
        Split large image into overlapping tiles.
        
        Args:
            image: Input image as numpy array (H, W) or (H, W, C)
        
        Returns:
            List of Tile objects
        """
        if len(image.shape) < 2:
            raise ValueError("Image must be at least 2D")
        
        height, width = image.shape[:2]
        tiles = []
        tile_id = 0
        
        # Calculate tile positions
        y_positions = self._calculate_positions(height, self.tile_size, self.overlap)
        x_positions = self._calculate_positions(width, self.tile_size, self.overlap)
        
        for y_start, y_end in y_positions:
            for x_start, x_end in x_positions:
                # Extract tile from image
                if len(image.shape) == 2:
                    tile_image = image[y_start:y_end, x_start:x_end]
                else:
                    tile_image = image[y_start:y_end, x_start:x_end, :]
                
                tile = Tile(
                    image=tile_image,
                    x_start=x_start,
                    y_start=y_start,
                    x_end=x_end,
                    y_end=y_end,
                    tile_id=tile_id
                )
                tiles.append(tile)
                tile_id += 1
        
        return tiles
    
    def stitch_masks(self, mask_tiles: List[Tuple[np.ndarray, Tile]], 
                     original_shape: Tuple[int, int]) -> np.ndarray:
        """
        Stitch segmentation masks from tiles back into full image.
        
        Args:
            mask_tiles: List of (mask, tile) pairs
            original_shape: Shape of original image (height, width)
        
        Returns:
            Stitched mask with globally unique cell IDs
        """
        height, width = original_shape
        stitched_mask = np.zeros((height, width), dtype=np.int32)
        next_cell_id = 1
        
        # Track which cells have been assigned in overlapping regions
        cell_mapping = {}
        
        for mask, tile in mask_tiles:
            # Get the region where this tile will be placed
            y_start, y_end = tile.y_start, tile.y_end
            x_start, x_end = tile.x_start, tile.x_end
            
            # Adjust mask to fit within original image bounds
            mask_y_start = max(0, y_start)
            mask_y_end = min(height, y_end)
            mask_x_start = max(0, x_start)
            mask_x_end = min(width, x_end)
            
            # Extract corresponding region from tile mask
            tile_y_start = mask_y_start - y_start
            tile_y_end = tile_y_start + (mask_y_end - mask_y_start)
            tile_x_start = mask_x_start - x_start
            tile_x_end = tile_x_start + (mask_x_end - mask_x_start)
            
            tile_mask_region = mask[tile_y_start:tile_y_end, tile_x_start:tile_x_end]
            
            # Relabel cells in this tile to avoid ID conflicts
            relabeled_mask = self._relabel_mask_region(
                tile_mask_region,
                stitched_mask[mask_y_start:mask_y_end, mask_x_start:mask_x_end],
                next_cell_id
            )
            
            # Update stitched mask
            stitched_mask[mask_y_start:mask_y_end, mask_x_start:mask_x_end] = relabeled_mask
            
            # Update next available cell ID
            if relabeled_mask.size > 0:
                max_id = np.max(relabeled_mask)
                if max_id > 0:
                    next_cell_id = max_id + 1
        
        return stitched_mask
    
    def _calculate_positions(self, dimension: int, tile_size: int, overlap: int) -> List[Tuple[int, int]]:
        """
        Calculate start/end positions for tiles along one dimension.
        
        Args:
            dimension: Size of the dimension (height or width)
            tile_size: Size of each tile
            overlap: Overlap between tiles
        
        Returns:
            List of (start, end) positions
        """
        positions = []
        effective_size = tile_size - overlap
        
        start = 0
        while start < dimension:
            end = min(start + tile_size, dimension)
            positions.append((start, end))
            
            # If this tile covers the rest of the dimension, we're done
            if end >= dimension:
                break
            
            # Move to next position
            start += effective_size
        
        return positions
    
    def _relabel_mask_region(self, new_mask: np.ndarray, existing_mask: np.ndarray, 
                           next_id: int) -> np.ndarray:
        """
        Relabel mask region to avoid conflicts with existing mask.
        
        Args:
            new_mask: New mask region to be added
            existing_mask: Existing mask in the same region
            next_id: Next available cell ID
        
        Returns:
            Relabeled mask that can be safely merged
        """
        if new_mask.shape != existing_mask.shape:
            raise ValueError("Mask shapes must match")
        
        result_mask = existing_mask.copy()
        
        # Get unique cell IDs in new mask (excluding background 0)
        unique_ids = np.unique(new_mask)
        unique_ids = unique_ids[unique_ids > 0]
        
        current_id = next_id
        
        for old_id in unique_ids:
            cell_pixels = new_mask == old_id
            
            # Check if this overlaps with existing cells
            overlapping_existing = existing_mask[cell_pixels]
            unique_overlaps = np.unique(overlapping_existing)
            
            if len(unique_overlaps) == 1 and unique_overlaps[0] == 0:
                # No overlap with existing cells, assign new ID
                result_mask[cell_pixels] = current_id
                current_id += 1
            else:
                # Overlap detected - use majority vote or skip
                # For simplicity, we'll keep the existing cell
                # In a more sophisticated approach, we could merge or split cells
                pass
        
        return result_mask
    
    def estimate_tile_count(self, image_shape: Tuple[int, int]) -> int:
        """
        Estimate number of tiles needed for an image.
        
        Args:
            image_shape: Shape of image (height, width)
        
        Returns:
            Estimated number of tiles
        """
        height, width = image_shape
        
        y_positions = self._calculate_positions(height, self.tile_size, self.overlap)
        x_positions = self._calculate_positions(width, self.tile_size, self.overlap)
        
        return len(y_positions) * len(x_positions)
    
    def get_tile_info(self, image_shape: Tuple[int, int]) -> dict:
        """
        Get information about tiling for an image.
        
        Args:
            image_shape: Shape of image (height, width)
        
        Returns:
            Dictionary with tiling information
        """
        height, width = image_shape
        tile_count = self.estimate_tile_count(image_shape)
        
        return {
            "image_shape": image_shape,
            "tile_size": self.tile_size,
            "overlap": self.overlap,
            "effective_tile_size": self.effective_tile_size,
            "estimated_tile_count": tile_count,
            "memory_reduction_factor": (height * width) / (self.tile_size * self.tile_size),
        }


def create_test_image(height: int, width: int, num_cells: int = 10) -> np.ndarray:
    """
    Create a synthetic test image with circular cells for testing.
    
    Args:
        height: Image height
        width: Image width
        num_cells: Number of cells to generate
    
    Returns:
        Synthetic image with cells
    """
    image = np.zeros((height, width), dtype=np.uint8)
    
    # Generate random cell positions and sizes
    np.random.seed(42)  # For reproducible tests
    
    for i in range(num_cells):
        # Random position
        center_y = np.random.randint(20, height - 20)
        center_x = np.random.randint(20, width - 20)
        
        # Random radius
        radius = np.random.randint(8, 20)
        
        # Draw circle
        y, x = np.ogrid[:height, :width]
        mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
        image[mask] = 255
    
    return image


def create_test_mask(height: int, width: int, num_cells: int = 10) -> np.ndarray:
    """
    Create a synthetic segmentation mask for testing.
    
    Args:
        height: Image height
        width: Image width
        num_cells: Number of cells to generate
    
    Returns:
        Segmentation mask with labeled cells
    """
    mask = np.zeros((height, width), dtype=np.int32)
    
    # Generate random cell positions and sizes
    np.random.seed(42)  # For reproducible tests
    
    for i in range(num_cells):
        # Random position
        center_y = np.random.randint(20, height - 20)
        center_x = np.random.randint(20, width - 20)
        
        # Random radius
        radius = np.random.randint(8, 20)
        
        # Draw circle with unique label
        y, x = np.ogrid[:height, :width]
        circle_mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
        
        # Only assign if not overlapping with existing cells
        if np.all(mask[circle_mask] == 0):
            mask[circle_mask] = i + 1
    
    return mask