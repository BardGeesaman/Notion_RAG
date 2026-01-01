"""Thumbnail generation and caching for microscopy images."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, TYPE_CHECKING
from uuid import UUID

import numpy as np
from PIL import Image

if TYPE_CHECKING:
    from amprenta_rag.imaging.models import MicroscopyImage
    from amprenta_rag.imaging.storage import ImageStorage

logger = logging.getLogger(__name__)


class ThumbnailService:
    """Service for generating and caching microscopy image thumbnails."""

    SUPPORTED_SIZES = [128, 256, 512]
    DEFAULT_SIZE = 256

    def __init__(
        self,
        storage: "ImageStorage",
        cache_dir: str = "data/thumbnails"
    ):
        """
        Initialize thumbnail service.

        Args:
            storage: ImageStorage instance for loading source images
            cache_dir: Directory for cached thumbnails
        """
        self.storage = storage
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_cache_path(self, image_id: UUID, size: int) -> Path:
        """Get cache file path for a thumbnail."""
        return self.cache_dir / f"{image_id.hex}_{size}.jpg"

    def exists(self, image_id: UUID, size: int = DEFAULT_SIZE) -> bool:
        """Check if thumbnail exists in cache."""
        return self.get_cache_path(image_id, size).exists()

    def get_or_generate(
        self,
        image: "MicroscopyImage",
        size: int = DEFAULT_SIZE,
        contrast_enhance: bool = True
    ) -> Optional[str]:
        """
        Get cached thumbnail or generate new one.

        Args:
            image: MicroscopyImage record
            size: Thumbnail size in pixels (128, 256, or 512)
            contrast_enhance: Apply contrast enhancement

        Returns:
            Path to thumbnail file, or None on failure
        """
        if size not in self.SUPPORTED_SIZES:
            size = self.DEFAULT_SIZE

        cache_path = self.get_cache_path(image.id, size)

        # Return cached if exists
        if cache_path.exists():
            return str(cache_path)

        # Generate new thumbnail
        try:
            thumbnail_path = self._generate_thumbnail(
                image, cache_path, size, contrast_enhance
            )
            return thumbnail_path
        except Exception as e:
            logger.error(f"Failed to generate thumbnail for {image.id}: {e}")
            return None

    def _generate_thumbnail(
        self,
        image: "MicroscopyImage",
        cache_path: Path,
        size: int,
        contrast_enhance: bool
    ) -> str:
        """
        Generate thumbnail from source image.

        Args:
            image: MicroscopyImage record
            cache_path: Path to save thumbnail
            size: Target size in pixels
            contrast_enhance: Apply auto-contrast

        Returns:
            Path to generated thumbnail
        """
        # Load source image
        img_array = self.storage.load_image(image.image_path)

        # Convert to PIL Image
        if img_array.dtype == np.uint16:
            # Scale 16-bit to 8-bit
            img_array = (img_array / 256).astype(np.uint8)
        
        pil_image = Image.fromarray(img_array)

        # Convert to RGB if grayscale
        if pil_image.mode == "L":
            pil_image = pil_image.convert("RGB")
        elif pil_image.mode not in ("RGB", "RGBA"):
            pil_image = pil_image.convert("RGB")

        # Resize maintaining aspect ratio
        pil_image.thumbnail((size, size), Image.Resampling.LANCZOS)

        # Apply contrast enhancement
        if contrast_enhance:
            from PIL import ImageOps
            pil_image = ImageOps.autocontrast(pil_image, cutoff=1)

        # Save as JPEG
        pil_image.save(cache_path, "JPEG", quality=85, optimize=True)
        
        logger.info(f"Generated thumbnail: {cache_path}")
        return str(cache_path)

    def delete(self, image_id: UUID) -> int:
        """
        Delete all cached thumbnails for an image.

        Args:
            image_id: Image UUID

        Returns:
            Number of thumbnails deleted
        """
        deleted = 0
        for size in self.SUPPORTED_SIZES:
            cache_path = self.get_cache_path(image_id, size)
            if cache_path.exists():
                cache_path.unlink()
                deleted += 1
        return deleted

    def clear_cache(self) -> int:
        """
        Clear all cached thumbnails.

        Returns:
            Number of thumbnails deleted
        """
        deleted = 0
        for cache_file in self.cache_dir.glob("*.jpg"):
            cache_file.unlink()
            deleted += 1
        logger.info(f"Cleared {deleted} cached thumbnails")
        return deleted


# Singleton instance (initialized lazily)
_thumbnail_service: Optional[ThumbnailService] = None


def get_thumbnail_service() -> ThumbnailService:
    """Get or create thumbnail service singleton."""
    global _thumbnail_service
    if _thumbnail_service is None:
        from amprenta_rag.imaging.storage import ImageStorage
        storage = ImageStorage.create_local()
        _thumbnail_service = ThumbnailService(storage)
    return _thumbnail_service
