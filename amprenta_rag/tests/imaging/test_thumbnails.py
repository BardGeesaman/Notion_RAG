"""Tests for thumbnail service."""

import pytest
import numpy as np
from pathlib import Path
from uuid import uuid4
from unittest.mock import MagicMock, patch

from amprenta_rag.imaging.thumbnails import ThumbnailService


class TestThumbnailService:
    """Test ThumbnailService functionality."""

    @pytest.fixture
    def mock_storage(self):
        """Create mock ImageStorage."""
        storage = MagicMock()
        # Return simple 100x100 grayscale image
        storage.load_image.return_value = np.random.randint(
            0, 255, (100, 100), dtype=np.uint8
        )
        return storage

    @pytest.fixture
    def thumbnail_service(self, mock_storage, tmp_path):
        """Create ThumbnailService with temp cache dir."""
        return ThumbnailService(mock_storage, cache_dir=str(tmp_path))

    def test_get_cache_path(self, thumbnail_service):
        """Test cache path generation."""
        image_id = uuid4()
        path = thumbnail_service.get_cache_path(image_id, 256)
        
        assert str(image_id.hex) in str(path)
        assert "_256.jpg" in str(path)

    def test_exists_returns_false_when_not_cached(self, thumbnail_service):
        """Test exists() returns False for uncached images."""
        image_id = uuid4()
        assert thumbnail_service.exists(image_id) is False

    def test_get_or_generate_creates_thumbnail(self, thumbnail_service, mock_storage):
        """Test thumbnail generation."""
        mock_image = MagicMock()
        mock_image.id = uuid4()
        mock_image.image_path = "/test/image.tiff"
        
        result = thumbnail_service.get_or_generate(mock_image, size=256)
        
        assert result is not None
        assert Path(result).exists()
        assert "_256.jpg" in result
        mock_storage.load_image.assert_called_once()

    def test_get_or_generate_returns_cached(self, thumbnail_service, mock_storage):
        """Test cached thumbnail is returned without regeneration."""
        mock_image = MagicMock()
        mock_image.id = uuid4()
        mock_image.image_path = "/test/image.tiff"
        
        # First call generates
        result1 = thumbnail_service.get_or_generate(mock_image)
        # Second call should return cached
        result2 = thumbnail_service.get_or_generate(mock_image)
        
        assert result1 == result2
        # Storage should only be called once
        assert mock_storage.load_image.call_count == 1

    def test_delete_removes_thumbnails(self, thumbnail_service, mock_storage):
        """Test thumbnail deletion."""
        mock_image = MagicMock()
        mock_image.id = uuid4()
        mock_image.image_path = "/test/image.tiff"
        
        # Generate thumbnail
        thumbnail_service.get_or_generate(mock_image)
        assert thumbnail_service.exists(mock_image.id)
        
        # Delete it
        deleted = thumbnail_service.delete(mock_image.id)
        
        assert deleted >= 1
        assert not thumbnail_service.exists(mock_image.id)

    def test_clear_cache(self, thumbnail_service, mock_storage):
        """Test cache clearing."""
        # Generate multiple thumbnails
        for _ in range(3):
            mock_image = MagicMock()
            mock_image.id = uuid4()
            mock_image.image_path = "/test/image.tiff"
            thumbnail_service.get_or_generate(mock_image)
        
        # Clear cache
        deleted = thumbnail_service.clear_cache()
        
        assert deleted == 3

    def test_handles_16bit_images(self, thumbnail_service):
        """Test 16-bit image conversion."""
        mock_image = MagicMock()
        mock_image.id = uuid4()
        mock_image.image_path = "/test/image.tiff"
        
        # Return 16-bit image
        thumbnail_service.storage.load_image.return_value = np.random.randint(
            0, 65535, (100, 100), dtype=np.uint16
        )
        
        result = thumbnail_service.get_or_generate(mock_image)
        
        assert result is not None
        assert Path(result).exists()
