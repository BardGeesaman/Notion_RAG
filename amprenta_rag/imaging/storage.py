"""Image storage backend for microscopy data."""

from __future__ import annotations

import os
import json
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, Union, Tuple
from abc import ABC, abstractmethod

try:
    import boto3
    from botocore.exceptions import ClientError
    HAS_S3 = True
except ImportError:
    HAS_S3 = False

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False


class StorageBackend(ABC):
    """Abstract base class for image storage backends."""

    @abstractmethod
    def save_file(self, data: bytes, path: str) -> str:
        """Save file data to storage and return the storage path."""
        pass

    @abstractmethod
    def load_file(self, path: str) -> bytes:
        """Load file data from storage."""
        pass

    @abstractmethod
    def exists(self, path: str) -> bool:
        """Check if file exists in storage."""
        pass

    @abstractmethod
    def delete_file(self, path: str) -> bool:
        """Delete file from storage. Returns True if successful."""
        pass


class LocalStorageBackend(StorageBackend):
    """Local filesystem storage backend."""

    def __init__(self, base_path: str = "data/imaging"):
        self.base_path = Path(base_path)
        self.base_path.mkdir(parents=True, exist_ok=True)

    def save_file(self, data: bytes, path: str) -> str:
        """Save file data to local filesystem."""
        full_path = self.base_path / path
        full_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(full_path, 'wb') as f:
            f.write(data)
        
        return str(full_path)

    def load_file(self, path: str) -> bytes:
        """Load file data from local filesystem."""
        full_path = self.base_path / path
        
        if not full_path.exists():
            raise FileNotFoundError(f"File not found: {full_path}")
        
        with open(full_path, 'rb') as f:
            return f.read()

    def exists(self, path: str) -> bool:
        """Check if file exists in local filesystem."""
        full_path = self.base_path / path
        return full_path.exists()

    def delete_file(self, path: str) -> bool:
        """Delete file from local filesystem."""
        try:
            full_path = self.base_path / path
            if full_path.exists():
                full_path.unlink()
                return True
            return False
        except Exception:
            return False


class S3StorageBackend(StorageBackend):
    """AWS S3 storage backend."""

    def __init__(self, bucket_name: str, prefix: str = "imaging/", 
                 aws_access_key_id: Optional[str] = None,
                 aws_secret_access_key: Optional[str] = None,
                 region_name: str = "us-east-1"):
        if not HAS_S3:
            raise ImportError("boto3 is required for S3 storage backend")
        
        self.bucket_name = bucket_name
        self.prefix = prefix
        
        # Initialize S3 client
        session = boto3.Session(
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
            region_name=region_name
        )
        self.s3_client = session.client('s3')

    def _get_s3_key(self, path: str) -> str:
        """Get S3 key from relative path."""
        return f"{self.prefix.rstrip('/')}/{path.lstrip('/')}"

    def save_file(self, data: bytes, path: str) -> str:
        """Save file data to S3."""
        s3_key = self._get_s3_key(path)
        
        try:
            self.s3_client.put_object(
                Bucket=self.bucket_name,
                Key=s3_key,
                Body=data
            )
            return f"s3://{self.bucket_name}/{s3_key}"
        except ClientError as e:
            raise RuntimeError(f"Failed to save to S3: {e}")

    def load_file(self, path: str) -> bytes:
        """Load file data from S3."""
        s3_key = self._get_s3_key(path)
        
        try:
            response = self.s3_client.get_object(
                Bucket=self.bucket_name,
                Key=s3_key
            )
            return response['Body'].read()
        except ClientError as e:
            if e.response['Error']['Code'] == 'NoSuchKey':
                raise FileNotFoundError(f"File not found in S3: s3://{self.bucket_name}/{s3_key}")
            raise RuntimeError(f"Failed to load from S3: {e}")

    def exists(self, path: str) -> bool:
        """Check if file exists in S3."""
        s3_key = self._get_s3_key(path)
        
        try:
            self.s3_client.head_object(
                Bucket=self.bucket_name,
                Key=s3_key
            )
            return True
        except ClientError as e:
            if e.response['Error']['Code'] == '404':
                return False
            raise RuntimeError(f"Failed to check S3 existence: {e}")

    def delete_file(self, path: str) -> bool:
        """Delete file from S3."""
        s3_key = self._get_s3_key(path)
        
        try:
            self.s3_client.delete_object(
                Bucket=self.bucket_name,
                Key=s3_key
            )
            return True
        except ClientError:
            return False


class ImageStorage:
    """High-level image storage interface for microscopy data."""

    def __init__(self, backend: StorageBackend):
        self.backend = backend

    @classmethod
    def create_local(cls, base_path: str = "data/imaging") -> "ImageStorage":
        """Create ImageStorage with local filesystem backend."""
        return cls(LocalStorageBackend(base_path))

    @classmethod
    def create_s3(cls, bucket_name: str, prefix: str = "imaging/", **kwargs) -> "ImageStorage":
        """Create ImageStorage with S3 backend."""
        return cls(S3StorageBackend(bucket_name, prefix, **kwargs))

    def save_image(self, image_data: Union[np.ndarray, bytes], 
                   well_id: str, channel: str, z_slice: int = 0, 
                   timepoint: int = 0, format: str = "tiff") -> str:
        """
        Save microscopy image data.
        
        Args:
            image_data: Image data as numpy array or bytes
            well_id: Well identifier
            channel: Channel name (e.g., "DAPI", "GFP")
            z_slice: Z-stack slice index
            timepoint: Time series index
            format: Image format ("tiff", "png", "jpg")
        
        Returns:
            Storage path of saved image
        """
        # Generate storage path
        path = f"images/{well_id}/{channel}_z{z_slice:03d}_t{timepoint:03d}.{format}"
        
        # Convert numpy array to bytes if needed
        if isinstance(image_data, np.ndarray):
            if not HAS_PIL:
                raise ImportError("PIL is required for numpy array conversion")
            
            # Convert to PIL Image and save to bytes
            if image_data.dtype != np.uint8 and image_data.dtype != np.uint16:
                # Normalize to 0-255 for 8-bit or 0-65535 for 16-bit
                if image_data.max() <= 1.0:
                    image_data = (image_data * 255).astype(np.uint8)
                else:
                    image_data = image_data.astype(np.uint16)
            
            pil_image = Image.fromarray(image_data)
            
            import io
            buffer = io.BytesIO()
            pil_image.save(buffer, format=format.upper())
            image_bytes = buffer.getvalue()
        else:
            image_bytes = image_data

        # Save to backend
        storage_path = self.backend.save_file(image_bytes, path)
        return path  # Return relative path for database storage

    def save_mask(self, mask_data: np.ndarray, segmentation_id: str) -> str:
        """
        Save segmentation mask as NPY file.
        
        Args:
            mask_data: Segmentation mask as numpy array
            segmentation_id: Segmentation identifier
        
        Returns:
            Storage path of saved mask
        """
        path = f"masks/{segmentation_id}_mask.npy"
        
        # Convert numpy array to bytes
        import io
        buffer = io.BytesIO()
        np.save(buffer, mask_data)
        mask_bytes = buffer.getvalue()
        
        # Save to backend
        storage_path = self.backend.save_file(mask_bytes, path)
        return path  # Return relative path for database storage

    def load_image(self, path: str) -> np.ndarray:
        """
        Load microscopy image from storage.
        
        Args:
            path: Storage path of image
        
        Returns:
            Image data as numpy array
        """
        if not HAS_PIL:
            raise ImportError("PIL is required for image loading")
        
        # Load bytes from backend
        image_bytes = self.backend.load_file(path)
        
        # Convert bytes to PIL Image and then to numpy array
        import io
        buffer = io.BytesIO(image_bytes)
        pil_image = Image.open(buffer)
        image_array = np.array(pil_image)
        
        return image_array

    def load_mask(self, path: str) -> np.ndarray:
        """
        Load segmentation mask from storage.
        
        Args:
            path: Storage path of mask
        
        Returns:
            Mask data as numpy array
        """
        # Load bytes from backend
        mask_bytes = self.backend.load_file(path)
        
        # Convert bytes to numpy array
        import io
        buffer = io.BytesIO(mask_bytes)
        mask_array = np.load(buffer)
        
        return mask_array

    def exists(self, path: str) -> bool:
        """Check if file exists in storage."""
        return self.backend.exists(path)

    def delete_image(self, path: str) -> bool:
        """Delete image from storage."""
        return self.backend.delete_file(path)

    def delete_mask(self, path: str) -> bool:
        """Delete mask from storage."""
        return self.backend.delete_file(path)

    def get_image_metadata(self, path: str) -> Dict[str, Any]:
        """
        Get metadata for stored image.
        
        Args:
            path: Storage path of image
        
        Returns:
            Dictionary with image metadata
        """
        if not self.exists(path):
            raise FileNotFoundError(f"Image not found: {path}")
        
        # Load image to get dimensions and properties
        try:
            image_array = self.load_image(path)
            
            metadata = {
                "width": image_array.shape[1] if len(image_array.shape) >= 2 else image_array.shape[0],
                "height": image_array.shape[0],
                "channels": image_array.shape[2] if len(image_array.shape) == 3 else 1,
                "dtype": str(image_array.dtype),
                "size_bytes": image_array.nbytes,
            }
            
            return metadata
        except Exception as e:
            return {"error": str(e)}
