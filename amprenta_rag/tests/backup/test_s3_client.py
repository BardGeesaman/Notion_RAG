"""Tests for S3 backup client functionality."""

import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.backup.s3_client import (
    BackupMetadata,
    LocalBackupClient,
    S3BackupClient,
    get_backup_client,
)


class TestBackupMetadata:
    """Test BackupMetadata class."""
    
    def test_metadata_creation(self):
        """Test backup metadata creation."""
        timestamp = datetime.now(timezone.utc)
        metadata = BackupMetadata(
            filename="test-backup.sql.gz",
            size=1024,
            timestamp=timestamp,
            checksum="abc123"
        )
        
        assert metadata.filename == "test-backup.sql.gz"
        assert metadata.size == 1024
        assert metadata.timestamp == timestamp
        assert metadata.checksum == "abc123"


class TestLocalBackupClient:
    """Test LocalBackupClient functionality."""
    
    def test_init_creates_backup_directory(self):
        """Test that initialization creates backup directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
                mock_config.return_value.backup.local_dir = temp_dir
                
                client = LocalBackupClient()
                assert client.backup_dir == Path(temp_dir)
                assert client.backup_dir.exists()
    
    def test_upload_file(self):
        """Test uploading (copying) file to local backup directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
                mock_config.return_value.backup.local_dir = temp_dir
                
                client = LocalBackupClient()
                
                # Create source file
                source_file = Path(temp_dir) / "source.sql"
                source_file.write_text("test backup data")
                
                # Upload file
                result = client.upload_file(str(source_file), "backups/test-backup.sql")
                
                # Check result
                assert "path" in result
                assert "size" in result
                assert "checksum" in result
                assert "last_modified" in result
                
                # Check file was copied
                backup_file = client.backup_dir / "backups" / "test-backup.sql"
                assert backup_file.exists()
                assert backup_file.read_text() == "test backup data"
    
    def test_upload_file_not_found(self):
        """Test upload fails when source file doesn't exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
                mock_config.return_value.backup.local_dir = temp_dir
                
                client = LocalBackupClient()
                
                with pytest.raises(FileNotFoundError, match="Backup file not found"):
                    client.upload_file("/nonexistent/file.sql", "backup.sql")
    
    def test_download_file(self):
        """Test downloading (copying) file from local backup directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
                mock_config.return_value.backup.local_dir = temp_dir
                
                client = LocalBackupClient()
                
                # Create backup file
                backup_file = client.backup_dir / "test-backup.sql"
                backup_file.parent.mkdir(parents=True, exist_ok=True)
                backup_file.write_text("backup content")
                
                # Download file
                download_path = Path(temp_dir) / "downloads" / "restored.sql"
                client.download_file("test-backup.sql", str(download_path))
                
                # Check file was copied
                assert download_path.exists()
                assert download_path.read_text() == "backup content"
    
    def test_list_backups(self):
        """Test listing available backups."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
                mock_config.return_value.backup.local_dir = temp_dir
                
                client = LocalBackupClient()
                
                # Create some backup files
                backup1 = client.backup_dir / "backup-1.sql"
                backup2 = client.backup_dir / "backup-2.sql"
                backup1.write_text("backup 1")
                backup2.write_text("backup 2")
                
                # List backups
                backups = client.list_backups()
                
                assert len(backups) == 2
                assert all(isinstance(b, BackupMetadata) for b in backups)
                
                # Check filenames (sorted by timestamp, newest first)
                filenames = [b.filename for b in backups]
                assert "backup-1.sql" in filenames
                assert "backup-2.sql" in filenames
    
    def test_delete_backup(self):
        """Test deleting backup file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
                mock_config.return_value.backup.local_dir = temp_dir
                
                client = LocalBackupClient()
                
                # Create backup file
                backup_file = client.backup_dir / "test-backup.sql"
                backup_file.write_text("backup to delete")
                assert backup_file.exists()
                
                # Delete backup
                client.delete_backup("test-backup.sql")
                
                # Check file was deleted
                assert not backup_file.exists()


class TestS3BackupClient:
    """Test S3BackupClient functionality."""
    
    def test_init_requires_s3_enabled(self):
        """Test that S3 client requires S3 to be enabled."""
        with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
            mock_config.return_value.backup.s3_enabled = False
            
            with pytest.raises(ValueError, match="S3 backup is disabled"):
                S3BackupClient()
    
    def test_init_requires_bucket_name(self):
        """Test that S3 client requires bucket name."""
        with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
            mock_config.return_value.backup.s3_enabled = True
            mock_config.return_value.backup.s3_bucket = ""
            
            with pytest.raises(ValueError, match="S3 bucket not configured"):
                S3BackupClient()
    
    @patch('boto3.client')
    def test_init_success(self, mock_boto3_client):
        """Test successful S3 client initialization."""
        with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
            mock_config.return_value.backup.s3_enabled = True
            mock_config.return_value.backup.s3_bucket = "test-bucket"
            mock_config.return_value.backup.kms_key_id = ""
            
            # Mock boto3 client
            mock_s3_client = MagicMock()
            mock_boto3_client.return_value = mock_s3_client
            
            # Mock successful head_bucket call
            mock_s3_client.head_bucket.return_value = {}
            
            client = S3BackupClient()
            
            assert client.bucket_name == "test-bucket"
            mock_s3_client.head_bucket.assert_called_once_with(Bucket="test-bucket")
    
    @patch('boto3.client')
    def test_upload_file_success(self, mock_boto3_client):
        """Test successful file upload to S3."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
                mock_config.return_value.backup.s3_enabled = True
                mock_config.return_value.backup.s3_bucket = "test-bucket"
                mock_config.return_value.backup.kms_key_id = ""
                
                # Mock boto3 client
                mock_s3_client = MagicMock()
                mock_boto3_client.return_value = mock_s3_client
                mock_s3_client.head_bucket.return_value = {}
                
                # Mock upload response
                mock_s3_client.head_object.return_value = {
                    'ContentLength': 1024,
                    'ETag': '"abc123"',
                    'LastModified': datetime.now(timezone.utc)
                }
                
                client = S3BackupClient()
                
                # Create test file
                test_file = Path(temp_dir) / "test.sql"
                test_file.write_text("test data")
                
                # Upload file
                result = client.upload_file(str(test_file), "backups/test.sql")
                
                # Check result
                assert result['bucket'] == "test-bucket"
                assert result['key'] == "backups/test.sql"
                assert result['size'] == 1024
                assert result['etag'] == "abc123"
                
                # Check upload was called
                mock_s3_client.upload_file.assert_called_once()


class TestGetBackupClient:
    """Test get_backup_client function."""
    
    def test_returns_s3_client_when_enabled(self):
        """Test returns S3 client when S3 is enabled and available."""
        with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
            mock_config.return_value.backup.s3_enabled = True
            mock_config.return_value.backup.s3_bucket = "test-bucket"
            
            with patch('amprenta_rag.backup.s3_client.S3BackupClient') as mock_s3_client:
                mock_instance = MagicMock()
                mock_s3_client.return_value = mock_instance
                
                client = get_backup_client()
                
                assert client == mock_instance
                mock_s3_client.assert_called_once()
    
    def test_falls_back_to_local_when_s3_fails(self):
        """Test falls back to local client when S3 fails."""
        with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
            mock_config.return_value.backup.s3_enabled = True
            mock_config.return_value.backup.local_dir = "/tmp/backups"
            
            with patch('amprenta_rag.backup.s3_client.S3BackupClient') as mock_s3_client:
                mock_s3_client.side_effect = ImportError("boto3 not installed")
                
                with patch('amprenta_rag.backup.s3_client.LocalBackupClient') as mock_local_client:
                    mock_instance = MagicMock()
                    mock_local_client.return_value = mock_instance
                    
                    client = get_backup_client()
                    
                    assert client == mock_instance
                    mock_local_client.assert_called_once()
    
    def test_returns_local_client_when_s3_disabled(self):
        """Test returns local client when S3 is disabled."""
        with patch('amprenta_rag.backup.s3_client.get_config') as mock_config:
            mock_config.return_value.backup.s3_enabled = False
            mock_config.return_value.backup.local_dir = "/tmp/backups"
            
            with patch('amprenta_rag.backup.s3_client.LocalBackupClient') as mock_local_client:
                mock_instance = MagicMock()
                mock_local_client.return_value = mock_instance
                
                client = get_backup_client()
                
                assert client == mock_instance
                mock_local_client.assert_called_once()
