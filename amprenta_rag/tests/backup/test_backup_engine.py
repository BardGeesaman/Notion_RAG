"""Tests for backup engine functionality."""

import gzip
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch, call
from uuid import uuid4

import pytest

from amprenta_rag.backup.backup_engine import BackupEngine, BackupEngineError
from amprenta_rag.database.models import BackupRecord


class TestBackupEngine:
    """Test BackupEngine functionality."""
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    def test_init_success(self, mock_which, mock_get_config, mock_get_backup_client):
        """Test successful backup engine initialization."""
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        
        mock_config = MagicMock()
        mock_get_config.return_value = mock_config
        
        mock_client = MagicMock()
        mock_get_backup_client.return_value = mock_client
        
        engine = BackupEngine()
        
        assert engine.config == mock_config
        assert engine.backup_client == mock_client
        mock_which.assert_any_call("pg_dump")
        mock_which.assert_any_call("pg_restore")
    
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    def test_init_missing_pg_dump(self, mock_which):
        """Test initialization fails when pg_dump is missing."""
        mock_which.side_effect = lambda cmd: None if cmd == "pg_dump" else f"/usr/bin/{cmd}"
        
        with pytest.raises(BackupEngineError, match="pg_dump not found"):
            BackupEngine()
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    @patch('amprenta_rag.backup.backup_engine.db_session')
    def test_create_backup_success(self, mock_db_session, mock_which, mock_get_config, mock_get_backup_client):
        """Test successful backup creation."""
        # Setup mocks
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        
        mock_config = MagicMock()
        mock_config.postgres.host = "localhost"
        mock_config.postgres.port = 5432
        mock_config.postgres.user = "test_user"
        mock_config.postgres.db = "test_db"
        mock_config.postgres.password = "test_pass"
        mock_get_config.return_value = mock_config
        
        mock_client = MagicMock()
        mock_client.upload_file.return_value = {"size": 1024, "checksum": "abc123"}
        mock_get_backup_client.return_value = mock_client
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock backup record
        mock_backup_record = MagicMock(spec=BackupRecord)
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_backup_record
        
        engine = BackupEngine()
        
        # Mock all the internal methods to avoid file operations
        with patch.object(engine, '_run_pg_dump') as mock_pg_dump:
            with patch.object(engine, '_compress_file') as mock_compress:
                with patch.object(engine, '_calculate_checksum') as mock_checksum:
                    with patch.object(engine, '_generate_manifest') as mock_manifest:
                        with patch('amprenta_rag.backup.backup_engine.tempfile.NamedTemporaryFile') as mock_temp_file:
                            with patch('amprenta_rag.backup.backup_engine.Path') as mock_path:
                                # Mock returns
                                mock_checksum.return_value = "test_checksum"
                                mock_manifest.return_value = {"tables": {}, "pg_version": "15.0"}
                                
                                # Mock path operations
                                mock_path_instance = MagicMock()
                                mock_path_instance.stat.return_value.st_size = 1024
                                mock_path_instance.exists.return_value = True
                                mock_path.return_value = mock_path_instance
                                
                                # Mock temp file
                                mock_temp_file.return_value.__enter__.return_value.name = "/tmp/test.sql"
                                
                                backup_id = engine.create_backup("full")
                                
                                assert isinstance(backup_id, type(uuid4()))
                                
                                # Verify methods were called
                                mock_pg_dump.assert_called_once()
                                mock_compress.assert_called_once()
                                mock_checksum.assert_called_once()
                                mock_manifest.assert_called_once()
                                mock_client.upload_file.assert_called_once()
                                
                                # Verify backup record was updated
                                assert mock_backup_record.status == "completed"
    
    def test_create_backup_invalid_type(self):
        """Test create_backup with invalid backup type."""
        with patch('amprenta_rag.backup.backup_engine.shutil.which') as mock_which:
            mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
            
            with patch('amprenta_rag.backup.backup_engine.get_backup_client'):
                with patch('amprenta_rag.backup.backup_engine.get_config'):
                    engine = BackupEngine()
                    
                    with pytest.raises(ValueError, match="backup_type must be"):
                        engine.create_backup("invalid")
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    @patch('amprenta_rag.backup.backup_engine.db_session')
    def test_restore_backup_success(self, mock_db_session, mock_which, mock_get_config, mock_get_backup_client):
        """Test successful backup restoration."""
        # Setup mocks
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        
        mock_config = MagicMock()
        mock_config.postgres.host = "localhost"
        mock_config.postgres.port = 5432
        mock_config.postgres.user = "test_user"
        mock_config.postgres.db = "test_db"
        mock_config.postgres.password = "test_pass"
        mock_get_config.return_value = mock_config
        
        mock_client = MagicMock()
        mock_get_backup_client.return_value = mock_client
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock backup record
        backup_id = uuid4()
        mock_backup_record = MagicMock(spec=BackupRecord)
        mock_backup_record.id = backup_id
        mock_backup_record.status = "completed"
        mock_backup_record.file_path = "backups/test.sql.gz"
        mock_backup_record.checksum_sha256 = "test_checksum"
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_backup_record
        
        engine = BackupEngine()
        
        # Mock all internal methods to avoid file operations
        with patch.object(engine, '_calculate_checksum') as mock_checksum:
            with patch.object(engine, '_decompress_file') as mock_decompress:
                with patch.object(engine, '_run_pg_restore') as mock_pg_restore:
                    with patch('amprenta_rag.backup.backup_engine.tempfile.TemporaryDirectory') as mock_temp_dir:
                        with patch('amprenta_rag.backup.backup_engine.Path') as mock_path:
                            # Mock checksum to match
                            mock_checksum.return_value = "test_checksum"
                            
                            # Mock temp directory
                            mock_temp_dir.return_value.__enter__.return_value = "/tmp/test_dir"
                            
                            # Mock path operations
                            mock_path_instance = MagicMock()
                            mock_path.return_value = mock_path_instance
                            
                            engine.restore_backup(backup_id)
                            
                            # Verify download was called
                            mock_client.download_file.assert_called_once()
                            
                            # Verify methods were called
                            mock_checksum.assert_called_once()
                            mock_decompress.assert_called_once()
                            mock_pg_restore.assert_called_once()
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    @patch('amprenta_rag.backup.backup_engine.db_session')
    def test_restore_backup_not_found(self, mock_db_session, mock_which, mock_get_config, mock_get_backup_client):
        """Test restore backup when backup record not found."""
        # Setup mocks
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        mock_get_config.return_value = MagicMock()
        mock_get_backup_client.return_value = MagicMock()
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = None
        
        engine = BackupEngine()
        backup_id = uuid4()
        
        with pytest.raises(BackupEngineError, match="Backup .* not found"):
            engine.restore_backup(backup_id)
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    @patch('amprenta_rag.backup.backup_engine.db_session')
    @patch('amprenta_rag.backup.backup_engine.tempfile.NamedTemporaryFile')
    def test_verify_backup_success(self, mock_temp_file, mock_db_session, mock_which, mock_get_config, mock_get_backup_client):
        """Test successful backup verification."""
        # Setup mocks
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        mock_get_config.return_value = MagicMock()
        
        mock_client = MagicMock()
        mock_get_backup_client.return_value = mock_client
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock backup record
        backup_id = uuid4()
        mock_backup_record = MagicMock(spec=BackupRecord)
        mock_backup_record.status = "completed"
        mock_backup_record.file_path = "backups/test.sql.gz"
        mock_backup_record.checksum_sha256 = "expected_checksum"
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_backup_record
        
        # Mock temporary file
        with tempfile.NamedTemporaryFile() as real_temp:
            mock_temp_file.return_value.__enter__.return_value.name = real_temp.name
            
            engine = BackupEngine()
            
            # Mock checksum calculation to match
            with patch.object(engine, '_calculate_checksum') as mock_checksum:
                mock_checksum.return_value = "expected_checksum"
                
                result = engine.verify_backup(backup_id)
                
                assert result is True
                mock_client.download_file.assert_called_once()
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    @patch('amprenta_rag.backup.backup_engine.db_session')
    def test_generate_manifest(self, mock_db_session, mock_which, mock_get_config, mock_get_backup_client):
        """Test manifest generation."""
        # Setup mocks
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        mock_get_config.return_value = MagicMock()
        mock_get_backup_client.return_value = MagicMock()
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock database queries
        mock_db.execute.side_effect = [
            MagicMock(fetchone=lambda: ("PostgreSQL 15.0",)),  # version
            MagicMock(fetchall=lambda: [("public", "users", 100), ("public", "posts", 50)]),  # tables
            MagicMock(fetchone=lambda: ("10 MB",)),  # database size
        ]
        
        engine = BackupEngine()
        backup_id = uuid4()
        
        manifest = engine.generate_manifest(backup_id)
        
        assert manifest["backup_id"] == str(backup_id)
        assert "timestamp" in manifest
        assert "pg_version" in manifest
        assert "tables" in manifest
        assert "database_size" in manifest
        assert "public.users" in manifest["tables"]
        assert "public.posts" in manifest["tables"]
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    @patch('amprenta_rag.backup.backup_engine.db_session')
    def test_get_backup_status(self, mock_db_session, mock_which, mock_get_config, mock_get_backup_client):
        """Test getting backup status."""
        # Setup mocks
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        mock_get_config.return_value = MagicMock()
        mock_get_backup_client.return_value = MagicMock()
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock backup record
        backup_id = uuid4()
        mock_backup_record = MagicMock(spec=BackupRecord)
        mock_backup_record.id = backup_id
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_backup_record
        
        engine = BackupEngine()
        
        result = engine.get_backup_status(backup_id)
        
        assert result == mock_backup_record
        mock_db.expunge.assert_called_once_with(mock_backup_record)
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    @patch('amprenta_rag.backup.backup_engine.db_session')
    def test_list_backups(self, mock_db_session, mock_which, mock_get_config, mock_get_backup_client):
        """Test listing backups."""
        # Setup mocks
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        mock_get_config.return_value = MagicMock()
        mock_get_backup_client.return_value = MagicMock()
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock backup records
        mock_record1 = MagicMock(spec=BackupRecord)
        mock_record2 = MagicMock(spec=BackupRecord)
        mock_records = [mock_record1, mock_record2]
        
        mock_query = mock_db.query.return_value
        mock_query.filter_by.return_value = mock_query
        mock_query.order_by.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = mock_records
        
        engine = BackupEngine()
        
        result = engine.list_backups(limit=10, backup_type="full")
        
        assert result == mock_records
        mock_db.query.assert_called_once()
        mock_query.filter_by.assert_called_once_with(backup_type="full")
        mock_query.limit.assert_called_once_with(10)
        
        # Verify records were detached
        mock_db.expunge.assert_any_call(mock_record1)
        mock_db.expunge.assert_any_call(mock_record2)
    
    @patch('amprenta_rag.backup.backup_engine.get_backup_client')
    @patch('amprenta_rag.backup.backup_engine.get_config')
    @patch('amprenta_rag.backup.backup_engine.shutil.which')
    @patch('amprenta_rag.backup.backup_engine.subprocess.run')
    def test_pg_dump_subprocess_error(self, mock_subprocess, mock_which, mock_get_config, mock_get_backup_client):
        """Test handling of pg_dump subprocess errors."""
        # Setup mocks
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
        
        mock_config = MagicMock()
        mock_config.postgres.host = "localhost"
        mock_config.postgres.port = 5432
        mock_config.postgres.user = "test_user"
        mock_config.postgres.db = "test_db"
        mock_config.postgres.password = "test_pass"
        mock_get_config.return_value = mock_config
        
        mock_get_backup_client.return_value = MagicMock()
        
        # Mock subprocess failure
        mock_subprocess.side_effect = subprocess.CalledProcessError(1, "pg_dump", stderr="Connection failed")
        
        engine = BackupEngine()
        
        with tempfile.NamedTemporaryFile() as temp_file:
            temp_path = Path(temp_file.name)
            
            with pytest.raises(BackupEngineError, match="pg_dump failed"):
                engine._run_pg_dump(temp_path, "full")
    
    def test_checksum_calculation(self):
        """Test SHA256 checksum calculation."""
        with patch('amprenta_rag.backup.backup_engine.shutil.which') as mock_which:
            mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ("pg_dump", "pg_restore") else None
            
            with patch('amprenta_rag.backup.backup_engine.get_backup_client'):
                with patch('amprenta_rag.backup.backup_engine.get_config'):
                    engine = BackupEngine()
                    
                    with tempfile.NamedTemporaryFile() as temp_file:
                        temp_path = Path(temp_file.name)
                        temp_path.write_text("test content for checksum")
                        
                        checksum = engine._calculate_checksum(temp_path)
                        
                        # Verify it's a valid SHA256 hex string
                        assert len(checksum) == 64
                        assert all(c in "0123456789abcdef" for c in checksum)
