"""Database backup engine with pg_dump/restore functionality."""

import gzip
import hashlib
import logging
import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Any
from uuid import UUID, uuid4

from amprenta_rag.backup.s3_client import get_backup_client
from amprenta_rag.config import get_config
from amprenta_rag.database.models import BackupRecord
from amprenta_rag.database.session import db_session

logger = logging.getLogger(__name__)


class BackupEngineError(Exception):
    """Base exception for backup engine errors."""
    pass


class BackupEngine:
    """Database backup engine for PostgreSQL databases."""
    
    def __init__(self):
        self.config = get_config()
        self.postgres_config = self.config.postgres
        self.backup_client = get_backup_client()
        
        # Verify pg_dump and pg_restore are available
        if not shutil.which("pg_dump"):
            raise BackupEngineError("pg_dump not found in PATH")
        if not shutil.which("pg_restore"):
            raise BackupEngineError("pg_restore not found in PATH")
    
    def create_backup(self, backup_type: str = "full") -> UUID:
        """Create a database backup.
        
        Args:
            backup_type: Type of backup ('full' or 'incremental')
            
        Returns:
            UUID of the created backup record
        """
        if backup_type not in ("full", "incremental"):
            raise ValueError("backup_type must be 'full' or 'incremental'")
        
        backup_id = uuid4()
        
        # Create backup record
        with db_session() as db:
            backup_record = BackupRecord(
                id=backup_id,
                backup_type=backup_type,
                status="pending"
            )
            db.add(backup_record)
            db.commit()
        
        try:
            # Update status to running
            with db_session() as db:
                backup_record = db.query(BackupRecord).filter_by(id=backup_id).first()
                if not backup_record:
                    raise BackupEngineError(f"Backup record {backup_id} not found")
                
                backup_record.status = "running"
                backup_record.started_at = datetime.now(timezone.utc)
                db.add(backup_record)
                db.commit()
            
            # Create temporary file for backup
            with tempfile.NamedTemporaryFile(suffix=".sql", delete=False) as temp_file:
                temp_path = Path(temp_file.name)
            
            try:
                # Run pg_dump
                self._run_pg_dump(temp_path, backup_type)
                
                # Compress backup
                compressed_path = Path(str(temp_path) + ".gz")
                self._compress_file(temp_path, compressed_path)
                
                # Calculate checksum
                checksum = self._calculate_checksum(compressed_path)
                
                # Generate backup key for storage
                timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
                backup_key = f"backups/{backup_type}/{timestamp}_{backup_id}.sql.gz"
                
                # Upload to storage
                self.backup_client.upload_file(
                    str(compressed_path),
                    backup_key
                )
                
                # Generate manifest
                manifest = self._generate_manifest(backup_id)
                
                # Update backup record with success
                with db_session() as db:
                    backup_record = db.query(BackupRecord).filter_by(id=backup_id).first()
                    if backup_record:
                        backup_record.status = "completed"
                        backup_record.completed_at = datetime.now(timezone.utc)
                        backup_record.file_path = backup_key
                        backup_record.file_size_bytes = compressed_path.stat().st_size
                        backup_record.checksum_sha256 = checksum
                        backup_record.backup_metadata = manifest
                        db.add(backup_record)
                        db.commit()
                
                logger.info(f"Backup {backup_id} completed successfully")
                return backup_id
                
            finally:
                # Clean up temporary files
                if temp_path.exists():
                    temp_path.unlink()
                if compressed_path.exists():
                    compressed_path.unlink()
                    
        except Exception as e:
            logger.exception(f"Backup {backup_id} failed: {e}")
            
            # Update backup record with failure
            with db_session() as db:
                backup_record = db.query(BackupRecord).filter_by(id=backup_id).first()
                if backup_record:
                    backup_record.status = "failed"
                    backup_record.completed_at = datetime.now(timezone.utc)
                    backup_record.error_message = str(e)
                    db.add(backup_record)
                    db.commit()
            
            raise BackupEngineError(f"Backup failed: {e}") from e
    
    def restore_backup(self, backup_id: UUID, target_db: Optional[str] = None) -> None:
        """Restore a database backup.
        
        Args:
            backup_id: UUID of the backup to restore
            target_db: Target database name (defaults to current database)
        """
        with db_session() as db:
            backup_record = db.query(BackupRecord).filter_by(id=backup_id).first()
            if not backup_record:
                raise BackupEngineError(f"Backup {backup_id} not found")
            
            if backup_record.status != "completed":
                raise BackupEngineError(f"Backup {backup_id} is not completed (status: {backup_record.status})")
            
            if not backup_record.file_path:
                raise BackupEngineError(f"Backup {backup_id} has no file path")
        
        # Create temporary directory for restore
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Download backup file
            compressed_file = temp_path / "backup.sql.gz"
            self.backup_client.download_file(backup_record.file_path, str(compressed_file))
            
            # Verify checksum
            downloaded_checksum = self._calculate_checksum(compressed_file)
            if downloaded_checksum != backup_record.checksum_sha256:
                raise BackupEngineError(
                    f"Checksum mismatch for backup {backup_id}: "
                    f"expected {backup_record.checksum_sha256}, got {downloaded_checksum}"
                )
            
            # Decompress backup
            sql_file = temp_path / "backup.sql"
            self._decompress_file(compressed_file, sql_file)
            
            # Run pg_restore
            self._run_pg_restore(sql_file, target_db)
            
        logger.info(f"Backup {backup_id} restored successfully")
    
    def verify_backup(self, backup_id: UUID) -> bool:
        """Verify backup integrity by checksum validation.
        
        Args:
            backup_id: UUID of the backup to verify
            
        Returns:
            True if backup is valid, False otherwise
        """
        with db_session() as db:
            backup_record = db.query(BackupRecord).filter_by(id=backup_id).first()
            if not backup_record:
                raise BackupEngineError(f"Backup {backup_id} not found")
            
            if backup_record.status != "completed":
                return False
            
            if not backup_record.file_path or not backup_record.checksum_sha256:
                return False
        
        try:
            # Create temporary file for verification
            with tempfile.NamedTemporaryFile(suffix=".sql.gz") as temp_file:
                temp_path = Path(temp_file.name)
                
                # Download backup file
                self.backup_client.download_file(backup_record.file_path, str(temp_path))
                
                # Verify checksum
                actual_checksum = self._calculate_checksum(temp_path)
                expected_checksum = backup_record.checksum_sha256
                
                if actual_checksum == expected_checksum:
                    logger.info(f"Backup {backup_id} verification successful")
                    return True
                else:
                    logger.error(
                        f"Backup {backup_id} verification failed: "
                        f"expected {expected_checksum}, got {actual_checksum}"
                    )
                    return False
                    
        except Exception as e:
            logger.exception(f"Backup {backup_id} verification error: {e}")
            return False
    
    def generate_manifest(self, backup_id: UUID) -> Dict[str, Any]:
        """Generate manifest for a backup.
        
        Args:
            backup_id: UUID of the backup
            
        Returns:
            Manifest dictionary with backup metadata
        """
        return self._generate_manifest(backup_id)
    
    def get_backup_status(self, backup_id: UUID) -> Optional[BackupRecord]:
        """Get backup record by ID.
        
        Args:
            backup_id: UUID of the backup
            
        Returns:
            BackupRecord if found, None otherwise
        """
        with db_session() as db:
            backup_record = db.query(BackupRecord).filter_by(id=backup_id).first()
            if backup_record:
                db.expunge(backup_record)  # Detach from session
            return backup_record
    
    def list_backups(self, limit: int = 50, backup_type: Optional[str] = None) -> List[BackupRecord]:
        """List backup records.
        
        Args:
            limit: Maximum number of records to return
            backup_type: Filter by backup type ('full' or 'incremental')
            
        Returns:
            List of BackupRecord objects
        """
        with db_session() as db:
            query = db.query(BackupRecord)
            
            if backup_type:
                query = query.filter_by(backup_type=backup_type)
            
            backup_records = query.order_by(BackupRecord.created_at.desc()).limit(limit).all()
            
            # Detach from session
            for record in backup_records:
                db.expunge(record)
            
            return backup_records
    
    def _run_pg_dump(self, output_path: Path, backup_type: str) -> None:
        """Run pg_dump to create database backup."""
        cmd = [
            "pg_dump",
            "--host", self.postgres_config.host,
            "--port", str(self.postgres_config.port),
            "--username", self.postgres_config.user,
            "--dbname", self.postgres_config.db,
            "--no-password",  # Use PGPASSWORD env var
            "--verbose",
            "--file", str(output_path),
        ]
        
        # Add backup-type specific options
        if backup_type == "full":
            cmd.extend([
                "--clean",  # Include DROP statements
                "--create",  # Include CREATE DATABASE
                "--if-exists",  # Use IF EXISTS for drops
            ])
        
        # Set environment for authentication
        env = os.environ.copy()
        if self.postgres_config.password:
            env["PGPASSWORD"] = self.postgres_config.password
        
        try:
            result = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stderr:
                logger.info(f"pg_dump stderr: {result.stderr}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"pg_dump failed: {e.stderr}")
            raise BackupEngineError(f"pg_dump failed: {e}")
    
    def _run_pg_restore(self, sql_file: Path, target_db: Optional[str] = None) -> None:
        """Run pg_restore to restore database backup."""
        db_name = target_db or self.postgres_config.db
        
        cmd = [
            "psql",
            "--host", self.postgres_config.host,
            "--port", str(self.postgres_config.port),
            "--username", self.postgres_config.user,
            "--dbname", db_name,
            "--no-password",
            "--file", str(sql_file),
        ]
        
        # Set environment for authentication
        env = os.environ.copy()
        if self.postgres_config.password:
            env["PGPASSWORD"] = self.postgres_config.password
        
        try:
            result = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stderr:
                logger.info(f"psql stderr: {result.stderr}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"psql failed: {e.stderr}")
            raise BackupEngineError(f"Database restore failed: {e}")
    
    def _compress_file(self, input_path: Path, output_path: Path) -> None:
        """Compress file using gzip."""
        with open(input_path, 'rb') as f_in:
            with gzip.open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    
    def _decompress_file(self, input_path: Path, output_path: Path) -> None:
        """Decompress gzip file."""
        with gzip.open(input_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    
    def _calculate_checksum(self, file_path: Path) -> str:
        """Calculate SHA256 checksum of file."""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()
    
    def _generate_manifest(self, backup_id: UUID) -> Dict[str, Any]:
        """Generate backup manifest with database metadata."""
        manifest = {
            "backup_id": str(backup_id),
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "database": self.postgres_config.db,
            "tables": {},
            "pg_version": None,
        }
        
        try:
            with db_session() as db:
                # Get PostgreSQL version
                version_result = db.execute("SELECT version()").fetchone()
                if version_result:
                    manifest["pg_version"] = version_result[0]
                
                # Get table information
                tables_result = db.execute("""
                    SELECT 
                        schemaname,
                        tablename,
                        n_tup_ins + n_tup_upd + n_tup_del as total_rows
                    FROM pg_stat_user_tables
                    ORDER BY schemaname, tablename
                """).fetchall()
                
                for row in tables_result:
                    schema, table, row_count = row
                    table_key = f"{schema}.{table}"
                    manifest["tables"][table_key] = {
                        "row_count": row_count or 0
                    }
                
                # Get database size
                size_result = db.execute(
                    "SELECT pg_size_pretty(pg_database_size(current_database()))"
                ).fetchone()
                if size_result:
                    manifest["database_size"] = size_result[0]
                    
        except Exception as e:
            logger.warning(f"Failed to generate complete manifest: {e}")
            manifest["error"] = str(e)
        
        return manifest
