"""S3 backup client with local filesystem fallback."""

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Dict, Any

from amprenta_rag.config import get_config

logger = logging.getLogger(__name__)


class BackupMetadata:
    """Backup file metadata."""
    
    def __init__(self, filename: str, size: int, timestamp: datetime, checksum: str):
        self.filename = filename
        self.size = size
        self.timestamp = timestamp
        self.checksum = checksum


class S3BackupClient:
    """S3-based backup client with encryption support."""
    
    def __init__(self):
        self.config = get_config().backup
        
        if not self.config.s3_enabled:
            raise ValueError("S3 backup is disabled. Set BACKUP_S3_ENABLED=true")
        
        if not self.config.s3_bucket:
            raise ValueError("S3 bucket not configured. Set BACKUP_S3_BUCKET")
        
        try:
            import boto3
            from botocore.exceptions import ClientError
            self.boto3 = boto3
            self.ClientError = ClientError
        except ImportError as e:
            raise ImportError("boto3 is required for S3 backup. Install with: pip install boto3") from e
        
        self.s3_client = self.boto3.client('s3')
        self.bucket_name = self.config.s3_bucket
        
        # Test S3 connection
        try:
            self.s3_client.head_bucket(Bucket=self.bucket_name)
            logger.info(f"S3 backup client initialized for bucket: {self.bucket_name}")
        except self.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == '404':
                raise ValueError(f"S3 bucket '{self.bucket_name}' does not exist")
            elif error_code == '403':
                raise ValueError(f"Access denied to S3 bucket '{self.bucket_name}'")
            else:
                raise ValueError(f"S3 connection failed: {e}")
    
    def upload_file(self, local_path: str, backup_key: str) -> Dict[str, Any]:
        """Upload backup file to S3."""
        local_file = Path(local_path)
        if not local_file.exists():
            raise FileNotFoundError(f"Backup file not found: {local_path}")
        
        extra_args = {}
        
        # Add KMS encryption if configured
        if self.config.kms_key_id:
            extra_args['ServerSideEncryption'] = 'aws:kms'
            extra_args['SSEKMSKeyId'] = self.config.kms_key_id
        
        # Add metadata
        extra_args['Metadata'] = {
            'backup-timestamp': datetime.now(timezone.utc).isoformat(),
            'backup-type': 'database',
            'original-filename': local_file.name
        }
        
        try:
            self.s3_client.upload_file(
                Filename=str(local_file),
                Bucket=self.bucket_name,
                Key=backup_key,
                ExtraArgs=extra_args
            )
            
            # Get file info
            response = self.s3_client.head_object(Bucket=self.bucket_name, Key=backup_key)
            
            logger.info(f"Backup uploaded to S3: s3://{self.bucket_name}/{backup_key}")
            
            return {
                'bucket': self.bucket_name,
                'key': backup_key,
                'size': response['ContentLength'],
                'etag': response['ETag'].strip('"'),
                'last_modified': response['LastModified'],
                'encrypted': 'ServerSideEncryption' in response
            }
            
        except self.ClientError as e:
            logger.error(f"Failed to upload backup to S3: {e}")
            raise
    
    def download_file(self, backup_key: str, local_path: str) -> None:
        """Download backup file from S3."""
        local_file = Path(local_path)
        local_file.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            self.s3_client.download_file(
                Bucket=self.bucket_name,
                Key=backup_key,
                Filename=str(local_file)
            )
            
            logger.info(f"Backup downloaded from S3: {backup_key} -> {local_path}")
            
        except self.ClientError as e:
            if e.response['Error']['Code'] == 'NoSuchKey':
                raise FileNotFoundError(f"Backup not found in S3: {backup_key}")
            logger.error(f"Failed to download backup from S3: {e}")
            raise
    
    def list_backups(self, prefix: str = "") -> List[BackupMetadata]:
        """List available backups in S3."""
        try:
            response = self.s3_client.list_objects_v2(
                Bucket=self.bucket_name,
                Prefix=prefix
            )
            
            backups = []
            for obj in response.get('Contents', []):
                # Get metadata
                head_response = self.s3_client.head_object(
                    Bucket=self.bucket_name,
                    Key=obj['Key']
                )
                
                backups.append(BackupMetadata(
                    filename=obj['Key'],
                    size=obj['Size'],
                    timestamp=obj['LastModified'],
                    checksum=head_response['ETag'].strip('"')
                ))
            
            # Sort by timestamp, newest first
            backups.sort(key=lambda x: x.timestamp, reverse=True)
            return backups
            
        except self.ClientError as e:
            logger.error(f"Failed to list backups from S3: {e}")
            raise
    
    def delete_backup(self, backup_key: str) -> None:
        """Delete backup file from S3."""
        try:
            self.s3_client.delete_object(
                Bucket=self.bucket_name,
                Key=backup_key
            )
            
            logger.info(f"Backup deleted from S3: {backup_key}")
            
        except self.ClientError as e:
            logger.error(f"Failed to delete backup from S3: {e}")
            raise


class LocalBackupClient:
    """Local filesystem backup client as fallback."""
    
    def __init__(self):
        self.config = get_config().backup
        self.backup_dir = Path(self.config.local_dir)
        self.backup_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Local backup client initialized: {self.backup_dir}")
    
    def upload_file(self, local_path: str, backup_key: str) -> Dict[str, Any]:
        """Copy backup file to local backup directory."""
        source_file = Path(local_path)
        if not source_file.exists():
            raise FileNotFoundError(f"Backup file not found: {local_path}")
        
        # Create subdirectories if needed
        backup_file = self.backup_dir / backup_key
        backup_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Copy file
        import shutil
        shutil.copy2(source_file, backup_file)
        
        # Calculate checksum
        import hashlib
        sha256_hash = hashlib.sha256()
        with open(backup_file, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                sha256_hash.update(chunk)
        
        logger.info(f"Backup copied locally: {backup_key}")
        
        return {
            'path': str(backup_file),
            'size': backup_file.stat().st_size,
            'checksum': sha256_hash.hexdigest(),
            'last_modified': datetime.fromtimestamp(backup_file.stat().st_mtime, tz=timezone.utc)
        }
    
    def download_file(self, backup_key: str, local_path: str) -> None:
        """Copy backup file from local backup directory."""
        backup_file = self.backup_dir / backup_key
        if not backup_file.exists():
            raise FileNotFoundError(f"Backup not found: {backup_key}")
        
        local_file = Path(local_path)
        local_file.parent.mkdir(parents=True, exist_ok=True)
        
        import shutil
        shutil.copy2(backup_file, local_file)
        
        logger.info(f"Backup copied from local storage: {backup_key} -> {local_path}")
    
    def list_backups(self, prefix: str = "") -> List[BackupMetadata]:
        """List available backups in local directory."""
        backups = []
        
        # Find all files matching prefix
        if prefix:
            pattern = f"{prefix}*"
        else:
            pattern = "*"
        
        for backup_file in self.backup_dir.rglob(pattern):
            if backup_file.is_file():
                # Calculate relative path from backup_dir
                relative_path = backup_file.relative_to(self.backup_dir)
                
                # Calculate checksum
                import hashlib
                sha256_hash = hashlib.sha256()
                with open(backup_file, "rb") as f:
                    for chunk in iter(lambda: f.read(4096), b""):
                        sha256_hash.update(chunk)
                
                backups.append(BackupMetadata(
                    filename=str(relative_path),
                    size=backup_file.stat().st_size,
                    timestamp=datetime.fromtimestamp(backup_file.stat().st_mtime, tz=timezone.utc),
                    checksum=sha256_hash.hexdigest()
                ))
        
        # Sort by timestamp, newest first
        backups.sort(key=lambda x: x.timestamp, reverse=True)
        return backups
    
    def delete_backup(self, backup_key: str) -> None:
        """Delete backup file from local directory."""
        backup_file = self.backup_dir / backup_key
        if not backup_file.exists():
            raise FileNotFoundError(f"Backup not found: {backup_key}")
        
        backup_file.unlink()
        
        # Remove empty parent directories
        try:
            backup_file.parent.rmdir()
        except OSError:
            pass  # Directory not empty, which is fine
        
        logger.info(f"Backup deleted from local storage: {backup_key}")


def get_backup_client():
    """Get appropriate backup client based on configuration."""
    config = get_config().backup
    
    if config.s3_enabled:
        try:
            return S3BackupClient()
        except (ImportError, ValueError) as e:
            logger.warning(f"S3 backup client failed, falling back to local: {e}")
            return LocalBackupClient()
    else:
        return LocalBackupClient()