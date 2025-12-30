# S3 bucket for automated database backups
resource "aws_s3_bucket" "backup_bucket" {
  bucket = "${var.app_name}-backups-${var.environment}"

  tags = {
    Name        = "${var.app_name} Backups"
    Environment = var.environment
    Purpose     = "Database backups and disaster recovery"
  }
}

# Enable versioning for backup files
resource "aws_s3_bucket_versioning" "backup_versioning" {
  bucket = aws_s3_bucket.backup_bucket.id
  versioning_configuration {
    status = "Enabled"
  }
}

# Server-side encryption with KMS
resource "aws_s3_bucket_server_side_encryption_configuration" "backup_encryption" {
  bucket = aws_s3_bucket.backup_bucket.id

  rule {
    apply_server_side_encryption_by_default {
      kms_master_key_id = aws_kms_key.backup_key.arn
      sse_algorithm     = "aws:kms"
    }
    bucket_key_enabled = true
  }
}

# KMS key for backup encryption
resource "aws_kms_key" "backup_key" {
  description             = "${var.app_name} backup encryption key"
  deletion_window_in_days = 7

  tags = {
    Name        = "${var.app_name}-backup-key"
    Environment = var.environment
  }
}

resource "aws_kms_alias" "backup_key_alias" {
  name          = "alias/${var.app_name}-backup-${var.environment}"
  target_key_id = aws_kms_key.backup_key.key_id
}

# Lifecycle configuration for cost optimization
resource "aws_s3_bucket_lifecycle_configuration" "backup_lifecycle" {
  bucket = aws_s3_bucket.backup_bucket.id

  rule {
    id     = "backup_lifecycle"
    status = "Enabled"

    # Transition to Glacier after 30 days
    transition {
      days          = 30
      storage_class = "GLACIER"
    }

    # Transition to Deep Archive after 90 days
    transition {
      days          = 90
      storage_class = "DEEP_ARCHIVE"
    }

    # Delete after retention period
    expiration {
      days = var.backup_retention_days
    }

    # Clean up incomplete multipart uploads
    abort_incomplete_multipart_upload {
      days_after_initiation = 7
    }
  }
}

# Block public access
resource "aws_s3_bucket_public_access_block" "backup_pab" {
  bucket = aws_s3_bucket.backup_bucket.id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

# IAM policy for ECS task backup access
resource "aws_iam_policy" "backup_policy" {
  name        = "${var.app_name}-backup-policy-${var.environment}"
  description = "Policy for backup operations on S3 bucket"

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "s3:GetObject",
          "s3:PutObject",
          "s3:DeleteObject",
          "s3:ListBucket"
        ]
        Resource = [
          aws_s3_bucket.backup_bucket.arn,
          "${aws_s3_bucket.backup_bucket.arn}/*"
        ]
      },
      {
        Effect = "Allow"
        Action = [
          "kms:Decrypt",
          "kms:GenerateDataKey"
        ]
        Resource = aws_kms_key.backup_key.arn
      }
    ]
  })
}

# Attach backup policy to ECS task role
resource "aws_iam_role_policy_attachment" "ecs_backup_policy" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = aws_iam_policy.backup_policy.arn
}

# Output bucket name for application configuration
output "backup_bucket_name" {
  description = "Name of the S3 bucket for backups"
  value       = aws_s3_bucket.backup_bucket.bucket
}

output "backup_kms_key_id" {
  description = "KMS key ID for backup encryption"
  value       = aws_kms_key.backup_key.key_id
}