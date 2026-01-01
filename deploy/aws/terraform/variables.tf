variable "aws_region" {
  description = "AWS region"
  type        = string
  default     = "us-east-1"
}

variable "project_name" {
  description = "Project name prefix"
  type        = string
  default     = "amprenta"
}

variable "app_name" {
  description = "Application name for resource naming"
  type        = string
  default     = "amprenta-rag"
}

variable "vpc_cidr" {
  description = "CIDR block for the ECS VPC"
  type        = string
  default     = "10.0.0.0/16"
}

variable "desired_count" {
  description = "Desired task count per service"
  type        = number
  default     = 1
}

variable "api_cpu" {
  description = "CPU units for API task"
  type        = number
  default     = 256
}

variable "api_memory" {
  description = "Memory (MiB) for API task"
  type        = number
  default     = 512
}

variable "dashboard_cpu" {
  description = "CPU units for dashboard task"
  type        = number
  default     = 256
}

variable "dashboard_memory" {
  description = "Memory (MiB) for dashboard task"
  type        = number
  default     = 512
}

variable "db_name" {
  description = "RDS database name"
  type        = string
}

variable "db_username" {
  description = "RDS master username"
  type        = string
  sensitive   = true
}

# NOTE: db_password variable removed - database passwords are now managed 
# via AWS Secrets Manager for enhanced security. Populate the password
# directly in the 'database' secret after running terraform apply.

variable "db_instance_class" {
  description = "RDS instance class"
  type        = string
  default     = "db.t3.micro"
}

variable "db_allowed_cidrs" {
  description = "CIDR blocks allowed to reach the RDS instance"
  type        = list(string)
}

variable "environment" {
  description = "Environment tag (e.g., dev, staging, prod)"
  type        = string
  default     = "dev"
}

variable "lightsail_blueprint_id" {
  description = "Lightsail blueprint ID"
  type        = string
  default     = "ubuntu_22_04"
}

variable "lightsail_bundle_id" {
  description = "Lightsail bundle ID"
  type        = string
  default     = "nano_2_0"
}

variable "allowed_ssh_cidrs" {
  description = "CIDR blocks allowed for SSH to Lightsail"
  type        = list(string)
  default     = ["0.0.0.0/0"]
}

# Backup configuration
variable "backup_retention_period" {
  description = "Days to retain RDS automated backups (0 to disable)"
  type        = number
  default     = 7
}

variable "backup_retention_days" {
  description = "Days to retain S3 backup files before expiration"
  type        = number
  default     = 365
}

