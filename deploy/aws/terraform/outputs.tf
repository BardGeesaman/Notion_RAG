output "rds_endpoint" {
  value       = aws_db_instance.postgres.address
  description = "RDS endpoint"
}

output "rds_port" {
  value       = aws_db_instance.postgres.port
  description = "RDS port"
}

output "lightsail_ip" {
  value       = aws_lightsail_static_ip.app_ip.ip_address
  description = "Lightsail public IP"
}

output "connection_string" {
  value       = "postgresql://${var.db_username}:${var.db_password}@${aws_db_instance.postgres.address}:${aws_db_instance.postgres.port}/${var.db_name}"
  description = "Connection string for convenience"
  sensitive   = true
}

output "environment" {
  value       = var.environment
  description = "Deployment environment"
}

output "ecr_api_repository_url" {
  value       = aws_ecr_repository.api.repository_url
  description = "ECR repository URL for the FastAPI backend image"
}

output "ecr_dashboard_repository_url" {
  value       = aws_ecr_repository.dashboard.repository_url
  description = "ECR repository URL for the Streamlit dashboard image"
}

