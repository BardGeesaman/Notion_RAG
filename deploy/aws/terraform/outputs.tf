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

