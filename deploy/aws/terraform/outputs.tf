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

output "vpc_id" {
  value       = aws_vpc.main.id
  description = "VPC ID for ECS resources"
}

output "ecs_cluster_arn" {
  value       = aws_ecs_cluster.main.arn
  description = "ECS cluster ARN"
}

output "api_service_name" {
  value       = aws_ecs_service.api.name
  description = "ECS service name for FastAPI"
}

output "dashboard_service_name" {
  value       = aws_ecs_service.dashboard.name
  description = "ECS service name for Streamlit dashboard"
}

output "alb_dns_name" {
  value       = aws_lb.main.dns_name
  description = "ALB DNS name for accessing the application"
}

output "alb_arn" {
  value       = aws_lb.main.arn
  description = "ALB ARN"
}

# NOTE: Connection string output removed - database password is now managed
# via AWS Secrets Manager. The connection string will be constructed by the
# application using the password from the 'database' secret group.

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

