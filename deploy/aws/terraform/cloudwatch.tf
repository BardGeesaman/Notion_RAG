# CloudWatch log groups for ECS tasks.

resource "aws_cloudwatch_log_group" "api" {
  name              = "/ecs/amprenta-api"
  retention_in_days = 30

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
  }
}

resource "aws_cloudwatch_log_group" "dashboard" {
  name              = "/ecs/amprenta-dashboard"
  retention_in_days = 30

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
  }
}


