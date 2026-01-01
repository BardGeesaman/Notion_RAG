# IAM roles for ECS tasks (execution + task role).

# ECS Task Execution Role (for pulling images, writing logs)
resource "aws_iam_role" "ecs_task_execution" {
  name = "${var.project_name}-ecs-task-execution"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Action = "sts:AssumeRole"
      Effect = "Allow"
      Principal = {
        Service = "ecs-tasks.amazonaws.com"
      }
    }]
  })

  tags = {
    Project     = var.project_name
    Environment = var.environment
  }
}

# Attach AWS managed policy for ECS task execution
resource "aws_iam_role_policy_attachment" "ecs_task_execution" {
  role       = aws_iam_role.ecs_task_execution.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}

# Policy for Secrets Manager access (execution role)
resource "aws_iam_role_policy" "ecs_secrets_execution" {
  name = "${var.project_name}-ecs-secrets-execution"
  role = aws_iam_role.ecs_task_execution.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect = "Allow"
      Action = [
        "secretsmanager:GetSecretValue",
        "secretsmanager:DescribeSecret"
      ]
      Resource = [
        aws_secretsmanager_secret.database.arn,
        aws_secretsmanager_secret.api_keys.arn,
        aws_secretsmanager_secret.auth.arn,
        aws_secretsmanager_secret.integrations.arn,
        aws_secretsmanager_secret.backup.arn
      ]
    }]
  })
}

# Policy for Secrets Manager access (task role)
resource "aws_iam_role_policy" "ecs_secrets_task" {
  name = "${var.project_name}-ecs-secrets-task"
  role = aws_iam_role.ecs_task.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect = "Allow"
      Action = [
        "secretsmanager:GetSecretValue",
        "secretsmanager:DescribeSecret"
      ]
      Resource = [
        aws_secretsmanager_secret.database.arn,
        aws_secretsmanager_secret.api_keys.arn,
        aws_secretsmanager_secret.auth.arn,
        aws_secretsmanager_secret.integrations.arn,
        aws_secretsmanager_secret.backup.arn
      ]
    }]
  })
}

# ECS Task Role (for app-level AWS access)
resource "aws_iam_role" "ecs_task" {
  name = "${var.project_name}-ecs-task"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Action = "sts:AssumeRole"
      Effect = "Allow"
      Principal = {
        Service = "ecs-tasks.amazonaws.com"
      }
    }]
  })

  tags = {
    Project     = var.project_name
    Environment = var.environment
  }
}


