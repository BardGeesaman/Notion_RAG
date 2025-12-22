# Secrets Manager secrets for ECS services.

# Database URL secret (Postgres connection string)
resource "aws_secretsmanager_secret" "db_url" {
  name = "${var.project_name}/${var.environment}/database-url"

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
  }
}

resource "aws_secretsmanager_secret_version" "db_url" {
  secret_id = aws_secretsmanager_secret.db_url.id

  # NOTE: this stores the DB password in Terraform state (since the input var is in state).
  # This is acceptable for now to unblock ECS wiring; we can move to external secret
  # injection later (e.g., SSM/Secrets managed outside Terraform).
  secret_string = "postgresql://${var.db_username}:${var.db_password}@${aws_db_instance.postgres.address}:${aws_db_instance.postgres.port}/${var.db_name}"
}

# API keys secret (JSON), can be empty initially.
resource "aws_secretsmanager_secret" "api_keys" {
  name = "${var.project_name}/${var.environment}/api-keys"

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
  }
}

resource "aws_secretsmanager_secret_version" "api_keys" {
  secret_id = aws_secretsmanager_secret.api_keys.id
  secret_string = jsonencode({
    OPENAI_API_KEY   = ""
    PINECONE_API_KEY = ""
  })
}


