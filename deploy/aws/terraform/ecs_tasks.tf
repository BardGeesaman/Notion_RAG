# ECS task definitions (FastAPI + Streamlit).

resource "aws_ecs_task_definition" "api" {
  family                   = "${var.project_name}-${var.environment}-api"
  requires_compatibilities = ["FARGATE"]
  network_mode             = "awsvpc"
  cpu                      = tostring(var.api_cpu)
  memory                   = tostring(var.api_memory)
  execution_role_arn       = aws_iam_role.ecs_task_execution.arn
  task_role_arn            = aws_iam_role.ecs_task.arn

  container_definitions = jsonencode([
    {
      name      = "api"
      image     = "${aws_ecr_repository.api.repository_url}:latest"
      essential = true
      environment = [
        { name = "ENVIRONMENT", value = var.environment }
      ]
      secrets = [
        # Database credentials
        { name = "POSTGRES_URL", valueFrom = "${aws_secretsmanager_secret.database.arn}:POSTGRES_URL::" },
        { name = "POSTGRES_PASSWORD", valueFrom = "${aws_secretsmanager_secret.database.arn}:POSTGRES_PASSWORD::" },

        # API keys for external services
        { name = "OPENAI_API_KEY", valueFrom = "${aws_secretsmanager_secret.api_keys.arn}:OPENAI_API_KEY::" },
        { name = "ZOTERO_API_KEY", valueFrom = "${aws_secretsmanager_secret.api_keys.arn}:ZOTERO_API_KEY::" },
        { name = "NOTION_API_KEY", valueFrom = "${aws_secretsmanager_secret.api_keys.arn}:NOTION_API_KEY::" },
        { name = "GEO_API_KEY", valueFrom = "${aws_secretsmanager_secret.api_keys.arn}:GEO_API_KEY::" },

        # Authentication credentials
        { name = "NCBI_EMAIL", valueFrom = "${aws_secretsmanager_secret.auth.arn}:NCBI_EMAIL::" },

        # Integration identifiers (Notion database IDs)
        { name = "NOTION_EXP_DATA_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_EXP_DATA_DB_ID::" },
        { name = "NOTION_METABOLITE_FEATURES_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_METABOLITE_FEATURES_DB_ID::" },
        { name = "NOTION_PROTEIN_FEATURES_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_PROTEIN_FEATURES_DB_ID::" },
        { name = "NOTION_GENE_FEATURES_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_GENE_FEATURES_DB_ID::" },
        { name = "NOTION_SIGNATURE_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_SIGNATURE_DB_ID::" },
        { name = "NOTION_SIGNATURE_COMPONENT_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_SIGNATURE_COMPONENT_DB_ID::" },
        { name = "NOTION_LIPID_SPECIES_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_LIPID_SPECIES_DB_ID::" },
        { name = "NOTION_PROGRAMS_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_PROGRAMS_DB_ID::" },
        { name = "NOTION_EXPERIMENTS_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_EXPERIMENTS_DB_ID::" },
        { name = "NOTION_COMPOUND_FEATURES_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_COMPOUND_FEATURES_DB_ID::" },
        { name = "NOTION_HTS_CAMPAIGNS_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_HTS_CAMPAIGNS_DB_ID::" },
        { name = "NOTION_BIOCHEMICAL_HITS_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_BIOCHEMICAL_HITS_DB_ID::" },
        { name = "NOTION_PATHWAYS_DB_ID", valueFrom = "${aws_secretsmanager_secret.integrations.arn}:NOTION_PATHWAYS_DB_ID::" },

        # Backup configuration
        { name = "BACKUP_KMS_KEY_ID", valueFrom = "${aws_secretsmanager_secret.backup.arn}:BACKUP_KMS_KEY_ID::" }
      ]
      portMappings = [
        {
          containerPort = 8000
          hostPort      = 8000
          protocol      = "tcp"
        }
      ]
      logConfiguration = {
        logDriver = "awslogs"
        options = {
          awslogs-group         = aws_cloudwatch_log_group.api.name
          awslogs-region        = var.aws_region
          awslogs-stream-prefix = "ecs"
        }
      }
      healthCheck = {
        command     = ["CMD-SHELL", "curl -fsS http://localhost:8000/health || exit 1"]
        interval    = 30
        timeout     = 5
        retries     = 3
        startPeriod = 20
      }
    }
  ])

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
  }
}

resource "aws_ecs_task_definition" "dashboard" {
  family                   = "${var.project_name}-${var.environment}-dashboard"
  requires_compatibilities = ["FARGATE"]
  network_mode             = "awsvpc"
  cpu                      = tostring(var.dashboard_cpu)
  memory                   = tostring(var.dashboard_memory)
  execution_role_arn       = aws_iam_role.ecs_task_execution.arn
  task_role_arn            = aws_iam_role.ecs_task.arn

  container_definitions = jsonencode([
    {
      name      = "dashboard"
      image     = "${aws_ecr_repository.dashboard.repository_url}:latest"
      essential = true
      environment = [
        { name = "ENVIRONMENT", value = var.environment },
        { name = "API_URL", value = "http://${aws_lb.main.dns_name}" }
      ]
      secrets = [
        # Dashboard may need database access for direct queries
        { name = "POSTGRES_URL", valueFrom = "${aws_secretsmanager_secret.database.arn}:POSTGRES_URL::" },
        { name = "POSTGRES_PASSWORD", valueFrom = "${aws_secretsmanager_secret.database.arn}:POSTGRES_PASSWORD::" },

        # API keys that dashboard might need
        { name = "OPENAI_API_KEY", valueFrom = "${aws_secretsmanager_secret.api_keys.arn}:OPENAI_API_KEY::" },
        { name = "NOTION_API_KEY", valueFrom = "${aws_secretsmanager_secret.api_keys.arn}:NOTION_API_KEY::" }
      ]
      portMappings = [
        {
          containerPort = 8501
          hostPort      = 8501
          protocol      = "tcp"
        }
      ]
      logConfiguration = {
        logDriver = "awslogs"
        options = {
          awslogs-group         = aws_cloudwatch_log_group.dashboard.name
          awslogs-region        = var.aws_region
          awslogs-stream-prefix = "ecs"
        }
      }
      healthCheck = {
        command     = ["CMD-SHELL", "curl -fsS http://localhost:8501/ || exit 1"]
        interval    = 30
        timeout     = 5
        retries     = 3
        startPeriod = 20
      }
    }
  ])

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
  }
}


