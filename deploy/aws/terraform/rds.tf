resource "aws_security_group" "rds_sg" {
  name        = "${var.project_name}-rds-sg"
  description = "RDS access"
  vpc_id      = aws_vpc.main.id

  ingress {
    from_port   = 5432
    to_port     = 5432
    protocol    = "tcp"
    cidr_blocks = var.db_allowed_cidrs
  }

  ingress {
    from_port       = 5432
    to_port         = 5432
    protocol        = "tcp"
    security_groups = [aws_security_group.ecs_sg.id]
    description     = "Postgres from ECS tasks"
  }

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_db_instance" "postgres" {
  identifier             = "${var.project_name}-postgres"
  allocated_storage      = 20
  engine                 = "postgres"
  engine_version         = "15"
  instance_class         = var.db_instance_class
  username               = var.db_username
  password               = var.db_password
  db_name                = var.db_name
  publicly_accessible    = false
  skip_final_snapshot    = var.environment != "prod"
  deletion_protection    = var.environment == "prod"
  vpc_security_group_ids = [aws_security_group.rds_sg.id]
  db_subnet_group_name   = aws_db_subnet_group.rds.name

  tags = {
    Name        = "${var.project_name}-postgres"
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
  }
}

