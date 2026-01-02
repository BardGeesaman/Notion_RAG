# Secrets Manager secrets for ECS services.
# Based on inventory from docs/SECRETS_INVENTORY.md

# Group 1: Database credentials and connection strings
resource "aws_secretsmanager_secret" "database" {
  name = "${var.project_name}/${var.environment}/database"

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
    SecretGroup = "database"
  }
}

resource "aws_secretsmanager_secret_version" "database" {
  secret_id = aws_secretsmanager_secret.database.id

  # NOTE: Secret values must be populated manually in AWS Console after terraform apply
  # This creates the secret structure but does not store actual passwords in Terraform state
  secret_string = jsonencode({
    POSTGRES_PASSWORD = "" # Required: Database password
    POSTGRES_URL      = "" # Optional: Full connection string override
  })
}

# Group 2: External API keys for third-party services
resource "aws_secretsmanager_secret" "api_keys" {
  name = "${var.project_name}/${var.environment}/api-keys"

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
    SecretGroup = "api-keys"
  }
}

resource "aws_secretsmanager_secret_version" "api_keys" {
  secret_id = aws_secretsmanager_secret.api_keys.id

  # NOTE: Secret values must be populated manually in AWS Console after terraform apply
  secret_string = jsonencode({
    OPENAI_API_KEY = "" # Required: OpenAI API key for LLM operations
    ZOTERO_API_KEY = "" # Required: Zotero API key for literature management
    NOTION_API_KEY = "" # Optional: Notion API key (only if ENABLE_NOTION_SYNC=true)
    GEO_API_KEY    = "" # Optional: NCBI GEO API key for higher rate limits
  })
}

# Group 3: Authentication credentials
resource "aws_secretsmanager_secret" "auth" {
  name = "${var.project_name}/${var.environment}/auth"

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
    SecretGroup = "auth"
  }
}

resource "aws_secretsmanager_secret_version" "auth" {
  secret_id = aws_secretsmanager_secret.auth.id

  # NOTE: Secret values must be populated manually in AWS Console after terraform apply
  secret_string = jsonencode({
    NCBI_EMAIL           = "" # Required: Email required by NCBI for API access
    SIGNATURE_SECRET_KEY = "" # HMAC key for electronic signatures
    JWT_SECRET_KEY       = "" # JWT signing key for JupyterHub
  })
}

# Group 4: Integration identifiers (Notion workspace IDs)
resource "aws_secretsmanager_secret" "integrations" {
  name = "${var.project_name}/${var.environment}/integrations"

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
    SecretGroup = "integrations"
  }
}

resource "aws_secretsmanager_secret_version" "integrations" {
  secret_id = aws_secretsmanager_secret.integrations.id

  # NOTE: Secret values must be populated manually in AWS Console after terraform apply
  secret_string = jsonencode({
    NOTION_EXP_DATA_DB_ID            = "" # Notion Experimental Data Assets DB ID
    NOTION_METABOLITE_FEATURES_DB_ID = "" # Notion Metabolite Features DB ID
    NOTION_PROTEIN_FEATURES_DB_ID    = "" # Notion Protein Features DB ID
    NOTION_GENE_FEATURES_DB_ID       = "" # Notion Gene Features DB ID
    NOTION_SIGNATURE_DB_ID           = "" # Notion Lipid Signatures DB ID
    NOTION_SIGNATURE_COMPONENT_DB_ID = "" # Notion Signature Components DB ID
    NOTION_LIPID_SPECIES_DB_ID       = "" # Notion Lipid Species DB ID
    NOTION_PROGRAMS_DB_ID            = "" # Notion Programs DB ID
    NOTION_EXPERIMENTS_DB_ID         = "" # Notion Experiments DB ID
    NOTION_COMPOUND_FEATURES_DB_ID   = "" # Notion Compound Features DB ID
    NOTION_HTS_CAMPAIGNS_DB_ID       = "" # Notion HTS Campaigns DB ID
    NOTION_BIOCHEMICAL_HITS_DB_ID    = "" # Notion Biochemical Hits DB ID
    NOTION_PATHWAYS_DB_ID            = "" # Notion Pathways DB ID
  })
}

# Group 5: Backup and encryption keys
resource "aws_secretsmanager_secret" "backup" {
  name = "${var.project_name}/${var.environment}/backup"

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
    SecretGroup = "backup"
  }
}

resource "aws_secretsmanager_secret_version" "backup" {
  secret_id = aws_secretsmanager_secret.backup.id

  # NOTE: Secret values must be populated manually in AWS Console after terraform apply
  secret_string = jsonencode({
    BACKUP_KMS_KEY_ID = "" # Optional: AWS KMS key ID for backup encryption
  })
}

# Group 6: Email service credentials
resource "aws_secretsmanager_secret" "email" {
  name = "${var.project_name}/${var.environment}/email"

  tags = {
    Project     = var.project_name
    Environment = var.environment
    ManagedBy   = "terraform"
    SecretGroup = "email"
  }
}

resource "aws_secretsmanager_secret_version" "email" {
  secret_id = aws_secretsmanager_secret.email.id

  # NOTE: Secret values must be populated manually in AWS Console after terraform apply
  secret_string = jsonencode({
    SMTP_HOST     = "smtp.gmail.com" # Default SMTP host
    SMTP_PORT     = "587"            # Default SMTP port
    SMTP_USER     = ""               # Required: SMTP username
    SMTP_PASSWORD = ""               # Required: SMTP password/app password
    FROM_EMAIL    = ""               # Optional: Override from address
  })
}


