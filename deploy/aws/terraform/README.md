# Terraform AWS Infrastructure

Terraform configuration for deploying Amprenta RAG infrastructure to AWS.

**Components**:
- RDS PostgreSQL database
- Lightsail compute instance
- VPC networking and security groups

---

## Prerequisites

- **Terraform** >= 1.0 ([Install](https://www.terraform.io/downloads))
- **AWS CLI** configured with valid credentials ([Setup](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html))
- AWS account with appropriate permissions (RDS, Lightsail, VPC, EC2)

Verify your setup:

```bash
terraform version
aws sts get-caller-identity
```

---

## Quick Start

### 1. Initialize Terraform

```bash
cd deploy/aws/terraform
terraform init
```

### 2. Configure Variables

Copy the example file and edit with your values:

```bash
cp terraform.tfvars.example terraform.tfvars
```

**Required variables** (must be set):
- `db_name` - Database name
- `db_username` - Database master username
- `db_password` - Database master password (use strong password)
- `db_allowed_cidrs` - CIDR blocks allowed to access RDS (⚠️ **required for security**)

Example `terraform.tfvars`:

```hcl
aws_region       = "us-east-1"
project_name     = "amprenta"
db_name          = "amprenta"
db_username      = "amprenta_user"
db_password      = "your-secure-password-here"
db_allowed_cidrs = ["10.0.0.0/16"]  # Restrict to your VPC
environment      = "dev"
```

### 3. Plan Deployment

Preview what will be created:

```bash
terraform plan
```

### 4. Apply Configuration

Deploy the infrastructure:

```bash
terraform apply
```

Type `yes` when prompted to confirm.

---

## Variable Reference

| Variable | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `aws_region` | string | No | `us-east-1` | AWS region for all resources |
| `project_name` | string | No | `amprenta` | Prefix for resource names |
| `db_name` | string | **Yes** | - | PostgreSQL database name |
| `db_username` | string | **Yes** | - | RDS master username (sensitive) |
| `db_password` | string | **Yes** | - | RDS master password (sensitive) |
| `db_instance_class` | string | No | `db.t3.micro` | RDS instance type |
| `db_allowed_cidrs` | list(string) | **Yes** | - | CIDR blocks allowed to access RDS |
| `environment` | string | No | `dev` | Environment tag (`dev`, `staging`, `prod`) |
| `lightsail_blueprint_id` | string | No | `ubuntu_22_04` | Lightsail OS image |
| `lightsail_bundle_id` | string | No | `nano_2_0` | Lightsail instance size |
| `allowed_ssh_cidrs` | list(string) | No | `["0.0.0.0/0"]` | CIDR blocks allowed for SSH |

---

## Outputs Reference

After successful apply, Terraform provides these outputs:

| Output | Description | Sensitive |
|--------|-------------|-----------|
| `rds_endpoint` | RDS PostgreSQL hostname | No |
| `rds_port` | RDS port (usually 5432) | No |
| `lightsail_ip` | Public IP of Lightsail instance | No |
| `connection_string` | Full PostgreSQL connection string | Yes |
| `environment` | Deployment environment | No |

View outputs anytime:

```bash
terraform output
terraform output rds_endpoint
terraform output -json
```

View sensitive outputs:

```bash
terraform output connection_string
```

---

## Security Notes

### ⚠️ Required: Configure `db_allowed_cidrs`

**Do not use `0.0.0.0/0`** for `db_allowed_cidrs`. Restrict database access to:
- Your VPC CIDR block
- Specific application server IPs
- VPN/bastion host ranges

Example (restrict to VPC):

```hcl
db_allowed_cidrs = ["10.0.0.0/16"]
```

### Best Practices

1. **Use AWS Secrets Manager** for production credentials instead of `terraform.tfvars`
2. **Restrict SSH access** by setting `allowed_ssh_cidrs` to your IP or VPN range
3. **Use strong passwords** (20+ characters, mixed case, symbols)
4. **Enable encryption** - RDS encryption at rest is configured by default
5. **Review security groups** - Check `rds.tf` and `lightsail.tf` for firewall rules

---

## Environment Usage

The `environment` variable affects tagging and resource naming:

### Development (`dev`)

```hcl
environment = "dev"
```

- Lower-cost instance types
- Relaxed retention policies
- Quick iteration

### Production (`prod`)

```hcl
environment      = "prod"
db_instance_class = "db.t3.medium"  # Larger instance
```

- Higher availability
- Stricter security
- Backup retention enabled

Adjust `db_instance_class` and `lightsail_bundle_id` based on environment needs.

---

## Cleanup Instructions

### Destroy All Resources

```bash
terraform destroy
```

Type `yes` when prompted.

### Destroy Specific Resources

```bash
# Preview what will be destroyed
terraform plan -destroy

# Destroy only RDS
terraform destroy -target=aws_db_instance.postgres

# Destroy only Lightsail
terraform destroy -target=aws_lightsail_instance.app
```

### Important Notes

- **RDS deletion protection** may prevent immediate deletion in production
- **Snapshots** are not automatically deleted (check AWS Console)
- **State file** contains sensitive data - store securely (e.g., S3 backend with encryption)

---

## Troubleshooting

### Error: "No valid credential sources found"

AWS credentials not configured. Run:

```bash
aws configure
```

### Error: "db_allowed_cidrs must be set"

Add to `terraform.tfvars`:

```hcl
db_allowed_cidrs = ["10.0.0.0/16"]
```

### Connection Issues After Apply

1. Check security group rules allow your IP
2. Verify `db_allowed_cidrs` includes your source IP
3. Test connection:

```bash
psql $(terraform output -raw connection_string)
```

---

## Advanced Usage

### Remote State (Recommended for Teams)

Configure S3 backend in `main.tf`:

```hcl
terraform {
  backend "s3" {
    bucket = "my-terraform-state"
    key    = "amprenta/terraform.tfstate"
    region = "us-east-1"
    encrypt = true
  }
}
```

### Import Existing Resources

```bash
terraform import aws_db_instance.postgres my-existing-rds
```

---

## Related Documentation

- [AWS Lightsail Pricing](https://aws.amazon.com/lightsail/pricing/)
- [RDS PostgreSQL Guide](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/CHAP_PostgreSQL.html)
- [Terraform AWS Provider](https://registry.terraform.io/providers/hashicorp/aws/latest/docs)

---

## Support

For issues or questions:
- Check [Troubleshooting](#troubleshooting) section
- Review AWS CloudWatch logs
- Validate `terraform plan` output before applying

