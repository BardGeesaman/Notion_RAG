# Terraform AWS ECS Infrastructure

Terraform configuration for deploying Amprenta RAG infrastructure to AWS using ECS Fargate.

**Infrastructure Overview:**
- **VPC** with public and private subnets across 2 availability zones
- **Application Load Balancer (ALB)** (internet-facing) with path-based routing
- **ECS Fargate Cluster** with API and Dashboard services
- **RDS PostgreSQL** in private subnets with automated backups
- **AWS Secrets Manager** for secure credential management
- **CloudWatch Logs** for centralized logging
- **ECR Repositories** for Docker images

---

## Prerequisites

- **Terraform** >= 1.0 ([Install](https://www.terraform.io/downloads))
- **AWS CLI** configured with valid credentials ([Setup](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html))
- **Docker** for building and pushing images ([Install](https://docs.docker.com/get-docker/))
- AWS account with appropriate permissions (ECS, RDS, VPC, ALB, Secrets Manager, ECR)

Verify your setup:

```bash
terraform version
aws sts get-caller-identity
docker --version
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
- **Note**: `db_password` is managed by AWS Secrets Manager (not set in variables)
- `vpc_cidr` - VPC CIDR block (e.g., `10.0.0.0/16`)

Example `terraform.tfvars`:

```hcl
aws_region       = "us-east-1"
project_name     = "amprenta"
environment      = "dev"

# VPC Configuration
vpc_cidr         = "10.0.0.0/16"

# Database Configuration
db_name          = "amprenta"
db_username      = "amprenta_user"
# db_password is managed by AWS Secrets Manager - populate after terraform apply
db_instance_class = "db.t3.micro"

# ECS Service Configuration
desired_count    = 1
api_cpu          = 512
api_memory       = 1024
dashboard_cpu    = 512
dashboard_memory = 1024
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

### 5. Build and Push Docker Images

After infrastructure is deployed, build and push your Docker images to ECR:

```bash
# Get ECR repository URLs from Terraform outputs
API_REPO=$(terraform output -raw api_ecr_repository_url)
DASHBOARD_REPO=$(terraform output -raw dashboard_ecr_repository_url)

# Authenticate Docker to ECR
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin $API_REPO

# Build and push API image
docker build -t amprenta-api -f deploy/docker/api/Dockerfile .
docker tag amprenta-api:latest $API_REPO:latest
docker push $API_REPO:latest

# Build and push Dashboard image
docker build -t amprenta-dashboard -f deploy/docker/dashboard/Dockerfile .
docker tag amprenta-dashboard:latest $DASHBOARD_REPO:latest
docker push $DASHBOARD_REPO:latest

# Trigger ECS service update
aws ecs update-service --cluster amprenta-dev-cluster --service amprenta-dev-api --force-new-deployment
aws ecs update-service --cluster amprenta-dev-cluster --service amprenta-dev-dashboard --force-new-deployment
```

### 6. Update Secrets Manager Values

After deployment, update API keys and other secrets:

```bash
# Update OpenAI API key
aws secretsmanager put-secret-value \
  --secret-id amprenta-dev-openai-key \
  --secret-string "your-openai-api-key"

# Update Pinecone API key
aws secretsmanager put-secret-value \
  --secret-id amprenta-dev-pinecone-key \
  --secret-string "your-pinecone-api-key"
```

Or use the AWS Console: Services → Secrets Manager → Select secret → Retrieve secret value → Edit

### 7. Access the Application

Get the ALB DNS name:

```bash
terraform output alb_dns_name
```

Access your application:
- **Dashboard**: `http://<alb_dns_name>/`
- **API**: `http://<alb_dns_name>/api/`
- **API Docs**: `http://<alb_dns_name>/api/docs`

---

## Variable Reference

| Variable | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `aws_region` | string | No | `us-east-1` | AWS region for all resources |
| `project_name` | string | No | `amprenta` | Prefix for resource names |
| `environment` | string | No | `dev` | Environment tag (`dev`, `staging`, `prod`) |
| `vpc_cidr` | string | **Yes** | - | VPC CIDR block (e.g., `10.0.0.0/16`) |
| `db_name` | string | **Yes** | - | PostgreSQL database name |
| `db_username` | string | **Yes** | - | RDS master username (sensitive) |
| ~~`db_password`~~ | ~~string~~ | **Removed** | - | **Now managed by AWS Secrets Manager** |
| `db_instance_class` | string | No | `db.t3.micro` | RDS instance type |
| `desired_count` | number | No | `1` | Number of ECS tasks to run per service |
| `api_cpu` | number | No | `512` | CPU units for API container (1024 = 1 vCPU) |
| `api_memory` | number | No | `1024` | Memory (MB) for API container |
| `dashboard_cpu` | number | No | `512` | CPU units for Dashboard container |
| `dashboard_memory` | number | No | `1024` | Memory (MB) for Dashboard container |

---

## Outputs Reference

After successful apply, Terraform provides these outputs:

| Output | Description | Sensitive |
|--------|-------------|-----------|
| `vpc_id` | VPC ID | No |
| `alb_dns_name` | ALB public DNS name | No |
| `api_ecr_repository_url` | ECR repository URL for API images | No |
| `dashboard_ecr_repository_url` | ECR repository URL for Dashboard images | No |
| `rds_endpoint` | RDS PostgreSQL hostname | No |
| `rds_port` | RDS port (usually 5432) | No |
| `ecs_cluster_name` | ECS cluster name | No |
| `api_service_name` | API ECS service name | No |
| `dashboard_service_name` | Dashboard ECS service name | No |

View outputs anytime:

```bash
terraform output
terraform output alb_dns_name
terraform output -json
```

---

## Architecture Diagram

```
Internet
    │
    ▼
[Application Load Balancer]
    │
    ├─► /api/*  ──────► [ECS Fargate: API Service]
    │                       │
    │                       ├─► [API Task 1] (Private Subnet)
    │                       └─► [API Task 2] (Private Subnet)
    │
    └─► /*      ──────► [ECS Fargate: Dashboard Service]
                            │
                            ├─► [Dashboard Task 1] (Private Subnet)
                            └─► [Dashboard Task 2] (Private Subnet)

[RDS PostgreSQL] (Private Subnet)
    ▲
    │
    └─────────────── [Both Services]

[AWS Secrets Manager]
    ▲
    │
    └─────────────── [Both Services]
```

---

## Secret Management

### Secrets Created Automatically

Terraform creates these secrets in AWS Secrets Manager:
- `{project_name}-{environment}-db-url` - Database connection string (auto-populated)
- `{project_name}-{environment}-openai-key` - OpenAI API key (empty, must set manually)
- `{project_name}-{environment}-pinecone-key` - Pinecone API key (empty, must set manually)

### Setting Secret Values

**Via AWS CLI:**

```bash
aws secretsmanager put-secret-value \
  --secret-id amprenta-dev-openai-key \
  --secret-string "sk-..."
```

**Via AWS Console:**
1. Navigate to: Services → Secrets Manager
2. Click on the secret name
3. Click "Retrieve secret value"
4. Click "Edit"
5. Enter new value
6. Click "Save"

### Using Secrets in Code

ECS tasks automatically inject secrets as environment variables:
- `DATABASE_URL` - Full PostgreSQL connection string
- `OPENAI_API_KEY` - OpenAI API key
- `PINECONE_API_KEY` - Pinecone API key

No code changes needed - environment variables are available at runtime.

---

## Security Notes

### Network Security

- **Public Subnets**: ALB only
- **Private Subnets**: ECS tasks and RDS (no direct internet access)
- **NAT Gateway**: Allows private subnets to reach internet (for Docker pulls, API calls)
- **Security Groups**: Restrict traffic to minimum required ports

### Database Security

- RDS is in **private subnets only** (not internet-accessible)
- Security group allows connections only from ECS tasks
- Encryption at rest enabled by default
- Automated backups enabled

### Best Practices

1. **Use strong passwords** (20+ characters, mixed case, symbols)
2. **Rotate secrets regularly** via Secrets Manager rotation
3. **Enable MFA** for AWS Console access
4. **Review CloudWatch logs** for suspicious activity
5. **Use AWS WAF** with ALB for production (not included in this config)
6. **Enable GuardDuty** for threat detection

---

## Scaling

### Horizontal Scaling (More Tasks)

Update `desired_count` in `terraform.tfvars`:

```hcl
desired_count = 3  # Run 3 tasks per service
```

Apply changes:

```bash
terraform apply
```

### Vertical Scaling (Larger Tasks)

Update CPU and memory in `terraform.tfvars`:

```hcl
api_cpu       = 1024  # 1 vCPU
api_memory    = 2048  # 2 GB
dashboard_cpu = 1024
dashboard_memory = 2048
```

Apply changes:

```bash
terraform apply
```

### Auto-Scaling (Future Enhancement)

Add ECS auto-scaling based on:
- CPU utilization
- Memory utilization
- Request count (ALB metrics)

---

## Monitoring and Logging

### CloudWatch Logs

View logs for each service:

```bash
# API logs
aws logs tail /ecs/amprenta-dev-api --follow

# Dashboard logs
aws logs tail /ecs/amprenta-dev-dashboard --follow
```

### CloudWatch Metrics

Monitor in AWS Console:
- ECS → Clusters → [cluster-name] → Metrics
- Key metrics: CPU, Memory, Network, Task count

### Alarms (Recommended)

Create CloudWatch Alarms for:
- High CPU usage (> 80%)
- High memory usage (> 80%)
- Task failure rate
- ALB 5xx errors

---

## Deployment Updates

### Update Application Code

1. Build and push new Docker images:

```bash
docker build -t amprenta-api:v2 -f deploy/docker/api/Dockerfile .
docker tag amprenta-api:v2 $API_REPO:latest
docker push $API_REPO:latest
```

2. Force ECS service update:

```bash
aws ecs update-service \
  --cluster amprenta-dev-cluster \
  --service amprenta-dev-api \
  --force-new-deployment
```

### Update Infrastructure

1. Modify `terraform.tfvars` or `*.tf` files
2. Run `terraform plan` to preview changes
3. Run `terraform apply` to apply changes

---

## Cleanup Instructions

### Destroy All Resources

```bash
terraform destroy
```

Type `yes` when prompted.

**Important**: This will delete:
- All ECS services and tasks
- RDS database (data will be lost unless snapshots exist)
- VPC and networking components
- ALB
- ECR repositories (images will be deleted)
- Secrets Manager secrets

### Selective Cleanup

```bash
# Preview what will be destroyed
terraform plan -destroy

# Destroy only ECS resources
terraform destroy -target=aws_ecs_cluster.main
```

---

## Troubleshooting

### Error: "No valid credential sources found"

AWS credentials not configured. Run:

```bash
aws configure
```

### ECS Tasks Failing to Start

1. Check CloudWatch logs:

```bash
aws logs tail /ecs/amprenta-dev-api --follow
```

2. Common issues:
   - Missing secrets (OpenAI, Pinecone keys not set)
   - Image pull errors (ECR authentication failed)
   - Database connection errors (security group misconfigured)

### Cannot Connect to Database

1. Verify security groups allow ECS → RDS traffic
2. Check RDS is in private subnets
3. Verify DATABASE_URL secret is correct

### ALB Health Check Failures

1. Ensure API/Dashboard containers expose correct ports
2. Verify health check path returns 200 OK
3. Check target group settings in AWS Console

### Docker Push Fails

Re-authenticate to ECR:

```bash
aws ecr get-login-password --region us-east-1 | \
  docker login --username AWS --password-stdin \
  $(terraform output -raw api_ecr_repository_url | cut -d'/' -f1)
```

---

## Cost Estimation

**Monthly costs** (approximate, us-east-1):

| Resource | Configuration | Cost |
|----------|---------------|------|
| ECS Fargate (API) | 1 task, 0.5 vCPU, 1GB RAM | ~$15 |
| ECS Fargate (Dashboard) | 1 task, 0.5 vCPU, 1GB RAM | ~$15 |
| RDS PostgreSQL | db.t3.micro | ~$15 |
| ALB | 1 ALB, low traffic | ~$20 |
| NAT Gateway | 1 NAT, low data transfer | ~$35 |
| CloudWatch Logs | 10 GB/month | ~$5 |
| **Total** | | **~$105/month** |

**Cost optimization tips:**
- Use Fargate Spot for non-production
- Reduce `desired_count` during off-hours
- Use RDS Reserved Instances for production
- Enable S3 VPC endpoint to avoid NAT costs

---

## Related Documentation

- [AWS ECS Fargate Guide](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/AWS_Fargate.html)
- [RDS PostgreSQL Guide](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/CHAP_PostgreSQL.html)
- [Terraform AWS Provider](https://registry.terraform.io/providers/hashicorp/aws/latest/docs)
- [AWS Secrets Manager](https://docs.aws.amazon.com/secretsmanager/latest/userguide/intro.html)

---

## Support

For issues or questions:
- Check [Troubleshooting](#troubleshooting) section
- Review AWS CloudWatch logs
- Validate `terraform plan` output before applying
- Check ECS service events in AWS Console
