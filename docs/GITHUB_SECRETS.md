# GitHub Secrets Configuration

## Overview

This document outlines the GitHub repository secrets required for CI/CD pipelines, testing, and deployment. GitHub Secrets provide secure storage for sensitive data used in GitHub Actions workflows.

## Required Secrets for CI/CD

Configure these secrets in GitHub repository settings:
**Settings → Secrets and variables → Actions → New repository secret**

### Required for AWS Deployment

| Secret Name | Description | Example | Notes |
|-------------|-------------|---------|-------|
| `AWS_ACCESS_KEY_ID` | AWS IAM access key for deployment | `AKIA...` | Dedicated GitHub Actions IAM user |
| `AWS_SECRET_ACCESS_KEY` | AWS IAM secret key | `wJal...` | Rotate quarterly |
| `AWS_REGION` | AWS region for deployment | `us-east-1` | Must match Terraform configuration |

**IAM Policy for GitHub Actions User:**
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "ecr:GetAuthorizationToken",
        "ecr:BatchCheckLayerAvailability",
        "ecr:GetDownloadUrlForLayer",
        "ecr:BatchGetImage",
        "ecr:InitiateLayerUpload",
        "ecr:UploadLayerPart",
        "ecr:CompleteLayerUpload",
        "ecr:PutImage",
        "ecs:UpdateService",
        "ecs:DescribeServices",
        "ecs:DescribeTasks",
        "ecs:DescribeTaskDefinition",
        "ecs:RegisterTaskDefinition"
      ],
      "Resource": "*"
    }
  ]
}
```

### Required for Testing (Environment-Specific)

| Secret Name | Description | Required For | Can Mock? |
|-------------|-------------|--------------|-----------|
| `TEST_DATABASE_URL` | PostgreSQL URL for integration tests | Database tests | No - use test DB |
| `OPENAI_API_KEY` | OpenAI API key | LLM integration tests | Yes - recommended |
| `POSTGRES_PASSWORD` | Test database password | Test database setup | No - test-specific |

**Test Database Configuration:**
```bash
# Example TEST_DATABASE_URL format
postgresql://test_user:test_password@localhost:5432/amprenta_test
```

### Optional Secrets (Enhanced CI/CD)

| Secret Name | Description | Required For | Priority |
|-------------|-------------|--------------|----------|
| `CODECOV_TOKEN` | Codecov upload token | Coverage reports | Medium |
| `SLACK_WEBHOOK_URL` | Slack notifications | Deploy alerts | Low |
| `SONAR_TOKEN` | SonarQube analysis token | Code quality | Medium |
| `DOCKER_HUB_USERNAME` | Docker Hub username | Public image push | Low |
| `DOCKER_HUB_TOKEN` | Docker Hub access token | Public image push | Low |

## Workflow Usage Examples

### Basic Test Workflow
```yaml
name: Test
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    
    services:
      postgres:
        image: postgres:15
        env:
          POSTGRES_PASSWORD: test_password
          POSTGRES_DB: amprenta_test
        options: >-
          --health-cmd pg_isready
          --health-interval 10s
          --health-timeout 5s
          --health-retries 5
    
    env:
      # Use GitHub secrets for external services
      OPENAI_API_KEY: ${{ secrets.OPENAI_API_KEY }}
      
      # Use service configuration for test database
      TEST_DATABASE_URL: postgresql://postgres:test_password@localhost:5432/amprenta_test
      
      # Test-specific configuration
      TESTING: true
      MOCK_EXTERNAL_APIS: true
    
    steps:
      - uses: actions/checkout@v4
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'
      
      - name: Install dependencies
        run: |
          pip install -r requirements.txt
          pip install -r requirements-test.txt
      
      - name: Run tests
        run: |
          pytest --cov=amprenta_rag --cov-report=xml
      
      - name: Upload coverage
        if: env.CODECOV_TOKEN
        uses: codecov/codecov-action@v3
        with:
          token: ${{ secrets.CODECOV_TOKEN }}
```

### Deployment Workflow
```yaml
name: Deploy
on:
  push:
    branches: [main]
  workflow_dispatch:

jobs:
  deploy:
    runs-on: ubuntu-latest
    if: github.ref == 'refs/heads/main'
    
    environment: production  # Requires manual approval
    
    steps:
      - uses: actions/checkout@v4
      
      - name: Configure AWS credentials
        uses: aws-actions/configure-aws-credentials@v4
        with:
          aws-access-key-id: ${{ secrets.AWS_ACCESS_KEY_ID }}
          aws-secret-access-key: ${{ secrets.AWS_SECRET_ACCESS_KEY }}
          aws-region: ${{ secrets.AWS_REGION }}
      
      - name: Login to Amazon ECR
        uses: aws-actions/amazon-ecr-login@v2
      
      - name: Build and push Docker image
        env:
          ECR_REGISTRY: ${{ steps.login-ecr.outputs.registry }}
          ECR_REPOSITORY: amprenta-rag-api
          IMAGE_TAG: ${{ github.sha }}
        run: |
          docker build -t $ECR_REGISTRY/$ECR_REPOSITORY:$IMAGE_TAG .
          docker push $ECR_REGISTRY/$ECR_REPOSITORY:$IMAGE_TAG
      
      - name: Deploy to ECS
        run: |
          aws ecs update-service \
            --cluster amprenta-prod \
            --service amprenta-api \
            --force-new-deployment
```

### Secret Scanning Workflow
```yaml
name: Security
on: [push, pull_request]

jobs:
  secret-scan:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0  # Full history for comprehensive scanning
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'
      
      - name: Install detect-secrets
        run: pip install detect-secrets
      
      - name: Check for secrets
        run: |
          detect-secrets scan --baseline .secrets.baseline
      
      - name: Audit secrets baseline
        run: |
          detect-secrets audit .secrets.baseline
```

## Security Best Practices

### 1. Least Privilege IAM

**Create dedicated GitHub Actions IAM user:**
```bash
# AWS CLI commands to create dedicated user
aws iam create-user --user-name github-actions-amprenta
aws iam attach-user-policy --user-name github-actions-amprenta --policy-arn arn:aws:iam::ACCOUNT:policy/GitHubActionsPolicy
aws iam create-access-key --user-name github-actions-amprenta
```

**Permissions checklist:**
- ✅ ECR: Push/pull container images
- ✅ ECS: Update services and task definitions
- ❌ IAM: No user management permissions
- ❌ S3: No broad bucket access
- ❌ Secrets Manager: No access (handled by ECS)

### 2. Secret Rotation

**Quarterly rotation schedule:**
- **AWS Keys**: Rotate every 90 days
- **API Keys**: Check provider recommendations
- **Database Passwords**: Rotate during maintenance windows
- **Tokens**: Regenerate when team members leave

**Rotation procedure:**
1. Generate new credentials
2. Update GitHub secrets
3. Test deployment pipeline
4. Revoke old credentials
5. Document rotation date

### 3. Access Control

**Branch Protection Rules:**
```yaml
# .github/branch-protection.yml (if using probot/settings)
branches:
  main:
    protection:
      required_status_checks:
        strict: true
        contexts:
          - test
          - secret-scan
      enforce_admins: true
      required_pull_request_reviews:
        required_approving_review_count: 2
        dismiss_stale_reviews: true
        restrict_pushes: true
```

**Environment Protection:**
- **Production**: Require manual approval
- **Staging**: Automatic deployment from main
- **Development**: Automatic deployment from feature branches

### 4. Audit and Monitoring

**Regular audit checklist:**
- [ ] Review Actions logs for secret exposure
- [ ] Check for hardcoded credentials in workflows
- [ ] Verify least privilege IAM permissions
- [ ] Confirm secret rotation dates
- [ ] Test emergency credential revocation

**Monitoring alerts:**
- Failed deployments due to credential issues
- Unusual API usage patterns
- New workflow files added
- Changes to secret configurations

## Environment-Specific Configurations

### Development Environment
```yaml
# .github/workflows/deploy-dev.yml
env:
  ENVIRONMENT: dev
  AWS_REGION: us-east-1
  ECR_REPOSITORY: amprenta-rag-dev
```

### Staging Environment
```yaml
# .github/workflows/deploy-staging.yml
env:
  ENVIRONMENT: staging
  AWS_REGION: us-east-1
  ECR_REPOSITORY: amprenta-rag-staging
```

### Production Environment
```yaml
# .github/workflows/deploy-prod.yml
env:
  ENVIRONMENT: prod
  AWS_REGION: us-east-1
  ECR_REPOSITORY: amprenta-rag-prod

# Requires manual approval
environment: production
```

## Testing Without Real Secrets

### Mock External APIs
```python
# In test files
import pytest
from unittest.mock import patch, MagicMock

@pytest.fixture
def mock_openai():
    with patch('openai.ChatCompletion.create') as mock:
        mock.return_value = MagicMock(
            choices=[MagicMock(message=MagicMock(content="Test response"))]
        )
        yield mock

def test_llm_integration(mock_openai):
    # Test runs without real OpenAI API key
    result = my_llm_function("test input")
    assert result == "Test response"
    mock_openai.assert_called_once()
```

### Conditional Test Execution
```python
import pytest
import os

@pytest.mark.skipif(
    not os.getenv("OPENAI_API_KEY"), 
    reason="OpenAI API key not available"
)
def test_real_openai_integration():
    # Only runs if real API key is available
    pass

@pytest.mark.integration
def test_database_integration():
    # Requires TEST_DATABASE_URL
    if not os.getenv("TEST_DATABASE_URL"):
        pytest.skip("Test database not configured")
```

## Troubleshooting

### Common Issues

#### 1. "Secret not found" in Actions
**Symptoms:** Workflow fails with empty environment variables
**Solutions:**
- Verify secret name matches exactly (case-sensitive)
- Check repository vs organization level secrets
- Ensure workflow has permission to access secrets

#### 2. AWS Authentication Failures
**Symptoms:** "Unable to locate credentials" in deployment
**Solutions:**
- Verify AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY are set
- Check IAM user permissions
- Ensure AWS_REGION matches deployment region

#### 3. Database Connection Failures in Tests
**Symptoms:** "Connection refused" to PostgreSQL
**Solutions:**
- Verify PostgreSQL service is configured in workflow
- Check TEST_DATABASE_URL format
- Ensure database service is healthy before tests

### Debug Commands

**Test secret availability:**
```yaml
- name: Debug secrets
  run: |
    echo "AWS Region: ${{ secrets.AWS_REGION }}"
    echo "Has OpenAI key: ${{ secrets.OPENAI_API_KEY != '' }}"
    echo "Has DB URL: ${{ secrets.TEST_DATABASE_URL != '' }}"
```

**Test AWS credentials:**
```yaml
- name: Test AWS access
  run: |
    aws sts get-caller-identity
    aws ecr describe-repositories --repository-names amprenta-rag-api
```

## Migration from Existing Setup

### From Hardcoded Values
1. **Identify hardcoded credentials** in workflow files
2. **Create GitHub secrets** for each credential
3. **Update workflows** to use `${{ secrets.SECRET_NAME }}`
4. **Test thoroughly** in staging environment
5. **Remove hardcoded values** from repository

### From Environment Variables
1. **Review current .env usage** in CI
2. **Map to GitHub secrets** where appropriate
3. **Keep non-sensitive config** as environment variables
4. **Update documentation** for new secret requirements

## Security Compliance

### SOC 2 / ISO 27001 Requirements
- **Access Logging**: GitHub provides audit logs for secret access
- **Encryption**: Secrets encrypted at rest and in transit
- **Access Control**: Role-based access via repository permissions
- **Rotation**: Documented quarterly rotation procedures

### GDPR Considerations
- **No PII in secrets**: Use UUIDs or tokens, not personal data
- **Data Residency**: Consider GitHub's data center locations
- **Access Rights**: Team members can request secret access logs

---

## Quick Reference

### Setup Checklist
- [ ] Create GitHub Actions IAM user with minimal permissions
- [ ] Add required secrets to GitHub repository settings
- [ ] Configure branch protection rules
- [ ] Set up environment-specific deployments
- [ ] Test workflows with mock data
- [ ] Document secret rotation schedule

### Emergency Procedures
1. **Suspected secret leak**: Immediately rotate all credentials
2. **Deployment failure**: Check AWS credentials and permissions
3. **Test failures**: Verify database and API key availability
4. **Access issues**: Review IAM policies and GitHub permissions

### Regular Maintenance
- **Monthly**: Review Actions logs for anomalies
- **Quarterly**: Rotate AWS credentials and API keys
- **Annually**: Audit IAM permissions and update documentation
