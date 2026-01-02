# GitHub Secrets Configuration

This document describes the GitHub secrets required for CI/CD workflows.

## Required Secrets

| Secret Name | Description | Where Used |
|-------------|-------------|------------|
| `OPENAI_API_KEY` | OpenAI API key for LLM tests | Integration tests |
| `POSTGRES_PASSWORD` | Test database password | Integration tests |
| `SIGNATURE_SECRET_KEY` | HMAC key for signatures | All tests |
| `JWT_SECRET_KEY` | JWT signing key | Auth tests |

## CD Deployment Secrets

These secrets are required for the CD workflow to deploy to AWS Lightsail:

| Secret Name | Description | Where Used |
|-------------|-------------|------------|
| `LIGHTSAIL_HOST` | Lightsail instance IP or hostname | CD workflow SSH connection |
| `LIGHTSAIL_USER` | SSH username (default: ubuntu) | CD workflow SSH connection |
| `LIGHTSAIL_SSH_KEY` | Private SSH key for deployment | CD workflow authentication |

### Setting Up CD Secrets

1. Generate SSH key pair (if not exists):
   ```bash
   ssh-keygen -t ed25519 -C "github-actions-deploy" -f deploy_key
   ```

2. Add public key to Lightsail instance:
   ```bash
   cat deploy_key.pub >> ~/.ssh/authorized_keys
   ```

3. Add private key to GitHub:
   - Go to Repository → Settings → Secrets → Actions
   - Create `LIGHTSAIL_SSH_KEY` with contents of `deploy_key`

4. Add host and user:
   - `LIGHTSAIL_HOST`: Your Lightsail instance public IP
   - `LIGHTSAIL_USER`: `ubuntu` (default for Ubuntu instances)

## Setting Up GitHub Secrets

1. Navigate to repository Settings → Secrets and variables → Actions
2. Click "New repository secret"
3. Add each required secret

## Generating Secure Keys

For cryptographic secrets, generate with:

```bash
python -c "import secrets; print(secrets.token_hex(32))"
```

## CI/CD Workflow Integration

Secrets are injected into workflows via:

```yaml
env:
  OPENAI_API_KEY: ${{ secrets.OPENAI_API_KEY }}
  SIGNATURE_SECRET_KEY: ${{ secrets.SIGNATURE_SECRET_KEY }}
  JWT_SECRET_KEY: ${{ secrets.JWT_SECRET_KEY }}
```

## Security Best Practices

1. **Never log secrets** - GitHub automatically masks them, but avoid `echo $SECRET`
2. **Rotate regularly** - Update secrets quarterly or after team changes
3. **Minimize scope** - Use environment-specific secrets when possible
4. **Audit access** - Review who has access to repository secrets

## Environment Variables vs Secrets

| Type | Use For | Example |
|------|---------|---------|
| Secrets | API keys, passwords, tokens | `OPENAI_API_KEY` |
| Variables | Non-sensitive config | `ENVIRONMENT=test` |

## Troubleshooting

### "Secret not found" errors
- Verify secret name matches exactly (case-sensitive)
- Check secret is available to the workflow environment

### "Permission denied" on Secrets Manager
- For AWS integration tests, ensure IAM role has `secretsmanager:GetSecretValue`