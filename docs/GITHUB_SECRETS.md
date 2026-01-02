# GitHub Secrets Configuration

This document describes the GitHub secrets required for CI/CD workflows.

## Required Secrets

| Secret Name | Description | Where Used |
|-------------|-------------|------------|
| `OPENAI_API_KEY` | OpenAI API key for LLM tests | Integration tests |
| `POSTGRES_PASSWORD` | Test database password | Integration tests |
| `SIGNATURE_SECRET_KEY` | HMAC key for signatures | All tests |
| `JWT_SECRET_KEY` | JWT signing key | Auth tests |

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