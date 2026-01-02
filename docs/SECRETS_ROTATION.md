# Secrets Rotation Runbook

## Overview

This document provides procedures for rotating secrets used by Amprenta RAG.

## Rotation Schedule

| Secret Type | Rotation Frequency | Trigger Events |
|-------------|-------------------|----------------|
| Database passwords | Quarterly | Team changes, suspected breach |
| API keys (OpenAI, Zotero) | Annually | Key compromise, billing changes |
| JWT/HMAC keys | Quarterly | Security audit, team changes |
| SSH keys (deployment) | Annually | Team changes, key compromise |

## Rotation Procedures

### Database Password (POSTGRES_PASSWORD)

**Local Development:**
1. Update password in PostgreSQL
2. Update `.env` file
3. Restart services

**AWS Production:**
1. Update in AWS Secrets Manager (`amprenta/prod/database`)
2. ECS tasks will pick up new value on next deployment
3. Force redeployment: `aws ecs update-service --force-new-deployment`

### API Keys (OPENAI_API_KEY, ZOTERO_API_KEY)

1. Generate new key in provider dashboard
2. Update in AWS Secrets Manager (`amprenta/{env}/api-keys`)
3. Verify new key works: `curl -H "Authorization: Bearer $NEW_KEY" ...`
4. Revoke old key in provider dashboard
5. Force ECS redeployment

### JWT/HMAC Keys (JWT_SECRET_KEY, SIGNATURE_SECRET_KEY)

**Zero-Downtime Rotation:**
1. Generate new key: `python -c "import secrets; print(secrets.token_hex(32))"`
2. Update in AWS Secrets Manager
3. Deploy new version (accepts both old and new tokens during transition)
4. After 24h grace period, remove old key acceptance

**Simple Rotation (with brief downtime):**
1. Generate new key
2. Update in AWS Secrets Manager
3. Force ECS redeployment
4. Note: Existing sessions will be invalidated

### SSH Keys (LIGHTSAIL_SSH_KEY)

1. Generate new key pair on secure workstation
2. Add new public key to Lightsail instance
3. Update `LIGHTSAIL_SSH_KEY` in GitHub Secrets
4. Test deployment with manual workflow trigger
5. Remove old public key from Lightsail instance

## Emergency Rotation (Suspected Breach)

**Immediate Actions:**
1. Identify compromised secrets
2. Generate new values immediately
3. Update all locations (AWS Secrets Manager, GitHub Secrets, local .env)
4. Force redeployment of all services
5. Revoke old credentials at source (API providers, database)
6. Review access logs for unauthorized usage
7. Document incident in security log

**Post-Incident:**
1. Conduct root cause analysis
2. Update access controls if needed
3. Consider enabling AWS CloudTrail if not already
4. Review and tighten IAM permissions

## Verification

After any rotation, verify:
- [ ] API health check passes: `curl http://localhost:8000/health`
- [ ] Dashboard loads without errors
- [ ] Database queries succeed
- [ ] LLM features work (if OPENAI_API_KEY rotated)
- [ ] Celery workers connect successfully

## AWS Secrets Manager Rotation (Future)

For automated rotation, AWS Lambda can be configured:
- Database: Use AWS RDS Secrets Manager integration
- API keys: Custom Lambda for each provider
- See: https://docs.aws.amazon.com/secretsmanager/latest/userguide/rotating-secrets.html
