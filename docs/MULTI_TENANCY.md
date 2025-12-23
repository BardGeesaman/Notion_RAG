# Multi-Tenancy (Companies + Postgres RLS)

## Overview (RLS Architecture)

Amprenta uses Postgres Row-Level Security (RLS) to isolate tenant data by `company_id`.

Key components:

1. **Tenant column**: each tenant-scoped table has `company_id UUID NULL REFERENCES companies(id)`
2. **Request-scoped DB context**: API requests set:

```sql
SET LOCAL app.current_company_id = '<uuid>';
```

3. **RLS policies**: each tenant table has an isolation policy:

```sql
CREATE POLICY tenant_isolation ON <table>
USING (company_id = current_setting('app.current_company_id', true)::uuid);
```

When RLS is enabled, queries from a session only “see” rows where `company_id` matches the current company id setting.

## Company setup

1. Run migrations for:
   - companies table + users.company_id / users.company_role
   - company_id columns on tenant tables
   - RLS policies

2. Seed a default company (for existing installations):

```bash
python scripts/seed_default_company.py
```

3. Backfill tenant data:

```bash
python scripts/backfill_company_data.py
```

By default, the backfill:

- Uses `{table}.created_by_id` (or `created_by` / `uploaded_by`) to copy `users.company_id`
- Otherwise uses `{table}.dataset_id -> datasets.company_id`
- Otherwise falls back to the default company

## Company user management (API)

All endpoints below are **superadmin only** (global `User.role == "admin"`) in the MVP.

Auth context is header-based in this MVP:

- `X-User-Id: <uuid>`

### Create a company

```bash
curl -X POST "http://localhost:8000/api/companies" \
  -H "Content-Type: application/json" \
  -H "X-User-Id: <admin_user_uuid>" \
  -d '{"name":"Acme","subdomain":"acme"}'
```

### Get / update a company

```bash
curl -X GET "http://localhost:8000/api/companies/<company_id>" \
  -H "X-User-Id: <admin_user_uuid>"

curl -X PATCH "http://localhost:8000/api/companies/<company_id>" \
  -H "Content-Type: application/json" \
  -H "X-User-Id: <admin_user_uuid>" \
  -d '{"primary_color":"#2E86AB"}'
```

### List company users

```bash
curl -X GET "http://localhost:8000/api/companies/<company_id>/users" \
  -H "X-User-Id: <admin_user_uuid>"
```

### Invite / assign a user

```bash
curl -X POST "http://localhost:8000/api/companies/<company_id>/users/invite" \
  -H "Content-Type: application/json" \
  -H "X-User-Id: <admin_user_uuid>" \
  -d '{"email":"user@example.com","role":"member"}'
```

### Change a user’s tenant role

```bash
curl -X PATCH "http://localhost:8000/api/companies/<company_id>/users/<user_id>" \
  -H "Content-Type: application/json" \
  -H "X-User-Id: <admin_user_uuid>" \
  -d '{"role":"admin"}'
```

### Remove user from company

```bash
curl -X DELETE "http://localhost:8000/api/companies/<company_id>/users/<user_id>" \
  -H "X-User-Id: <admin_user_uuid>"
```

## Dashboard

The Streamlit page **Company Settings** calls the Companies API and expects:

- `DASHBOARD_USER_ID` env var set to an **admin** user UUID, so it can send `X-User-Id`.


