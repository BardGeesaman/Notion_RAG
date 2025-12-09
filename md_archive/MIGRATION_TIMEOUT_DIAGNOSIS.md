# Migration Timeout Diagnosis

## Problem

The database migration command (`alembic upgrade head`) timed out or failed with "Target database is not up to date" error.

## Root Cause Analysis

### 1. Migration State
- **Current DB Revision**: `80008f433a68` (initial schema)
- **Pending Migration 1**: `0c9c72e35979` (Dataset/Experiment fields) - **NOT APPLIED**
- **Pending Migration 2**: `XXXX` (Signature fields) - **NOT APPLIED**

### 2. Migration Dependency Chain
```
<base> 
  → 80008f433a68 (Initial schema) ✅ APPLIED
    → 0c9c72e35979 (Dataset/Experiment fields) ⏳ PENDING
      → XXXX (Signature fields) ⏳ PENDING
```

### 3. Why It Failed

**Error: "Target database is not up to date"**

This happens when:
1. Alembic detects that the database is at an older revision than what's in the migration files
2. We tried to create a new migration (`XXXX`) before applying the previous one (`0c9c72e35979`)
3. Alembic requires migrations to be applied sequentially

**Timeout Possible Causes:**
1. Database connection was slow/unresponsive during schema comparison
2. Alembic's autogenerate was comparing large schemas
3. Network latency or database load

## Solution

### Step 1: Check Database Connection

Test basic connectivity:
```bash
python -c "from amprenta_rag.database.base import get_engine; engine = get_engine(); conn = engine.connect(); print('✅ Connected'); conn.close()"
```

### Step 2: Check Current Migration State

```bash
alembic current
# Should show: 80008f433a68
```

### Step 3: Apply Pending Migrations Sequentially

**Migration 1: Dataset/Experiment Fields**
```bash
# This should work now
alembic upgrade 0c9c72e35979
```

**Migration 2: Signature Fields**
After Migration 1 completes:
```bash
# Generate the Signature migration properly
alembic revision --autogenerate -m "Add signature metadata fields"

# Then apply it
alembic upgrade head
```

### Step 4: Verify

```bash
alembic current
# Should show the latest revision
```

## Common Issues & Solutions

### Issue 1: Database Not Running
**Symptoms**: Connection timeout, connection refused

**Solution**: 
- Check if Postgres is running: `pg_isready` or `psql -l`
- Start Postgres if needed
- Verify connection settings in `.env`

### Issue 2: Database Connection Pool Exhausted
**Symptoms**: Hangs during connection

**Solution**:
- Check for hanging connections: `SELECT * FROM pg_stat_activity;`
- Increase connection timeout in config
- Use connection pooling

### Issue 3: Large Schema Comparison
**Symptoms**: Very slow autogenerate, timeout

**Solution**:
- Add timeout to Alembic config
- Run migrations during low-load periods
- Manually create migration if autogenerate is too slow

### Issue 4: Permission Issues
**Symptoms**: Permission denied errors

**Solution**:
- Check database user permissions
- Ensure user can CREATE/ALTER tables
- Grant necessary privileges

## Quick Fix Commands

### Test Connection
```bash
python -c "from amprenta_rag.database.base import get_engine; engine = get_engine(); print('✅ Engine created')"
```

### Check Migration Status
```bash
alembic current
alembic history
```

### Apply Migrations (Safe)
```bash
# Apply one at a time
alembic upgrade 0c9c72e35979  # Dataset/Experiment fields
alembic upgrade head          # Signature fields (after generating)
```

### Dry Run (See SQL)
```bash
alembic upgrade head --sql > migration.sql
```

## Prevention

1. **Always check current revision** before creating new migrations
2. **Apply migrations immediately** after creating them
3. **Test connection** before running migrations
4. **Use dry-run** to preview changes
5. **Apply migrations sequentially** in development

## Recommended Workflow

```bash
# 1. Check status
alembic current

# 2. Apply pending migrations
alembic upgrade head

# 3. Create new migration (if needed)
alembic revision --autogenerate -m "Description"

# 4. Review generated migration file
# Edit if needed

# 5. Apply new migration
alembic upgrade head

# 6. Verify
alembic current
```

## Next Steps

1. ✅ Database connection works (verified)
2. ⏳ Apply Migration 1: `alembic upgrade 0c9c72e35979`
3. ⏳ Generate Migration 2 properly: `alembic revision --autogenerate -m "Add signature metadata fields"`
4. ⏳ Apply Migration 2: `alembic upgrade head`
5. ⏳ Verify all migrations applied

