# Migration Timeout - Resolution Summary

## Problem

Database migrations were timing out due to blocking connections.

## Root Cause

**Idle in Transaction Connection (PID 65365)**
- Left a SELECT query open in a transaction
- Was blocking all ALTER TABLE operations (migrations)
- Multiple migration attempts got stuck waiting for locks

**Blocking Chain:**
```
PID 65365 (idle in transaction) 
  → Blocks PID 79250 (migration attempt)
    → Blocks PID 77970 (migration attempt)
      → Blocks PID 79337 (migration attempt)
```

## Solution Applied

1. **Terminated idle connection**: Used `pg_terminate_backend(65365)` to kill the blocking connection
2. **Cleared stuck migrations**: All migration attempts were automatically released when the blocker was killed
3. **Verified clean state**: Confirmed no more blocking connections

## Prevention

### Immediate Actions
- ✅ Created `scripts/fix_database_locks.py` to diagnose and fix locks
- ✅ Documented the issue and solution

### Long-term Prevention
1. **Always use proper session cleanup**:
   ```python
   db = next(get_db())
   try:
       # ... use db ...
   finally:
       db.close()
   ```

2. **Use context managers** for database sessions
3. **Set statement timeouts** on connections
4. **Monitor for idle connections** regularly

## Status

- ✅ Blocking connections cleared
- ⏳ Ready to apply migrations
- ✅ Database is clean

## Next Steps

1. Apply Dataset/Experiment migration:
   ```bash
   alembic upgrade 0c9c72e35979
   ```

2. Generate Signature migration:
   ```bash
   alembic revision --autogenerate -m "Add signature metadata fields"
   ```

3. Apply Signature migration:
   ```bash
   alembic upgrade head
   ```

## Tools Created

- `scripts/fix_database_locks.py` - Diagnose and fix database locks
- `docs/MIGRATION_TIMEOUT_FIX.md` - Detailed troubleshooting guide

## Quick Reference

```bash
# Check for locks
python scripts/fix_database_locks.py --check

# Fix idle connections
python scripts/fix_database_locks.py --terminate-idle

# Terminate specific PID
python scripts/fix_database_locks.py --terminate-pid <PID>
```

