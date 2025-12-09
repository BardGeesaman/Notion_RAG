# Migration Timeout - Root Cause & Fix

## 🔍 Root Cause Identified

The migration is timing out because of **blocking database connections**:

### Problem Connections
1. **PID 65365**: `idle in transaction` - **BLOCKING OTHERS**
   - State: idle in transaction since 15:13:39
   - Query: SELECT features... (left open)
   - **This is the main blocker**

2. **PID 77970**: Blocked by PID 65365
3. **PID 79250**: Blocked by PID 65365

### Why This Happens

When a database connection is left in "idle in transaction" state:
- It holds locks on tables
- Other operations (like migrations) cannot acquire the same locks
- The migration times out waiting for locks to be released

This typically happens when:
- A previous script didn't properly close its database session
- An exception occurred during a transaction
- A connection pool didn't clean up properly

## ✅ Solution

### Option 1: Terminate Blocking Connections (Recommended)

Run the diagnostic script to see what's blocking:
```bash
python scripts/fix_database_locks.py --check
```

Terminate idle connections:
```bash
python scripts/fix_database_locks.py --terminate-idle
```

Or terminate a specific PID:
```bash
python scripts/fix_database_locks.py --terminate-pid 65365
```

### Option 2: Manual SQL Fix

Connect to Postgres and run:
```sql
-- Check idle connections
SELECT pid, usename, state, query_start, state_change, query
FROM pg_stat_activity
WHERE datname = current_database()
AND state = 'idle in transaction'
AND pid != pg_backend_pid();

-- Terminate specific connection
SELECT pg_terminate_backend(65365);

-- Or terminate all idle in transaction connections
SELECT pg_terminate_backend(pid)
FROM pg_stat_activity
WHERE datname = current_database()
AND state = 'idle in transaction'
AND pid != pg_backend_pid();
```

### Option 3: Restart Database (Nuclear Option)

If connections persist:
```bash
# Restart Postgres (macOS)
brew services restart postgresql

# Or if using Docker
docker restart <postgres_container>
```

## 🚀 After Fixing Locks

Once connections are cleared, apply migrations:

```bash
# Step 1: Apply Dataset/Experiment migration
alembic upgrade 0c9c72e35979

# Step 2: Generate Signature migration (after Step 1)
alembic revision --autogenerate -m "Add signature metadata fields"

# Step 3: Apply Signature migration
alembic upgrade head
```

## 🔒 Prevention

To prevent this in the future:

1. **Always close database sessions**:
   ```python
   db = next(get_db())
   try:
       # ... use db ...
   finally:
       db.close()
   ```

2. **Use context managers**:
   ```python
   from contextlib import contextmanager
   
   @contextmanager
   def get_db_session():
       db = next(get_db())
       try:
           yield db
       finally:
           db.close()
   ```

3. **Set transaction timeouts**:
   ```python
   # In database connection config
   connect_args={"options": "-c statement_timeout=30000"}  # 30 seconds
   ```

4. **Monitor connections**:
   - Regularly check for idle connections
   - Set up connection pool limits
   - Monitor connection duration

## 📋 Quick Fix Steps

1. **Identify blockers**:
   ```bash
   python scripts/fix_database_locks.py --check
   ```

2. **Clear blockers**:
   ```bash
   python scripts/fix_database_locks.py --terminate-idle
   ```

3. **Apply migrations**:
   ```bash
   alembic upgrade head
   ```

4. **Verify**:
   ```bash
   alembic current
   ```

## ⚠️ Important Notes

- **Terminating connections** will rollback any uncommitted transactions
- **Save work** before terminating connections
- **Check what's running** before terminating (may be important operations)
- **Use with caution** in production environments

