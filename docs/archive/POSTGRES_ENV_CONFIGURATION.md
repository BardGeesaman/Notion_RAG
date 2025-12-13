# Postgres Connection Configuration Guide

**How to configure Postgres database connection in `.env` file**

## Overview

The TIER 3 infrastructure supports two ways to configure Postgres connection:
1. **Individual components** (recommended for clarity)
2. **Full connection URL** (convenient for some setups)

## Configuration Options

### Option 1: Individual Components (Recommended)

Add these variables to your `.env` file:

```bash
# Postgres Database Configuration
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=postgres
POSTGRES_PASSWORD=your_password_here
POSTGRES_ECHO=false
```

**Variables Explained**:
- `POSTGRES_HOST` - Database server hostname (default: `localhost`)
- `POSTGRES_PORT` - Database server port (default: `5432`)
- `POSTGRES_DB` - Database name (default: `amprenta_rag`)
- `POSTGRES_USER` - Database username (default: `postgres`)
- `POSTGRES_PASSWORD` - Database password (**required**)
- `POSTGRES_ECHO` - Enable SQL query logging (default: `false`, set to `true` for debugging)

### Option 2: Full Connection URL

Alternatively, you can use a single connection URL:

```bash
# Postgres Database Configuration (Full URL)
POSTGRES_URL=postgresql://username:password@hostname:port/database_name

# Example:
POSTGRES_URL=postgresql://postgres:mypassword@localhost:5432/amprenta_rag
```

**Note**: If `POSTGRES_URL` is set, individual component variables are ignored.

## Step-by-Step Setup

### 1. Create or Edit `.env` File

Navigate to your project root directory and create/edit the `.env` file:

```bash
cd "/Users/bard/Documents/Notion RAG"
nano .env  # or use your preferred editor
```

### 2. Add Postgres Configuration

Add the following to your `.env` file (replace with your actual values):

```bash
# ============================================
# Postgres Database Configuration (TIER 3)
# ============================================

# Option 1: Individual components (recommended)
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=postgres
POSTGRES_PASSWORD=your_actual_password_here

# Optional: Enable SQL logging for debugging
# POSTGRES_ECHO=false

# OR

# Option 2: Full connection URL (alternative)
# POSTGRES_URL=postgresql://postgres:password@localhost:5432/amprenta_rag
```

### 3. Example Configurations

#### Local Development (Default)
```bash
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=postgres
POSTGRES_PASSWORD=postgres
```

#### Remote Server
```bash
POSTGRES_HOST=db.example.com
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=myuser
POSTGRES_PASSWORD=secure_password
```

#### Using Full URL
```bash
POSTGRES_URL=postgresql://myuser:secure_password@db.example.com:5432/amprenta_rag
```

#### Docker Compose
```bash
POSTGRES_HOST=postgres  # Service name in docker-compose.yml
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=postgres
POSTGRES_PASSWORD=postgres
```

#### Cloud Database (e.g., AWS RDS, Heroku)
```bash
# Heroku-style
POSTGRES_URL=postgresql://user:pass@host.amazonaws.com:5432/dbname

# Or with SSL
POSTGRES_URL=postgresql://user:pass@host.amazonaws.com:5432/dbname?sslmode=require
```

## Security Best Practices

### ⚠️ Important Security Notes

1. **Never commit `.env` to Git**
   - `.env` should be in `.gitignore`
   - Use `.env.example` for documentation (without real passwords)

2. **Use Strong Passwords**
   - Use complex passwords for production databases
   - Consider using password managers

3. **Environment-Specific Files**
   - Use `.env.local` for local development
   - Use `.env.production` for production (with proper access controls)

4. **Password in URL**
   - If using `POSTGRES_URL`, the password will be in the URL string
   - Be extra careful with file permissions

## Verifying Configuration

### 1. Test Configuration Loading

```bash
python -c "from amprenta_rag.config import get_config; cfg = get_config().postgres; print(f'Host: {cfg.host}, DB: {cfg.db}, User: {cfg.user}')"
```

### 2. Test Database Connection

```bash
# This will test if the connection works
python scripts/validate_postgres_setup.py
```

### 3. Check Environment Variables

```bash
# Check if variables are loaded
python -c "import os; from dotenv import load_dotenv; load_dotenv(); print('POSTGRES_HOST:', os.getenv('POSTGRES_HOST', 'NOT SET'))"
```

## Troubleshooting

### Issue: "Connection refused"

**Solution**: 
- Check if Postgres is running: `pg_isready` or `psql -h localhost -U postgres`
- Verify host and port are correct
- Check firewall settings

### Issue: "Authentication failed"

**Solution**:
- Verify username and password are correct
- Check `pg_hba.conf` configuration
- Ensure user has proper permissions

### Issue: "Database does not exist"

**Solution**:
- Create the database: `createdb amprenta_rag`
- Or use an existing database name

### Issue: Environment variables not loading

**Solution**:
- Ensure `.env` file is in project root
- Check file permissions
- Verify variable names are correct (case-sensitive)
- Restart your terminal/IDE

## Complete .env Example

Here's a complete `.env` example with all configurations:

```bash
# ============================================
# Core API Keys (Required)
# ============================================
OPENAI_API_KEY=sk-your-key-here
PINECONE_API_KEY=your-key-here
NOTION_API_KEY=secret_your-key-here
ZOTERO_API_KEY=your-key-here

# ============================================
# Postgres Database Configuration (TIER 3)
# ============================================
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=postgres
POSTGRES_PASSWORD=your_password_here
POSTGRES_ECHO=false

# Alternative: Full URL (uncomment to use instead)
# POSTGRES_URL=postgresql://postgres:password@localhost:5432/amprenta_rag

# ============================================
# Notion Database IDs (Required)
# ============================================
NOTION_LIT_DB_ID=your-db-id-here
NOTION_EMAIL_DB_ID=your-db-id-here
NOTION_RAG_DB_ID=your-db-id-here
NOTION_EXP_DATA_DB_ID=your-db-id-here
NOTION_SIGNATURE_DB_ID=your-db-id-here
NOTION_PROGRAMS_DB_ID=your-db-id-here
NOTION_EXPERIMENTS_DB_ID=your-db-id-here

# ============================================
# Optional Configuration
# ============================================
PINECONE_INDEX_NAME=amprenta-rag
PINECONE_NAMESPACE=default
SIGNATURE_OVERLAP_THRESHOLD=0.3
ENABLE_SIGNATURE_SCORING=true
ENABLE_LIPID_MAPPING=true
```

## Next Steps

After configuring the `.env` file:

1. **Verify Configuration**
   ```bash
   python scripts/validate_postgres_setup.py
   ```

2. **Run Migrations**
   ```bash
   python scripts/migrate_database.py
   ```

3. **Run Tests**
   ```bash
   python scripts/run_tier3_tests.py
   ```

## Quick Reference

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `POSTGRES_HOST` | No | `localhost` | Database hostname |
| `POSTGRES_PORT` | No | `5432` | Database port |
| `POSTGRES_DB` | No | `amprenta_rag` | Database name |
| `POSTGRES_USER` | No | `postgres` | Database username |
| `POSTGRES_PASSWORD` | **Yes** | - | Database password |
| `POSTGRES_ECHO` | No | `false` | Enable SQL logging |
| `POSTGRES_URL` | No | - | Full connection URL (overrides above) |

## Additional Resources

- See `docs/ALEMBIC_MIGRATIONS.md` for migration setup
- See `docs/POSTGRES_SOT_TRANSITION.md` for transition guide
- See `scripts/validate_postgres_setup.py` for validation script

