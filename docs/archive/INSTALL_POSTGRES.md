# Installing PostgreSQL on macOS

**Guide to install and configure PostgreSQL for TIER 3 architecture**

## What's Already Installed

✅ **Python packages** (already done):
- `psycopg2-binary` - Python package to connect to Postgres
- `FastAPI` - API framework
- `uvicorn` - API server

❌ **PostgreSQL server** (needs to be installed):
- The actual database server
- Database management tools

## Installation Options

### Option 1: Install via Homebrew (Recommended for macOS)

**If you have Homebrew installed:**

```bash
# Install PostgreSQL
brew install postgresql@15

# Start PostgreSQL service
brew services start postgresql@15

# Create a database
createdb amprenta

# Set a password for the postgres user (optional, but recommended)
psql postgres -c "ALTER USER postgres WITH PASSWORD 'your_password';"
```

**Default credentials after installation:**
- Username: `postgres` (or your macOS username)
- Password: Usually no password by default (or your macOS password)
- Host: `localhost`
- Port: `5432`

### Option 2: Install PostgreSQL.app (GUI Option)

1. Download from: https://postgresapp.com/
2. Install the app
3. Click "Initialize" to create a new server
4. Default settings:
   - Port: `5432`
   - Username: your macOS username
   - Password: none (by default)

### Option 3: Use Docker (No Installation Required)

If you have Docker installed, you can run Postgres in a container:

```bash
# Run Postgres in Docker
docker run --name amprenta-postgres \
  -e POSTGRES_PASSWORD=postgres \
  -e POSTGRES_DB=amprenta \
  -p 5432:5432 \
  -d postgres:15

# Default credentials:
# Username: postgres
# Password: postgres
# Host: localhost
# Port: 5432
```

## Finding Your Password

### If Postgres is Already Installed

**Common default scenarios:**

1. **No password** (common on macOS):
   - Try leaving `POSTGRES_PASSWORD` empty in `.env`
   - Or try your macOS user password

2. **Your macOS username as password**:
   - Some installations use your macOS username

3. **Reset the password**:
   ```bash
   psql postgres
   ALTER USER postgres WITH PASSWORD 'new_password';
   ```

### Check Current Postgres User

```bash
# Try to connect (will prompt for password or use trust authentication)
psql postgres

# Or check current user
whoami  # This might be your Postgres username
```

## Quick Setup Guide

### Step 1: Check if Postgres is Installed

```bash
which psql
psql --version
```

If these fail, Postgres is not installed.

### Step 2: Install Postgres (choose one method above)

### Step 3: Verify Installation

```bash
# Check if Postgres is running
pg_isready

# Or try connecting
psql postgres
```

### Step 4: Create Database

```bash
createdb amprenta
```

### Step 5: Configure `.env` File

Add to your `.env` file:

```bash
# If using default (no password)
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta
POSTGRES_USER=postgres
POSTGRES_PASSWORD=

# OR if you set a password
POSTGRES_PASSWORD=your_password_here
```

### Step 6: Test Connection

```bash
python scripts/validate_postgres_setup.py
```

## Troubleshooting

### "Command not found: psql"

**Solution**: Postgres is not installed. Install using one of the methods above.

### "Connection refused"

**Solution**: Postgres server is not running.
```bash
# Start Postgres service
brew services start postgresql@15

# Or if using PostgreSQL.app, make sure it's running
```

### "Password authentication failed"

**Solution**: 
- Check if password is required
- Try empty password: `POSTGRES_PASSWORD=`
- Try your macOS username/password
- Reset password (see above)

### "Database does not exist"

**Solution**: Create the database:
```bash
createdb amprenta
```

## Recommended Setup for Development

### Using Homebrew (Easiest)

```bash
# 1. Install
brew install postgresql@15

# 2. Start service
brew services start postgresql@15

# 3. Create database
createdb amprenta

# 4. Configure .env (use no password for local dev)
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta
POSTGRES_USER=postgres
POSTGRES_PASSWORD=
```

### Using Docker (No System Installation)

```bash
# 1. Start Postgres container
docker run --name amprenta-postgres \
  -e POSTGRES_PASSWORD=postgres \
  -e POSTGRES_DB=amprenta \
  -p 5432:5432 \
  -d postgres:15

# 2. Configure .env
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta
POSTGRES_USER=postgres
POSTGRES_PASSWORD=postgres

# 3. Stop container when done
docker stop amprenta-postgres
```

## Alternative: Skip Postgres for Now

If you don't want to install Postgres right now, you can:

1. **Continue using Notion as database** (current setup works)
2. **Test non-Postgres features** (most tests work without Postgres)
3. **Install Postgres later** when you're ready for TIER 3 migration

The system is designed to work with or without Postgres!

## Next Steps

After installing and configuring Postgres:

1. **Verify connection**: `python scripts/validate_postgres_setup.py`
2. **Run migrations**: `python scripts/migrate_database.py`
3. **Run tests**: `python scripts/run_tier3_tests.py`

## Questions?

- Check if Postgres is installed: `which psql`
- Check if it's running: `pg_isready`
- Test connection: `psql postgres`
- View help: `psql --help`

