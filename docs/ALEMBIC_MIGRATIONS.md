# Alembic Database Migrations

This document explains how to use Alembic for database migrations in the multi-omics platform.

## Overview

Alembic is used to manage database schema changes as we evolve from Notion-as-database to Postgres-as-database (TIER 3 architecture evolution).

## Setup

Alembic is configured with:
- `alembic.ini` - Main configuration file
- `alembic/env.py` - Environment configuration that connects to our database
- `alembic/versions/` - Directory containing migration scripts

## Creating Migrations

### Initial Migration

To create the initial migration from your SQLAlchemy models:

```bash
# Option 1: Use the helper script
python scripts/create_initial_migration.py

# Option 2: Use alembic directly
alembic revision --autogenerate -m "Initial schema"
```

This will:
1. Compare your SQLAlchemy models to the current database state
2. Generate a migration script with all the table definitions
3. Save it to `alembic/versions/`

### Subsequent Migrations

After modifying your SQLAlchemy models, create a new migration:

```bash
alembic revision --autogenerate -m "Description of changes"
```

**Important**: Always review the generated migration script before applying it!

## Applying Migrations

### Apply All Pending Migrations

```bash
# Option 1: Use the helper script
python scripts/migrate_database.py

# Option 2: Use alembic directly
alembic upgrade head
```

### Apply to Specific Revision

```bash
alembic upgrade <revision_id>
```

### Dry Run (View History)

```bash
python scripts/migrate_database.py --dry-run
```

## Migration Workflow

1. **Modify SQLAlchemy Models**: Update models in `amprenta_rag/database/models.py`

2. **Generate Migration**:
   ```bash
   alembic revision --autogenerate -m "Your change description"
   ```

3. **Review Migration**: Check the generated file in `alembic/versions/` to ensure it's correct

4. **Apply Migration**:
   ```bash
   alembic upgrade head
   ```

5. **Verify**: Check that tables were created/modified correctly

## Common Commands

```bash
# View migration history
alembic history

# View current database revision
alembic current

# Show SQL that would be executed (without running it)
alembic upgrade head --sql

# Downgrade one revision
alembic downgrade -1

# Downgrade to specific revision
alembic downgrade <revision_id>
```

## Configuration

Database connection is configured in `.env`:

```bash
POSTGRES_URL=postgresql://user:password@localhost:5432/amprenta_rag
# OR individual components:
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta
POSTGRES_USER=postgres
POSTGRES_PASSWORD=your_password
```

The connection is automatically read by `alembic/env.py` from `amprenta_rag/config.py`.

## Troubleshooting

### Migration Conflicts

If you have conflicting migrations, you may need to:
1. Check current revision: `alembic current`
2. View history: `alembic history`
3. Resolve conflicts manually or reset (be careful!)

### Database Connection Issues

Make sure:
- Postgres is running
- Connection credentials in `.env` are correct
- Database exists (create it first if needed)

### Models Not Detected

If Alembic doesn't detect your model changes:
1. Make sure models are imported in `alembic/env.py`
2. Check that models inherit from `Base`
3. Verify model definitions are correct

## Best Practices

1. **Always review generated migrations** before applying
2. **Test migrations** on a development database first
3. **Keep migrations small and focused** - one logical change per migration
4. **Never edit applied migrations** - create a new one instead
5. **Use descriptive migration messages** - they become your change log

## Next Steps

After setting up migrations:
1. Create the initial migration with `create_initial_migration.py`
2. Review the generated migration file
3. Apply it to create the database schema
4. Begin using Postgres as the database backend

