"""Add companies table and attach users to companies.

Revision ID: e1f2g3h4i5j6
Revises: d4e5f6a7b8c9
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "e1f2g3h4i5j6"
down_revision = "d4e5f6a7b8c9"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "companies",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("subdomain", sa.String(length=100), nullable=False),
        sa.Column("logo_url", sa.String(length=500), nullable=True),
        sa.Column("primary_color", sa.String(length=50), nullable=True),
        sa.Column("max_users", sa.Integer(), nullable=True),
        sa.Column("max_datasets", sa.Integer(), nullable=True),
        sa.Column("storage_quota_gb", sa.Float(), nullable=True),
        sa.Column("status", sa.String(length=50), nullable=True, server_default=sa.text("'active'")),
        sa.Column("trial_ends_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
        sa.Column("updated_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
        sa.UniqueConstraint("subdomain", name="uq_companies_subdomain"),
    )
    op.create_index("ix_companies_name", "companies", ["name"])
    op.create_index("ix_companies_subdomain", "companies", ["subdomain"])
    op.create_index("ix_companies_status", "companies", ["status"])

    op.add_column("users", sa.Column("company_id", UUID(as_uuid=True), sa.ForeignKey("companies.id"), nullable=True))
    op.add_column("users", sa.Column("company_role", sa.String(length=20), nullable=True))
    op.create_index("ix_users_company_id", "users", ["company_id"])


def downgrade() -> None:
    op.drop_index("ix_users_company_id", table_name="users")
    op.drop_column("users", "company_role")
    op.drop_column("users", "company_id")

    op.drop_index("ix_companies_status", table_name="companies")
    op.drop_index("ix_companies_subdomain", table_name="companies")
    op.drop_index("ix_companies_name", table_name="companies")
    op.drop_table("companies")


