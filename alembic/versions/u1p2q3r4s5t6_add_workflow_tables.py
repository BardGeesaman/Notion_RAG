"""Add workflow automation tables

Revision ID: u1p2q3r4s5t6
Revises: t0o1p2q3r4s5
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON

revision = 'u1p2q3r4s5t6'
down_revision = 't0o1p2q3r4s5'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'workflow_rules',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('name', sa.String(255), nullable=False),
        sa.Column('description', sa.Text(), nullable=True),
        sa.Column('trigger_type', sa.String(100), nullable=False),
        sa.Column('trigger_config', JSON, nullable=True),
        sa.Column('action_type', sa.String(100), nullable=False),
        sa.Column('action_config', JSON, nullable=True),
        sa.Column('is_active', sa.Boolean(), default=True),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index('ix_workflow_rules_trigger_type', 'workflow_rules', ['trigger_type'])
    op.create_index('ix_workflow_rules_is_active', 'workflow_rules', ['is_active'])

    op.create_table(
        'workflow_executions',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('rule_id', UUID(as_uuid=True), sa.ForeignKey('workflow_rules.id'), nullable=False),
        sa.Column('trigger_context', JSON, nullable=True),
        sa.Column('status', sa.String(50), nullable=False, default='pending'),
        sa.Column('result', JSON, nullable=True),
        sa.Column('triggered_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
        sa.Column('completed_at', sa.DateTime(timezone=True), nullable=True),
    )
    op.create_index('ix_workflow_executions_rule_id', 'workflow_executions', ['rule_id'])
    op.create_index('ix_workflow_executions_status', 'workflow_executions', ['status'])


def downgrade():
    op.drop_index('ix_workflow_executions_status')
    op.drop_index('ix_workflow_executions_rule_id')
    op.drop_table('workflow_executions')
    op.drop_index('ix_workflow_rules_is_active')
    op.drop_index('ix_workflow_rules_trigger_type')
    op.drop_table('workflow_rules')
