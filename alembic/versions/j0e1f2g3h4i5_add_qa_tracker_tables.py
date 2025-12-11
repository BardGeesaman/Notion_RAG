"""Add QA tracker tables

Revision ID: j0e1f2g3h4i5
Revises: i9d0e1f2g3h4
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON

revision = "j0e1f2g3h4i5"
down_revision = "i9d0e1f2g3h4"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "saved_questions",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("question_text", sa.Text, nullable=False),
        sa.Column("tags", JSON, nullable=True),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("created_at", sa.DateTime, nullable=False),
        sa.Column("updated_at", sa.DateTime, nullable=False),
        sa.Column("is_archived", sa.Boolean, nullable=False, server_default=sa.text("false")),
    )

    op.create_table(
        "saved_answers",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("question_id", UUID(as_uuid=True), sa.ForeignKey("saved_questions.id"), nullable=False),
        sa.Column("answer_text", sa.Text, nullable=False),
        sa.Column("evidence", JSON, nullable=True),
        sa.Column("model_used", sa.String(200), nullable=True),
        sa.Column("version", sa.Integer, nullable=False, server_default=sa.text("1")),
        sa.Column("confidence_score", sa.Float, nullable=True),
        sa.Column("created_at", sa.DateTime, nullable=False),
    )

    op.create_index("ix_saved_answers_question_id", "saved_answers", ["question_id"])


def downgrade() -> None:
    op.drop_index("ix_saved_answers_question_id")
    op.drop_table("saved_answers")
    op.drop_table("saved_questions")
