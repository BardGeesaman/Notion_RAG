"""add aromatic_rings to compounds"""

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "aa1b2c3d4e5f"
down_revision = "z6u7v8w9x0y1"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("compounds", sa.Column("aromatic_rings", sa.Integer(), nullable=True))


def downgrade() -> None:
    op.drop_column("compounds", "aromatic_rings")

