"""Unit tests for ReviewCycle and ReviewSLA models."""

from datetime import datetime, timezone, timedelta
from unittest.mock import MagicMock, patch
from uuid import uuid4

from amprenta_rag.database.models import ReviewCycle, ReviewSLA, Program
from amprenta_rag.models.auth import EntityReview


class TestReviewSLAModels:
    """Test ReviewCycle and ReviewSLA model functionality."""

    @patch("amprenta_rag.database.session.db_session")
    def test_review_cycle_create(self, mock_session):
        """Test ReviewCycle model creation with all fields."""
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db

        cycle_id = uuid4()
        created_by_id = uuid4()
        program_id = uuid4()
        
        # Create a ReviewCycle with all fields
        cycle = ReviewCycle(
            id=cycle_id,
            name="Weekly Dataset Reviews",
            entity_type="dataset",
            frequency="weekly",
            day_of_week=1,  # Tuesday
            day_of_month=None,
            next_run_at=datetime.now(timezone.utc) + timedelta(days=7),
            reviewer_pool=[str(uuid4()), str(uuid4())],
            is_active=True,
            program_id=program_id,
            created_by_id=created_by_id,
            created_at=datetime.now(timezone.utc),
        )

        mock_db.add(cycle)
        mock_db.commit()
        mock_db.refresh(cycle)

        assert cycle.id == cycle_id
        assert cycle.name == "Weekly Dataset Reviews"
        assert cycle.entity_type == "dataset"
        assert cycle.frequency == "weekly"
        assert cycle.day_of_week == 1
        assert cycle.day_of_month is None
        assert cycle.next_run_at is not None
        assert len(cycle.reviewer_pool) == 2
        assert cycle.is_active is True
        assert cycle.program_id == program_id
        assert cycle.created_by_id == created_by_id
        assert cycle.created_at is not None

        mock_db.add.assert_called_once_with(cycle)
        mock_db.commit.assert_called_once()
        mock_db.refresh.assert_called_once_with(cycle)

    @patch("amprenta_rag.database.session.db_session")
    def test_review_cycle_program_scoped(self, mock_session):
        """Test ReviewCycle with program_id scoping."""
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db

        program_id = uuid4()
        
        # Mock Program relationship
        mock_program = Program(id=program_id, name="ALS Research Program")
        
        cycle = ReviewCycle(
            name="Monthly Experiment Reviews",
            entity_type="experiment",
            frequency="monthly",
            day_of_month=15,
            program_id=program_id,
            is_active=True,
            created_at=datetime.now(timezone.utc),
        )

        # Mock the relationship
        cycle.program = mock_program

        mock_db.add(cycle)
        mock_db.commit()
        mock_db.refresh(cycle)

        assert cycle.entity_type == "experiment"
        assert cycle.frequency == "monthly"
        assert cycle.day_of_month == 15
        assert cycle.program_id == program_id
        assert cycle.program == mock_program
        assert cycle.program.name == "ALS Research Program"

        mock_db.add.assert_called_once_with(cycle)
        mock_db.commit.assert_called_once()

    @patch("amprenta_rag.database.session.db_session")
    def test_review_sla_create(self, mock_session):
        """Test ReviewSLA model creation with escalation chain."""
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db

        sla_id = uuid4()
        escalation_chain = [str(uuid4()), str(uuid4()), str(uuid4())]
        
        sla = ReviewSLA(
            id=sla_id,
            name="Critical Dataset SLA",
            entity_type="dataset",
            max_review_hours=48,  # 2 days
            warning_threshold_pct=80,
            escalation_chain=escalation_chain,
            is_default=False,
            is_active=True,
            created_at=datetime.now(timezone.utc),
        )

        mock_db.add(sla)
        mock_db.commit()
        mock_db.refresh(sla)

        assert sla.id == sla_id
        assert sla.name == "Critical Dataset SLA"
        assert sla.entity_type == "dataset"
        assert sla.max_review_hours == 48
        assert sla.warning_threshold_pct == 80
        assert sla.escalation_chain == escalation_chain
        assert len(sla.escalation_chain) == 3
        assert sla.is_default is False
        assert sla.is_active is True
        assert sla.created_at is not None

        mock_db.add.assert_called_once_with(sla)
        mock_db.commit.assert_called_once()
        mock_db.refresh.assert_called_once_with(sla)

    @patch("amprenta_rag.database.session.db_session")
    def test_review_sla_default(self, mock_session):
        """Test ReviewSLA as default for entity type."""
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db

        # Default SLA for all experiments
        sla = ReviewSLA(
            name="Default Experiment SLA",
            entity_type="experiment",
            max_review_hours=120,  # 5 days default
            warning_threshold_pct=75,
            escalation_chain=None,
            is_default=True,
            is_active=True,
            created_at=datetime.now(timezone.utc),
        )

        mock_db.add(sla)
        mock_db.commit()
        mock_db.refresh(sla)

        assert sla.name == "Default Experiment SLA"
        assert sla.entity_type == "experiment"
        assert sla.max_review_hours == 120
        assert sla.warning_threshold_pct == 75
        assert sla.escalation_chain is None
        assert sla.is_default is True
        assert sla.is_active is True

        mock_db.add.assert_called_once_with(sla)

    @patch("amprenta_rag.database.session.db_session")
    def test_entity_review_sla_fields(self, mock_session):
        """Test EntityReview with SLA tracking fields."""
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db

        review_id = uuid4()
        entity_id = uuid4()
        reviewer_id = uuid4()
        sla_id = uuid4()
        
        due_at = datetime.now(timezone.utc) + timedelta(hours=48)
        
        # Create EntityReview with SLA fields
        review = EntityReview(
            id=review_id,
            entity_type="dataset",
            entity_id=entity_id,
            reviewer_id=reviewer_id,
            status="pending",
            comments=None,
            reviewed_at=None,
            # SLA tracking fields
            due_at=due_at,
            sla_id=sla_id,
            cycle_id=None,
            reminder_sent_at=None,
            escalated_at=None,
            escalation_level=0,
        )

        mock_db.add(review)
        mock_db.commit()
        mock_db.refresh(review)

        assert review.id == review_id
        assert review.entity_type == "dataset"
        assert review.entity_id == entity_id
        assert review.reviewer_id == reviewer_id
        assert review.status == "pending"
        assert review.due_at == due_at
        assert review.sla_id == sla_id
        assert review.cycle_id is None
        assert review.reminder_sent_at is None
        assert review.escalated_at is None
        assert review.escalation_level == 0

        mock_db.add.assert_called_once_with(review)
        mock_db.commit.assert_called_once()

    @patch("amprenta_rag.database.session.db_session")
    def test_entity_review_cycle_relationship(self, mock_session):
        """Test EntityReview linked to ReviewCycle."""
        mock_db = MagicMock()
        mock_session.return_value.__enter__.return_value = mock_db

        review_id = uuid4()
        cycle_id = uuid4()
        sla_id = uuid4()
        entity_id = uuid4()
        reviewer_id = uuid4()
        
        # Mock ReviewCycle and ReviewSLA
        mock_cycle = ReviewCycle(
            id=cycle_id,
            name="Weekly Dataset Reviews",
            entity_type="dataset",
            frequency="weekly",
        )
        
        mock_sla = ReviewSLA(
            id=sla_id,
            name="Standard Dataset SLA",
            entity_type="dataset",
            max_review_hours=96,
        )
        
        review = EntityReview(
            id=review_id,
            entity_type="dataset",
            entity_id=entity_id,
            reviewer_id=reviewer_id,
            status="pending",
            cycle_id=cycle_id,
            sla_id=sla_id,
            due_at=datetime.now(timezone.utc) + timedelta(hours=96),
            escalation_level=0,
        )

        # Mock the relationships
        review.cycle = mock_cycle
        review.sla = mock_sla

        mock_db.add(review)
        mock_db.commit()
        mock_db.refresh(review)

        assert review.cycle_id == cycle_id
        assert review.sla_id == sla_id
        assert review.cycle == mock_cycle
        assert review.sla == mock_sla
        assert review.cycle.name == "Weekly Dataset Reviews"
        assert review.cycle.frequency == "weekly"
        assert review.sla.name == "Standard Dataset SLA"
        assert review.sla.max_review_hours == 96
        assert review.due_at is not None

        mock_db.add.assert_called_once_with(review)
        mock_db.commit.assert_called_once()
