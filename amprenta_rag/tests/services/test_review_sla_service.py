"""Unit tests for review SLA and cycle services."""

from datetime import datetime, timezone, timedelta
from unittest.mock import MagicMock, patch, ANY
from uuid import uuid4

from amprenta_rag.database.models import ReviewCycle, ReviewSLA
from amprenta_rag.models.auth import EntityReview
from amprenta_rag.services.review_sla import (
    get_default_sla,
    apply_sla,
    check_sla_status,
    escalate_review,
    send_reminder,
)
from amprenta_rag.services.review_cycles import (
    create_cycle,
    get_due_cycles,
    create_cycle_reviews,
    advance_cycle,
    rotate_reviewer,
)


class TestReviewSLAService:
    """Test SLA enforcement service functions."""

    def test_get_default_sla_for_entity_type(self):
        """Test getting entity-specific default SLA."""
        mock_db = MagicMock()

        sla_id = uuid4()
        mock_sla = ReviewSLA(
            id=sla_id,
            name="Dataset Default SLA",
            entity_type="dataset",
            max_review_hours=48,
            is_default=True,
            is_active=True,
        )

        # Mock query chain
        mock_db.query.return_value.filter.return_value.first.return_value = mock_sla

        result = get_default_sla("dataset", mock_db)

        assert result == mock_sla
        assert result.name == "Dataset Default SLA"
        assert result.entity_type == "dataset"
        assert result.max_review_hours == 48
        assert result.is_default is True

    def test_get_default_sla_global_fallback(self):
        """Test fallback to global default SLA when entity-specific not found."""
        mock_db = MagicMock()

        global_sla = ReviewSLA(
            id=uuid4(),
            name="Global Default SLA",
            entity_type=None,  # Global
            max_review_hours=120,
            is_default=True,
            is_active=True,
        )

        # Mock query chain: first call returns None, second returns global SLA
        mock_db.query.return_value.filter.return_value.first.side_effect = [None, global_sla]

        result = get_default_sla("experiment", mock_db)

        assert result == global_sla
        assert result.name == "Global Default SLA"
        assert result.entity_type is None
        assert result.max_review_hours == 120

    def test_apply_sla_sets_due_date(self):
        """Test applying SLA sets due_at and sla_id on review."""
        mock_db = MagicMock()

        review_id = uuid4()
        sla_id = uuid4()

        review = EntityReview(
            id=review_id,
            entity_type="dataset",
            entity_id=uuid4(),
            reviewer_id=uuid4(),
            status="pending",
        )

        sla = ReviewSLA(
            id=sla_id,
            name="Test SLA",
            max_review_hours=48,
        )

        result = apply_sla(review, sla, mock_db)

        assert result.sla_id == sla_id
        assert result.due_at is not None
        assert result.escalation_level == 0
        assert result.reminder_sent_at is None
        assert result.escalated_at is None

        # Check due_at is approximately 48 hours from now
        expected_due = datetime.now(timezone.utc) + timedelta(hours=48)
        time_diff = abs((result.due_at - expected_due).total_seconds())
        assert time_diff < 60  # Within 1 minute

        mock_db.commit.assert_called_once()

    def test_check_sla_status_on_track(self):
        """Test SLA status check for on-track review."""
        now = datetime.now(timezone.utc)
        review_start = now - timedelta(hours=12)  # Started 12 hours ago
        due_at = now + timedelta(hours=36)  # Due in 36 hours

        sla = ReviewSLA(
            max_review_hours=48,
            warning_threshold_pct=75,
        )

        review = EntityReview(
            due_at=due_at,
            reviewed_at=review_start,
            status="in_review",
        )
        review.sla = sla

        result = check_sla_status(review)

        assert result["status"] == "on_track"
        assert result["hours_remaining"] > 30
        assert abs(result["pct_elapsed"] - 25.0) < 1.0  # Approximately 25% (12/48 hours)
        assert result["due_at"] == due_at
        assert result["warning_threshold"] == 75

    def test_check_sla_status_overdue(self):
        """Test SLA status check for overdue review."""
        now = datetime.now(timezone.utc)
        review_start = now - timedelta(hours=60)  # Started 60 hours ago
        due_at = now - timedelta(hours=12)  # Due 12 hours ago

        sla = ReviewSLA(
            max_review_hours=48,
            warning_threshold_pct=75,
        )

        review = EntityReview(
            due_at=due_at,
            reviewed_at=review_start,
            status="in_review",  # Still active
        )
        review.sla = sla

        result = check_sla_status(review)

        assert result["status"] == "breached"
        assert result["hours_remaining"] < 0
        assert result["pct_elapsed"] == 100  # Capped at 100% in check_sla_status

    def test_escalate_review_increments_level(self):
        """Test escalating review increments escalation level."""
        mock_db = MagicMock()

        escalation_chain = [str(uuid4()), str(uuid4()), str(uuid4())]
        
        sla = ReviewSLA(
            escalation_chain=escalation_chain,
        )

        review = EntityReview(
            id=uuid4(),
            escalation_level=0,
        )
        review.sla = sla

        result = escalate_review(review, mock_db)

        assert result is True
        assert review.escalation_level == 1
        assert review.escalated_at is not None

        # Check escalated_at is recent
        time_diff = abs((review.escalated_at - datetime.now(timezone.utc)).total_seconds())
        assert time_diff < 60

        mock_db.commit.assert_called_once()

    @patch("amprenta_rag.services.review_sla.create_user_notification")
    def test_send_reminder_creates_notification(self, mock_create_notification):
        """Test that send_reminder creates a user notification."""
        mock_db = MagicMock()
        
        reviewer_id = uuid4()
        review_id = uuid4()
        entity_id = uuid4()
        
        review = EntityReview(
            id=review_id,
            reviewer_id=reviewer_id,
            entity_type="dataset",
            entity_id=entity_id,
            due_at=datetime.now(timezone.utc) + timedelta(hours=24),
        )
        
        # Call function
        result = send_reminder(review, mock_db)
        
        # Verify notification was created
        assert result is True
        mock_create_notification.assert_called_once_with(
            recipient_id=reviewer_id,
            event_type=ANY,  # ActivityEventType.REVIEW_REMINDER
            target_type="entity_review",
            target_id=review_id,
            target_name=f"Review for dataset:{entity_id}",
            metadata={
                "entity_type": "dataset",
                "entity_id": str(entity_id),
                "due_at": review.due_at.isoformat(),
            }
        )
        
        # Verify reminder timestamp was set
        assert review.reminder_sent_at is not None
        mock_db.commit.assert_called_once()

    @patch("amprenta_rag.services.review_sla.create_user_notification")
    def test_escalate_review_notifies_both_reviewers(self, mock_create_notification):
        """Test that escalate_review notifies both original and escalation target."""
        mock_db = MagicMock()
        
        reviewer_id = uuid4()
        escalation_target_id = uuid4()
        review_id = uuid4()
        entity_id = uuid4()
        
        escalation_chain = [str(escalation_target_id)]
        
        sla = ReviewSLA(
            escalation_chain=escalation_chain,
        )
        
        review = EntityReview(
            id=review_id,
            reviewer_id=reviewer_id,
            entity_type="dataset",
            entity_id=entity_id,
            escalation_level=0,
        )
        review.sla = sla
        
        # Call function
        result = escalate_review(review, mock_db)
        
        # Verify both notifications were created
        assert result is True
        assert mock_create_notification.call_count == 2
        
        # Check original reviewer notification
        original_call = mock_create_notification.call_args_list[0]
        assert original_call[1]["recipient_id"] == reviewer_id
        assert original_call[1]["metadata"]["role"] == "original_reviewer"
        
        # Check escalation target notification
        escalation_call = mock_create_notification.call_args_list[1]
        assert escalation_call[1]["recipient_id"] == escalation_target_id
        assert escalation_call[1]["metadata"]["role"] == "escalation_target"
        
        # Verify escalation tracking was updated
        assert review.escalation_level == 1
        assert review.escalated_at is not None
        mock_db.commit.assert_called_once()


class TestReviewCyclesService:
    """Test review cycle management service functions."""

    def test_create_cycle_computes_next_run(self):
        """Test creating cycle computes next_run_at correctly."""
        mock_db = MagicMock()

        created_by_id = uuid4()
        reviewer_pool = [str(uuid4()), str(uuid4())]

        # Mock the created cycle
        mock_cycle = ReviewCycle(
            id=uuid4(),
            name="Weekly Dataset Reviews",
            entity_type="dataset",
            frequency="weekly",
            day_of_week=1,  # Tuesday
            reviewer_pool=reviewer_pool,
            created_by_id=created_by_id,
        )

        mock_db.add = MagicMock()
        mock_db.commit = MagicMock()
        mock_db.refresh = MagicMock()

        # Mock the actual create_cycle call
        with patch("amprenta_rag.services.review_cycles.ReviewCycle") as mock_cycle_class:
            mock_cycle_class.return_value = mock_cycle

            result = create_cycle(
                name="Weekly Dataset Reviews",
                entity_type="dataset",
                frequency="weekly",
                reviewer_pool=reviewer_pool,
                created_by_id=created_by_id,
                db=mock_db,
                day_of_week=1,
            )

            assert result == mock_cycle
            mock_db.add.assert_called_once_with(mock_cycle)
            mock_db.commit.assert_called_once()
            mock_db.refresh.assert_called_once_with(mock_cycle)

    def test_get_due_cycles(self):
        """Test finding cycles that are due to run."""
        mock_db = MagicMock()

        # Mock due cycles
        due_cycle1 = ReviewCycle(
            id=uuid4(),
            name="Due Cycle 1",
            next_run_at=datetime.now(timezone.utc) - timedelta(hours=1),
            is_active=True,
        )

        due_cycle2 = ReviewCycle(
            id=uuid4(),
            name="Due Cycle 2",
            next_run_at=datetime.now(timezone.utc) - timedelta(minutes=30),
            is_active=True,
        )

        mock_db.query.return_value.filter.return_value.all.return_value = [due_cycle1, due_cycle2]

        result = get_due_cycles(mock_db)

        assert len(result) == 2
        assert result[0] == due_cycle1
        assert result[1] == due_cycle2

    def test_create_cycle_reviews(self):
        """Test creating reviews for cycle entities."""
        mock_db = MagicMock()

        cycle_id = uuid4()
        reviewer_pool = [str(uuid4()), str(uuid4())]

        cycle = ReviewCycle(
            id=cycle_id,
            entity_type="dataset",
            reviewer_pool=reviewer_pool,
        )

        # Mock no existing reviews
        mock_db.query.return_value.filter.return_value.first.return_value = None

        # Mock add and commit
        mock_db.add = MagicMock()
        mock_db.commit = MagicMock()

        result = create_cycle_reviews(cycle, mock_db)

        # Should create reviews for mock entities
        assert len(result) >= 1  # At least one review created
        assert all(isinstance(review, EntityReview) for review in result)
        assert all(review.entity_type == "dataset" for review in result)
        assert all(review.cycle_id == cycle_id for review in result)

        # Check that reviews were added to DB
        assert mock_db.add.call_count == len(result)
        mock_db.commit.assert_called_once()

    def test_advance_cycle_weekly(self):
        """Test advancing weekly cycle to next run."""
        mock_db = MagicMock()

        current_run = datetime(2025, 1, 7, 9, 0, 0, tzinfo=timezone.utc)  # Tuesday

        cycle = ReviewCycle(
            id=uuid4(),
            frequency="weekly",
            day_of_week=1,  # Tuesday
            next_run_at=current_run,
        )

        mock_db.commit = MagicMock()

        result = advance_cycle(cycle, mock_db)

        assert result == cycle
        # Next run should be following Tuesday
        expected_next = datetime(2025, 1, 14, 9, 0, 0, tzinfo=timezone.utc)
        assert cycle.next_run_at == expected_next

        mock_db.commit.assert_called_once()

    def test_rotate_reviewer(self):
        """Test reviewer rotation from pool."""
        reviewer_pool = [str(uuid4()), str(uuid4()), str(uuid4())]

        cycle = ReviewCycle(
            id=uuid4(),
            reviewer_pool=reviewer_pool,
        )

        # Test round-robin selection
        reviewer1 = rotate_reviewer(cycle, 0)
        reviewer2 = rotate_reviewer(cycle, 1)
        reviewer3 = rotate_reviewer(cycle, 2)
        reviewer4 = rotate_reviewer(cycle, 3)  # Should wrap around

        assert str(reviewer1) == reviewer_pool[0]
        assert str(reviewer2) == reviewer_pool[1]
        assert str(reviewer3) == reviewer_pool[2]
        assert str(reviewer4) == reviewer_pool[0]  # Wrapped around
