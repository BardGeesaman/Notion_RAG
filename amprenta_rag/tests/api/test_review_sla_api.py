"""Unit tests for SLA management API endpoints and Celery tasks."""

from datetime import datetime, timezone, timedelta
from unittest.mock import MagicMock, patch, ANY
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.api.main import app
from amprenta_rag.database.models import ReviewSLA, ReviewCycle
from amprenta_rag.models.auth import EntityReview, User
from amprenta_rag.jobs.tasks.review_sla import check_overdue_reviews_task, process_due_cycles_task


# Test fixtures
@pytest.fixture
def mock_admin_user():
    """Mock admin user for testing."""
    return User(
        id=uuid4(),
        username="admin_user",
        email="admin@test.com",
        role="admin"
    )


@pytest.fixture
def mock_regular_user():
    """Mock regular user for testing."""
    return User(
        id=uuid4(),
        username="regular_user",
        email="user@test.com",
        role="user"
    )


@pytest.fixture
def mock_db_session():
    """Mock database session."""
    return MagicMock()


@pytest.fixture
def client():
    """Test client."""
    return TestClient(app)


class TestSLAManagementAPI:
    """Test SLA management API endpoints."""

    @pytest.fixture(autouse=True)
    def setup(self, mock_admin_user, mock_db_session):
        """Setup test dependencies."""
        app.dependency_overrides[get_current_user] = lambda: mock_admin_user
        app.dependency_overrides[get_database_session] = lambda: mock_db_session
        self.client = TestClient(app)
        yield
        app.dependency_overrides.clear()

    def test_list_sla_rules(self, mock_db_session):
        """Test GET /api/v1/sla/rules endpoint."""
        # Mock SLA rules
        mock_rules = [
            ReviewSLA(
                id=uuid4(),
                name="Dataset SLA",
                entity_type="dataset",
                max_review_hours=48,
                warning_threshold_pct=75,
                is_default=True,
                is_active=True,
                created_at=datetime.now(timezone.utc),
            ),
            ReviewSLA(
                id=uuid4(),
                name="Global SLA",
                entity_type=None,
                max_review_hours=120,
                warning_threshold_pct=80,
                is_default=True,
                is_active=True,
                created_at=datetime.now(timezone.utc),
            ),
        ]

        mock_db_session.query.return_value.order_by.return_value.all.return_value = mock_rules

        response = self.client.get("/api/v1/sla/rules")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        assert data[0]["name"] == "Dataset SLA"
        assert data[0]["entity_type"] == "dataset"
        assert data[1]["name"] == "Global SLA"
        assert data[1]["entity_type"] is None

    def test_create_sla_rule_admin_only(self, mock_db_session):
        """Test POST /api/v1/sla/rules requires admin access."""
        rule_data = {
            "name": "Test SLA",
            "entity_type": "experiment",
            "max_review_hours": 72,
            "warning_threshold_pct": 85,
            "escalation_chain": [str(uuid4())],
            "is_default": False,
            "is_active": True,
        }

        mock_sla = ReviewSLA(
            id=uuid4(),
            **rule_data,
            created_at=datetime.now(timezone.utc),
        )

        mock_db_session.add = MagicMock()
        mock_db_session.commit = MagicMock()
        mock_db_session.refresh = MagicMock()

        with patch("amprenta_rag.api.routers.review_sla.ReviewSLA") as mock_sla_class:
            mock_sla_class.return_value = mock_sla

            response = self.client.post("/api/v1/sla/rules", json=rule_data)

            assert response.status_code == 201
            data = response.json()
            assert data["name"] == "Test SLA"
            assert data["entity_type"] == "experiment"
            assert data["max_review_hours"] == 72

            mock_db_session.add.assert_called_once()
            mock_db_session.commit.assert_called_once()

    def test_list_review_cycles(self, mock_db_session):
        """Test GET /api/v1/sla/review-cycles endpoint."""
        # Mock review cycles
        mock_cycles = [
            ReviewCycle(
                id=uuid4(),
                name="Weekly Dataset Reviews",
                entity_type="dataset",
                frequency="weekly",
                day_of_week=1,
                reviewer_pool=[str(uuid4()), str(uuid4())],
                is_active=True,
                next_run_at=datetime.now(timezone.utc) + timedelta(days=7),
                created_at=datetime.now(timezone.utc),
            ),
        ]

        with patch("amprenta_rag.api.routers.review_sla.get_cycles") as mock_get_cycles:
            mock_get_cycles.return_value = mock_cycles

            response = self.client.get("/api/v1/sla/review-cycles")

            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["name"] == "Weekly Dataset Reviews"
            assert data[0]["frequency"] == "weekly"

    def test_create_cycle_admin_only(self, mock_db_session, mock_admin_user):
        """Test POST /api/v1/sla/review-cycles requires admin access."""
        cycle_data = {
            "name": "Monthly Experiment Reviews",
            "entity_type": "experiment",
            "frequency": "monthly",
            "reviewer_pool": [str(uuid4()), str(uuid4())],
            "day_of_month": 15,
        }

        mock_cycle = ReviewCycle(
            id=uuid4(),
            **cycle_data,
            created_by_id=mock_admin_user.id,
            created_at=datetime.now(timezone.utc),
            next_run_at=datetime.now(timezone.utc) + timedelta(days=30),
            is_active=True,
        )

        with patch("amprenta_rag.api.routers.review_sla.create_cycle") as mock_create_cycle:
            mock_create_cycle.return_value = mock_cycle

            response = self.client.post("/api/v1/sla/review-cycles", json=cycle_data)

            assert response.status_code == 201
            data = response.json()
            assert data["name"] == "Monthly Experiment Reviews"
            assert data["frequency"] == "monthly"

            mock_create_cycle.assert_called_once_with(
                name=cycle_data["name"],
                entity_type=cycle_data["entity_type"],
                frequency=cycle_data["frequency"],
                reviewer_pool=cycle_data["reviewer_pool"],
                created_by_id=mock_admin_user.id,
                db=ANY,
                program_id=None,
                day_of_week=None,
                day_of_month=15,
            )

    def test_get_sla_status_summary(self, mock_db_session):
        """Test GET /api/v1/sla/status endpoint."""
        # Mock active reviews with SLAs
        mock_reviews = [
            EntityReview(
                id=uuid4(),
                status="in_review",
                sla_id=uuid4(),
                due_at=datetime.now(timezone.utc) + timedelta(hours=24),
            ),
            EntityReview(
                id=uuid4(),
                status="pending",
                sla_id=uuid4(),
                due_at=datetime.now(timezone.utc) - timedelta(hours=12),
            ),
        ]

        mock_db_session.query.return_value.filter.return_value.all.return_value = mock_reviews

        with patch("amprenta_rag.api.routers.review_sla.check_sla_status") as mock_check_status:
            mock_check_status.side_effect = [
                {"status": "on_track"},
                {"status": "overdue"},
            ]

            response = self.client.get("/api/v1/sla/status")

            assert response.status_code == 200
            data = response.json()
            assert data["on_track"] == 1
            assert data["overdue"] == 1
            assert data["warning"] == 0
            assert data["breached"] == 0
            assert data["total"] == 2

    def test_get_review_sla_status(self, mock_db_session):
        """Test GET /api/v1/sla/reviews/{review_id}/sla endpoint."""
        review_id = uuid4()
        sla_id = uuid4()

        mock_review = EntityReview(
            id=review_id,
            status="in_review",
            sla_id=sla_id,
            due_at=datetime.now(timezone.utc) + timedelta(hours=24),
        )
        mock_review.sla = ReviewSLA(
            id=sla_id,
            name="Test SLA",
            max_review_hours=48,
            warning_threshold_pct=75,
        )

        mock_db_session.query.return_value.filter.return_value.first.return_value = mock_review

        with patch("amprenta_rag.api.routers.review_sla.check_sla_status") as mock_check_status:
            mock_check_status.return_value = {
                "status": "on_track",
                "hours_remaining": 24.0,
                "pct_elapsed": 50.0,
                "due_at": mock_review.due_at,
                "warning_threshold": 75,
            }

            response = self.client.get(f"/api/v1/sla/reviews/{review_id}/sla")

            assert response.status_code == 200
            data = response.json()
            assert data["review_id"] == str(review_id)
            assert data["status"] == "on_track"
            assert data["hours_remaining"] == 24.0
            assert data["pct_elapsed"] == 50.0
            assert data["sla_name"] == "Test SLA"


class TestSLACeleryTasks:
    """Test SLA-related Celery tasks."""

    @patch("amprenta_rag.jobs.tasks.review_sla.db_session")
    def test_check_overdue_reviews_task_sends_reminders(self, mock_db_session):
        """Test check_overdue_reviews_task sends reminders."""
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db

        # Mock overdue reviews
        review1 = EntityReview(
            id=uuid4(),
            status="in_review",
            reminder_sent_at=None,  # No reminder sent
        )
        review1.sla = ReviewSLA(max_review_hours=48, escalation_chain=None)

        review2 = EntityReview(
            id=uuid4(),
            status="pending",
            reminder_sent_at=datetime.now(timezone.utc) - timedelta(hours=25),  # Old reminder
        )
        review2.sla = ReviewSLA(max_review_hours=48, escalation_chain=None)

        with patch("amprenta_rag.jobs.tasks.review_sla.get_overdue_reviews") as mock_get_overdue:
            with patch("amprenta_rag.jobs.tasks.review_sla.check_sla_status") as mock_check_status:
                with patch("amprenta_rag.jobs.tasks.review_sla.send_reminder") as mock_send_reminder:
                    
                    mock_get_overdue.return_value = [review1, review2]
                    mock_check_status.return_value = {"hours_remaining": -12}  # Not severely overdue
                    mock_send_reminder.return_value = True

                    result = check_overdue_reviews_task()

                    assert result["checked"] == 2
                    assert result["reminded"] == 2
                    assert result["escalated"] == 0
                    assert result["errors"] == 0

                    assert mock_send_reminder.call_count == 2

    @patch("amprenta_rag.jobs.tasks.review_sla.db_session")
    def test_check_overdue_reviews_task_escalates(self, mock_db_session):
        """Test check_overdue_reviews_task escalates severely overdue reviews."""
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db

        # Mock severely overdue review
        review = EntityReview(
            id=uuid4(),
            status="in_review",
        )
        review.sla = ReviewSLA(
            max_review_hours=48,
            escalation_chain=[str(uuid4()), str(uuid4())],
        )

        with patch("amprenta_rag.jobs.tasks.review_sla.get_overdue_reviews") as mock_get_overdue:
            with patch("amprenta_rag.jobs.tasks.review_sla.check_sla_status") as mock_check_status:
                with patch("amprenta_rag.jobs.tasks.review_sla.escalate_review") as mock_escalate:
                    
                    mock_get_overdue.return_value = [review]
                    mock_check_status.return_value = {"hours_remaining": -36}  # Severely overdue (>1.5x SLA)
                    mock_escalate.return_value = True

                    result = check_overdue_reviews_task()

                    assert result["checked"] == 1
                    assert result["reminded"] == 0
                    assert result["escalated"] == 1
                    assert result["errors"] == 0

                    mock_escalate.assert_called_once_with(review, mock_db)

    @patch("amprenta_rag.jobs.tasks.review_sla.db_session")
    def test_process_due_cycles_task(self, mock_db_session):
        """Test process_due_cycles_task processes due cycles."""
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db

        # Mock due cycles
        cycle1 = ReviewCycle(
            id=uuid4(),
            name="Cycle 1",
            next_run_at=datetime.now(timezone.utc) - timedelta(hours=1),
        )
        cycle2 = ReviewCycle(
            id=uuid4(),
            name="Cycle 2",
            next_run_at=datetime.now(timezone.utc) - timedelta(minutes=30),
        )

        # Mock created reviews
        mock_reviews = [
            EntityReview(id=uuid4(), entity_type="dataset"),
            EntityReview(id=uuid4(), entity_type="dataset"),
        ]

        with patch("amprenta_rag.jobs.tasks.review_sla.get_due_cycles") as mock_get_due:
            with patch("amprenta_rag.jobs.tasks.review_sla.create_cycle_reviews") as mock_create_reviews:
                with patch("amprenta_rag.jobs.tasks.review_sla.advance_cycle") as mock_advance:
                    
                    mock_get_due.return_value = [cycle1, cycle2]
                    mock_create_reviews.side_effect = [mock_reviews, []]  # First cycle creates 2, second creates 0
                    mock_advance.return_value = True

                    result = process_due_cycles_task()

                    assert result["cycles_processed"] == 2
                    assert result["reviews_created"] == 2
                    assert result["errors"] == 0

                    assert mock_create_reviews.call_count == 2
                    assert mock_advance.call_count == 2

    @patch("amprenta_rag.jobs.tasks.review_sla.db_session")
    def test_process_due_cycles_advances_schedule(self, mock_db_session):
        """Test process_due_cycles_task advances cycle schedules."""
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db

        cycle = ReviewCycle(
            id=uuid4(),
            name="Weekly Cycle",
            frequency="weekly",
            next_run_at=datetime.now(timezone.utc) - timedelta(hours=1),
        )

        with patch("amprenta_rag.jobs.tasks.review_sla.get_due_cycles") as mock_get_due:
            with patch("amprenta_rag.jobs.tasks.review_sla.create_cycle_reviews") as mock_create_reviews:
                with patch("amprenta_rag.jobs.tasks.review_sla.advance_cycle") as mock_advance:
                    
                    mock_get_due.return_value = [cycle]
                    mock_create_reviews.return_value = []
                    
                    # Mock advance_cycle to update next_run_at
                    def advance_side_effect(cycle_obj, db):
                        cycle_obj.next_run_at = datetime.now(timezone.utc) + timedelta(days=7)
                        return cycle_obj
                    
                    mock_advance.side_effect = advance_side_effect

                    result = process_due_cycles_task()

                    assert result["cycles_processed"] == 1
                    mock_advance.assert_called_once_with(cycle, mock_db)
                    
                    # Verify cycle's next_run_at was updated
                    assert cycle.next_run_at > datetime.now(timezone.utc)
