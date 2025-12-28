"""
Unit tests for activity service.
"""

import pytest
from uuid import uuid4
from unittest.mock import Mock, patch

from amprenta_rag.services.activity import (
    log_activity,
    get_activity_feed,
    get_user_notifications,
    mark_notification_read,
    mark_all_notifications_read,
    get_unread_count,
)
from amprenta_rag.database.models import (
    ActivityEvent,
    ActivityEventType,
    RepositoryNotification,
    Program,
)
from amprenta_rag.models.auth import User
from amprenta_rag.database.session import db_session


@pytest.fixture
def mock_db_session():
    """Provide a mock database session for testing."""
    mock_session = Mock()
    mock_db = Mock()
    mock_session.return_value.__enter__.return_value = mock_db
    mock_session.return_value.__exit__.return_value = None
    return mock_session, mock_db


class TestActivityService:
    """Test activity service functions."""

    def test_log_activity_creates_event(self):
        """Test that log_activity creates an ActivityEvent in the database."""
        target_id = uuid4()
        
        # Log activity
        event = log_activity(
            event_type=ActivityEventType.COMPOUND_ADDED,
            target_type="compound",
            target_id=target_id,
            target_name="Test Compound",
            actor_id=None,
            program_id=None,
            metadata={"smiles": "CCO"}
        )
        
        # Verify event was created
        assert event is not None
        assert event.event_type == "compound_added"
        assert event.target_type == "compound"
        assert event.target_id == target_id
        assert event.target_name == "Test Compound"
        assert event.actor_id is None
        assert event.program_id is None
        assert event.event_metadata == {"smiles": "CCO"}
        
        # Verify it's in the database
        with db_session() as db:
            db_event = db.query(ActivityEvent).filter(ActivityEvent.id == event.id).first()
            assert db_event is not None
            assert db_event.event_type == "compound_added"

    def test_log_activity_with_program_routes_notification(self):
        """Test that activity with program_id creates notification for program lead."""
        # Create a user to be program lead with unique credentials
        lead_id = uuid4()
        with db_session() as db:
            user = User(
                id=lead_id,
                username=f"program_lead_{uuid4().hex[:8]}",
                email=f"lead_{uuid4().hex[:8]}@test.com",
                password_hash="dummy_hash"
            )
            db.add(user)
            db.commit()
            
            # Create a program with the lead
            program_id = uuid4()
            program = Program(
                id=program_id,
                name="Test Program",
                created_by_id=lead_id,
                description="Test program"
            )
            db.add(program)
            db.commit()
        
        # Log activity with program (actor_id=None for system action)
        event = log_activity(
            event_type=ActivityEventType.EXPERIMENT_CREATED,
            target_type="experiment",
            target_id=uuid4(),
            target_name="Test Experiment",
            actor_id=None,  # System action (or could create another user)
            program_id=program_id,
            metadata={"experiment_type": "assay"}
        )
        
        # Verify notification was created
        assert event is not None
        with db_session() as db:
            notification = db.query(RepositoryNotification).filter(
                RepositoryNotification.activity_event_id == event.id
            ).first()
            assert notification is not None
            assert notification.notification_type == "activity"
            assert notification.is_read is False

    def test_route_notifications_excludes_actor(self):
        """Test that notifications are not sent to the actor of the action."""
        # Create a user who is both actor and program lead with unique credentials
        user_id = uuid4()
        with db_session() as db:
            user = User(
                id=user_id,
                username=f"testuser_{uuid4().hex[:8]}",
                email=f"user_{uuid4().hex[:8]}@test.com",
                password_hash="dummy_hash"
            )
            db.add(user)
            db.commit()
            
            # Create program where user is lead
            program_id = uuid4()
            program = Program(
                id=program_id,
                name="Test Program",
                created_by_id=user_id,
                description="Test program"
            )
            db.add(program)
            db.commit()
        
        # Log activity where actor is the program lead
        event = log_activity(
            event_type=ActivityEventType.MODEL_TRAINED,
            target_type="model",
            target_id=uuid4(),
            target_name="Test Model",
            actor_id=user_id,  # Same as program lead
            program_id=program_id,
            metadata={"model_type": "classifier"}
        )
        
        # Verify no notification was created (actor == lead)
        assert event is not None
        with db_session() as db:
            notification_count = db.query(RepositoryNotification).filter(
                RepositoryNotification.activity_event_id == event.id
            ).count()
            assert notification_count == 0

    def test_get_activity_feed_filters(self):
        """Test that get_activity_feed correctly filters by event type."""
        # Create multiple events with different types
        compound_event = log_activity(
            event_type=ActivityEventType.COMPOUND_ADDED,
            target_type="compound",
            target_id=uuid4(),
            target_name="Compound 1",
        )
        
        experiment_event = log_activity(
            event_type=ActivityEventType.EXPERIMENT_CREATED,
            target_type="experiment",
            target_id=uuid4(),
            target_name="Experiment 1",
        )
        
        model_event = log_activity(
            event_type=ActivityEventType.MODEL_TRAINED,
            target_type="model",
            target_id=uuid4(),
            target_name="Model 1",
        )
        
        # Test filtering by event type
        compound_events = get_activity_feed(event_type="compound_added")
        experiment_events = get_activity_feed(event_type="experiment_created")
        all_events = get_activity_feed()
        
        # Verify filtering
        compound_ids = [e.id for e in compound_events]
        experiment_ids = [e.id for e in experiment_events]
        all_ids = [e.id for e in all_events]
        
        assert compound_event.id in compound_ids
        assert compound_event.id not in experiment_ids
        assert experiment_event.id in experiment_ids
        assert experiment_event.id not in compound_ids
        assert compound_event.id in all_ids
        assert experiment_event.id in all_ids
        assert model_event.id in all_ids

    def test_get_user_notifications(self):
        """Test that get_user_notifications returns correct notifications for user."""
        # Create user and program with unique credentials
        user_id = uuid4()
        with db_session() as db:
            user = User(
                id=user_id,
                username=f"testuser_{uuid4().hex[:8]}",
                email=f"user_{uuid4().hex[:8]}@test.com",
                password_hash="dummy_hash"
            )
            db.add(user)
            db.commit()
            
            program_id = uuid4()
            program = Program(
                id=program_id,
                name="Test Program",
                created_by_id=user_id,
                description="Test program"
            )
            db.add(program)
            db.commit()
        
        # Create activity that should generate notification
        event = log_activity(
            event_type=ActivityEventType.HIT_CONFIRMED,
            target_type="compound",
            target_id=uuid4(),
            target_name="Hit Compound",
            actor_id=None,  # System action
            program_id=program_id,
        )
        
        # Get notifications for user
        notifications = get_user_notifications(user_id)
        
        # Verify notification is returned
        assert len(notifications) > 0
        notification = notifications[0]
        assert notification.activity_event_id == event.id
        assert notification.notification_type == "activity"

    def test_mark_notification_read(self):
        """Test that mark_notification_read updates notification status."""
        # Setup user, program, and notification with unique credentials
        user_id = uuid4()
        with db_session() as db:
            user = User(
                id=user_id,
                username=f"testuser_{uuid4().hex[:8]}",
                email=f"user_{uuid4().hex[:8]}@test.com",
                password_hash="dummy_hash"
            )
            db.add(user)
            db.commit()
            
            program_id = uuid4()
            program = Program(
                id=program_id,
                name="Test Program",
                created_by_id=user_id,
                description="Test program"
            )
            db.add(program)
            db.commit()
        
        # Create activity and notification
        log_activity(
            event_type=ActivityEventType.STATUS_CHANGED,
            target_type="experiment",
            target_id=uuid4(),
            target_name="Test Experiment",
            actor_id=None,  # System action
            program_id=program_id,
        )
        
        # Get the notification
        notifications = get_user_notifications(user_id)
        assert len(notifications) > 0
        notification = notifications[0]
        assert notification.is_read is False
        
        # Mark as read
        success = mark_notification_read(notification.id, user_id)
        assert success is True
        
        # Verify it's marked as read
        with db_session() as db:
            updated_notification = db.query(RepositoryNotification).filter(
                RepositoryNotification.id == notification.id
            ).first()
            assert updated_notification.is_read is True

    def test_get_unread_count(self):
        """Test that get_unread_count returns correct count of unread notifications."""
        # Setup user and program with unique credentials
        user_id = uuid4()
        with db_session() as db:
            user = User(
                id=user_id,
                username=f"testuser_{uuid4().hex[:8]}",
                email=f"user_{uuid4().hex[:8]}@test.com",
                password_hash="dummy_hash"
            )
            db.add(user)
            db.commit()
            
            program_id = uuid4()
            program = Program(
                id=program_id,
                name="Test Program",
                created_by_id=user_id,
                description="Test program"
            )
            db.add(program)
            db.commit()
        
        # Get initial count
        initial_count = get_unread_count(user_id)
        
        # Create multiple activities
        log_activity(
            event_type=ActivityEventType.COMPOUND_ADDED,
            target_type="compound",
            target_id=uuid4(),
            target_name="Compound 1",
            actor_id=None,  # System action
            program_id=program_id,
        )
        
        log_activity(
            event_type=ActivityEventType.EXPERIMENT_CREATED,
            target_type="experiment",
            target_id=uuid4(),
            target_name="Experiment 1",
            actor_id=None,  # System action
            program_id=program_id,
        )
        
        # Check count increased
        new_count = get_unread_count(user_id)
        assert new_count == initial_count + 2
        
        # Mark one as read
        notifications = get_user_notifications(user_id)
        if notifications:
            mark_notification_read(notifications[0].id, user_id)
        
        # Check count decreased
        final_count = get_unread_count(user_id)
        assert final_count == initial_count + 1

    def test_mark_all_notifications_read(self):
        """Test that mark_all_notifications_read updates all user notifications."""
        # Setup user and program with unique credentials
        user_id = uuid4()
        with db_session() as db:
            user = User(
                id=user_id,
                username=f"testuser_{uuid4().hex[:8]}",
                email=f"user_{uuid4().hex[:8]}@test.com",
                password_hash="dummy_hash"
            )
            db.add(user)
            db.commit()
            
            program_id = uuid4()
            program = Program(
                id=program_id,
                name="Test Program",
                created_by_id=user_id,
                description="Test program"
            )
            db.add(program)
            db.commit()
        
        # Create multiple activities
        for i in range(3):
            log_activity(
                event_type=ActivityEventType.MODEL_TRAINED,
                target_type="model",
                target_id=uuid4(),
                target_name=f"Model {i}",
                actor_id=None,  # System action
                program_id=program_id,
            )
        
        # Verify we have unread notifications
        initial_unread = get_unread_count(user_id)
        assert initial_unread >= 3
        
        # Mark all as read
        marked_count = mark_all_notifications_read(user_id)
        assert marked_count >= 3
        
        # Verify all are read
        final_unread = get_unread_count(user_id)
        assert final_unread == initial_unread - marked_count

    @patch('amprenta_rag.services.activity.db_session')
    def test_log_activity_with_mocks(self, mock_db_session):
        """Test log_activity with mocked database to avoid schema issues."""
        # Setup mock
        mock_db = Mock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db_session.return_value.__exit__.return_value = None
        
        # Create mock event
        mock_event = Mock()
        mock_event.id = uuid4()
        mock_event.event_type = "compound_added"
        mock_event.target_type = "compound"
        mock_event.target_id = uuid4()
        
        # Mock the database operations
        mock_db.add.return_value = None
        mock_db.flush.return_value = None
        mock_db.commit.return_value = None
        
        # Mock route_notifications to return empty list
        with patch('amprenta_rag.services.activity.route_notifications') as mock_route:
            mock_route.return_value = []
            
            # Call the function with enum
            log_activity(
                event_type=ActivityEventType.COMPOUND_ADDED,
                target_type="compound",
                target_id=uuid4(),
                target_name="Test Compound"
            )
            
            # Verify the function was called without errors
            mock_db.add.assert_called_once()
            mock_db.flush.assert_called_once()
            mock_db.commit.assert_called_once()

    @patch('amprenta_rag.services.activity.db_session')
    def test_log_activity_with_string_event_type(self, mock_db_session):
        """Test log_activity with string event type to verify enum handling."""
        # Setup mock
        mock_db = Mock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db_session.return_value.__exit__.return_value = None
        
        # Mock the database operations
        mock_db.add.return_value = None
        mock_db.flush.return_value = None
        mock_db.commit.return_value = None
        
        # Mock route_notifications to return empty list
        with patch('amprenta_rag.services.activity.route_notifications') as mock_route:
            mock_route.return_value = []
            
            # Call the function with string (should not raise AttributeError)
            log_activity(
                event_type="experiment_created",  # String instead of enum
                target_type="experiment",
                target_id=uuid4(),
                target_name="Test Experiment"
            )
            
            # Verify the function was called without errors
            mock_db.add.assert_called_once()
            mock_db.flush.assert_called_once()
            mock_db.commit.assert_called_once()
