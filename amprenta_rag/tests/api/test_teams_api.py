"""Tests for teams API endpoints."""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import MagicMock
from uuid import uuid4
from datetime import datetime, timezone

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.models.auth import User, Team, TeamMember


# Test constants
TEST_USER_ID = uuid4()


def mock_current_user():
    """Mock current user for authentication."""
    return User(
        id=TEST_USER_ID,
        username=f"test_user_{uuid4().hex[:8]}",
        email=f"test_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )


def mock_db():
    """Mock database session."""
    db = MagicMock()
    
    # Mock team memberships query
    mock_memberships = [
        TeamMember(
            id=uuid4(),
            team_id=uuid4(),
            user_id=TEST_USER_ID,
            role="owner",
            joined_at=datetime.now(timezone.utc)
        ),
        TeamMember(
            id=uuid4(),
            team_id=uuid4(),
            user_id=TEST_USER_ID,
            role="member",
            joined_at=datetime.now(timezone.utc)
        )
    ]
    
    # Mock teams
    mock_teams = [
        Team(
            id=mock_memberships[0].team_id,
            name="Team Alpha",
            description="First team",
            created_at=datetime.now(timezone.utc)
        ),
        Team(
            id=mock_memberships[1].team_id,
            name="Team Beta", 
            description="Second team",
            created_at=datetime.now(timezone.utc)
        )
    ]
    
    # Configure mock query chains
    def mock_query_side_effect(model):
        if model == TeamMember:
            query_mock = MagicMock()
            query_mock.filter.return_value.all.return_value = mock_memberships
            return query_mock
        elif model == Team:
            query_mock = MagicMock()
            # Return teams based on which one is being queried
            def filter_side_effect(*args):
                filter_mock = MagicMock()
                filter_mock.first.return_value = mock_teams[0]  # Default to first team
                return filter_mock
            query_mock.filter.side_effect = filter_side_effect
            return query_mock
        else:
            return MagicMock()
    
    db.query.side_effect = mock_query_side_effect
    
    # Mock count for member count
    def mock_count_side_effect():
        return 3  # Mock member count
    
    # Setup the count mock for the specific query chain
    count_mock = MagicMock()
    count_mock.count.return_value = 3
    db.query.return_value.filter.return_value = count_mock
    
    return db


@pytest.fixture
def client():
    """Test client with mocked dependencies."""
    app.dependency_overrides[get_current_user] = mock_current_user
    app.dependency_overrides[get_db] = mock_db
    
    try:
        with TestClient(app) as test_client:
            yield test_client
    finally:
        app.dependency_overrides.clear()


def test_my_teams_endpoint(client):
    """Test getting user's teams via API."""
    # Act
    response = client.get("/api/v1/teams/my-teams")
    
    # Assert
    assert response.status_code == 200
    data = response.json()
    assert "teams" in data
    # The mock setup should return teams data
    teams = data["teams"]
    assert isinstance(teams, list)
    
    # If teams are returned, verify structure
    if teams:
        team = teams[0]
        assert "id" in team
        assert "name" in team
        assert "member_count" in team
        assert "user_role" in team


def test_create_team_endpoint(client):
    """Test creating a team via API."""
    # Arrange
    team_data = {
        "name": f"Test Team {uuid4().hex[:8]}",
        "description": "A test team for collaboration"
    }
    
    # Mock successful team creation
    mock_team = Team(
        id=uuid4(),
        name=team_data["name"],
        description=team_data["description"],
        created_at=datetime.now(timezone.utc)
    )
    
    # Configure the mock database for team creation
    def mock_db_for_creation():
        db = MagicMock()
        
        # Mock no existing team with same name
        db.query.return_value.filter.return_value.first.return_value = None
        
        # Mock the team object after add/flush
        def mock_add(obj):
            if hasattr(obj, 'name') and obj.name == team_data["name"]:
                obj.id = mock_team.id
                obj.created_at = mock_team.created_at
        
        db.add.side_effect = mock_add
        db.flush.return_value = None
        db.commit.return_value = None
        
        return db
    
    # Override the mock for this specific test
    app.dependency_overrides[get_db] = mock_db_for_creation
    
    try:
        # Act
        response = client.post("/api/v1/teams/teams", json=team_data)
        
        # Assert
        assert response.status_code == 201
        data = response.json()
        assert data["name"] == team_data["name"]
        assert data["description"] == team_data["description"]
        assert data["member_count"] == 1
        assert data["user_role"] == "owner"
        
    finally:
        # Restore original mock
        app.dependency_overrides[get_db] = mock_db


def test_create_team_endpoint_duplicate_name(client):
    """Test creating team with duplicate name."""
    # Arrange
    team_data = {
        "name": "Existing Team",
        "description": "This should fail"
    }
    
    # Mock existing team with same name
    existing_team = Team(
        id=uuid4(),
        name=team_data["name"],
        description="Already exists",
        created_at=datetime.now(timezone.utc)
    )
    
    def mock_db_with_existing():
        db = MagicMock()
        db.query.return_value.filter.return_value.first.return_value = existing_team
        return db
    
    app.dependency_overrides[get_db] = mock_db_with_existing
    
    try:
        # Act
        response = client.post("/api/v1/teams/teams", json=team_data)
        
        # Assert
        assert response.status_code == 400
        assert "Team name already exists" in response.json()["detail"]
        
    finally:
        app.dependency_overrides[get_db] = mock_db


def test_get_team_members_endpoint(client):
    """Test getting team members via API."""
    # Arrange
    team_id = uuid4()
    
    mock_team = Team(
        id=team_id,
        name="Test Team",
        description="Test team",
        created_at=datetime.now(timezone.utc)
    )
    
    mock_user_membership = TeamMember(
        id=uuid4(),
        team_id=team_id,
        user_id=TEST_USER_ID,
        role="owner",
        joined_at=datetime.now(timezone.utc)
    )
    
    mock_other_membership = TeamMember(
        id=uuid4(),
        team_id=team_id,
        user_id=uuid4(),
        role="member",
        joined_at=datetime.now(timezone.utc)
    )
    
    mock_other_user = User(
        id=mock_other_membership.user_id,
        username=f"other_user_{uuid4().hex[:8]}",
        email=f"other_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )
    
    def mock_db_for_members():
        db = MagicMock()
        
        # Mock team existence check
        def team_query_side_effect(*args):
            query_mock = MagicMock()
            query_mock.filter.return_value.first.return_value = mock_team
            return query_mock
        
        # Mock user membership check
        def membership_query_side_effect(*args):
            query_mock = MagicMock()
            if "team_id" in str(args):
                # For membership queries
                if "user_id" in str(args):
                    query_mock.filter.return_value.first.return_value = mock_user_membership
                else:
                    query_mock.filter.return_value.all.return_value = [mock_user_membership, mock_other_membership]
            return query_mock
        
        # Mock user details lookup
        def user_query_side_effect(*args):
            query_mock = MagicMock()
            query_mock.filter.return_value.first.return_value = mock_other_user
            return query_mock
        
        # Configure query routing
        def query_router(model):
            if model == Team:
                return team_query_side_effect(model)
            elif model == TeamMember:
                return membership_query_side_effect(model)
            elif model == User:
                return user_query_side_effect(model)
            return MagicMock()
        
        db.query.side_effect = query_router
        return db
    
    app.dependency_overrides[get_db] = mock_db_for_members
    
    try:
        # Act
        response = client.get(f"/api/v1/teams/teams/{team_id}/members")
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        # Should have at least the mock member
        
    finally:
        app.dependency_overrides[get_db] = mock_db


def test_get_team_members_endpoint_not_member(client):
    """Test getting team members when user is not a member."""
    # Arrange
    team_id = uuid4()
    
    mock_team = Team(
        id=team_id,
        name="Private Team",
        description="User is not a member",
        created_at=datetime.now(timezone.utc)
    )
    
    def mock_db_no_membership():
        db = MagicMock()
        
        def query_side_effect(model):
            query_mock = MagicMock()
            if model == Team:
                query_mock.filter.return_value.first.return_value = mock_team
            elif model == TeamMember:
                query_mock.filter.return_value.first.return_value = None  # No membership
            return query_mock
        
        db.query.side_effect = query_side_effect
        return db
    
    app.dependency_overrides[get_db] = mock_db_no_membership
    
    try:
        # Act
        response = client.get(f"/api/v1/teams/teams/{team_id}/members")
        
        # Assert
        assert response.status_code == 403
        assert "Access denied" in response.json()["detail"]
        
    finally:
        app.dependency_overrides[get_db] = mock_db


def test_get_team_endpoint(client):
    """Test getting team details via API."""
    # Arrange
    team_id = uuid4()
    
    mock_team = Team(
        id=team_id,
        name="Test Team Details",
        description="Team for testing details endpoint",
        created_at=datetime.now(timezone.utc)
    )
    
    mock_membership = TeamMember(
        id=uuid4(),
        team_id=team_id,
        user_id=TEST_USER_ID,
        role="admin",
        joined_at=datetime.now(timezone.utc)
    )
    
    def mock_db_for_team_details():
        db = MagicMock()
        
        def query_side_effect(model):
            query_mock = MagicMock()
            if model == Team:
                query_mock.filter.return_value.first.return_value = mock_team
            elif model == TeamMember:
                # For membership check
                if hasattr(query_mock.filter.return_value, 'first'):
                    query_mock.filter.return_value.first.return_value = mock_membership
                # For member count
                query_mock.filter.return_value.count.return_value = 5
            return query_mock
        
        db.query.side_effect = query_side_effect
        return db
    
    app.dependency_overrides[get_db] = mock_db_for_team_details
    
    try:
        # Act
        response = client.get(f"/api/v1/teams/teams/{team_id}")
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(team_id)
        assert data["name"] == mock_team.name
        assert data["description"] == mock_team.description
        assert data["user_role"] == "admin"
        assert data["member_count"] == 5
        
    finally:
        app.dependency_overrides[get_db] = mock_db
