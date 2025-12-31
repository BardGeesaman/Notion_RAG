"""Collaboration session API endpoints for real-time collaborative editing."""

import os
from datetime import datetime, timedelta
from typing import List, Optional
from uuid import UUID, uuid4

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, ConfigDict, Field
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.database.models import User
from amprenta_rag.logging_utils import get_logger


# ============================================================================
# Collaboration Schemas
# ============================================================================


class CollaborationUser(BaseModel):
    """User information in collaboration context."""
    id: UUID
    username: str
    display_name: Optional[str] = None
    is_online: bool = True
    last_seen: Optional[datetime] = None


class CollaborationSession(BaseModel):
    """Collaboration session information."""
    document_id: str = Field(..., description="Unique document identifier")
    document_path: str = Field(..., description="Path to the collaborative document")
    document_type: str = Field(..., description="Type of document (notebook, text, etc.)")
    created_at: datetime = Field(..., description="Session creation timestamp")
    last_activity: datetime = Field(..., description="Last activity timestamp")
    participants: List[CollaborationUser] = Field(default_factory=list, description="Active participants")
    owner: CollaborationUser = Field(..., description="Session owner/creator")
    is_active: bool = Field(True, description="Whether session is currently active")
    document_title: Optional[str] = Field(None, description="Human-readable document title")


class SessionSummary(BaseModel):
    """Summary information for collaboration session listing."""
    document_id: str
    document_path: str
    document_type: str
    document_title: Optional[str] = None
    participant_count: int
    last_activity: datetime
    is_active: bool
    is_owner: bool = Field(..., description="Whether current user is the session owner")


class InviteRequest(BaseModel):
    """Request to invite a user to collaboration session."""
    username: str = Field(..., description="Username of user to invite")
    message: Optional[str] = Field(None, description="Optional invitation message")
    permissions: str = Field("read-write", description="Permissions for invited user")

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "username": "researcher123",
                "message": "Would you like to collaborate on this analysis?",
                "permissions": "read-write"
            }
        }
    )


class InviteResponse(BaseModel):
    """Response after sending collaboration invite."""
    success: bool
    invite_id: UUID = Field(..., description="Unique identifier for the invitation")
    invited_user: str = Field(..., description="Username of invited user")
    document_id: str = Field(..., description="Document being shared")
    message: str = Field(..., description="Status message")
    expires_at: datetime = Field(..., description="Invitation expiration time")


class SessionListResponse(BaseModel):
    """Response for listing collaboration sessions."""
    sessions: List[SessionSummary] = Field(default_factory=list)
    total_count: int = Field(..., description="Total number of sessions")
    active_count: int = Field(..., description="Number of currently active sessions")


class SessionDetailResponse(BaseModel):
    """Detailed response for specific collaboration session."""
    session: CollaborationSession
    user_permissions: str = Field(..., description="Current user's permissions in this session")
    connection_url: str = Field(..., description="WebSocket URL for joining the session")
    recent_activity: List[str] = Field(default_factory=list, description="Recent activity log")

logger = get_logger(__name__)

router = APIRouter(prefix="/collaboration", tags=["Collaboration"])


def _get_yjs_server_url() -> str:
    """Get Y.js WebSocket server URL from environment."""
    return os.environ.get("YJS_SERVER_URL", "ws://yjs-server:1234")


def _user_to_collaboration_user(user: User, is_online: bool = True) -> CollaborationUser:
    """Convert database User to CollaborationUser schema."""
    return CollaborationUser(
        id=user.id,
        username=user.username,
        display_name=getattr(user, 'display_name', None) or user.username,
        is_online=is_online,
        last_seen=datetime.utcnow() if is_online else None,
    )


def _mock_get_active_sessions(user_id: UUID) -> List[SessionSummary]:
    """Mock function to get active collaboration sessions for a user.
    
    In a real implementation, this would query the Y.js server or a
    collaboration service to get active sessions.
    """
    # Mock data for demonstration
    mock_sessions = [
        SessionSummary(
            document_id="notebook_analysis_001",
            document_path="/home/jovyan/work/data_analysis.ipynb",
            document_type="notebook",
            document_title="Customer Data Analysis",
            participant_count=3,
            last_activity=datetime.utcnow() - timedelta(minutes=5),
            is_active=True,
            is_owner=True,
        ),
        SessionSummary(
            document_id="script_ml_model_002",
            document_path="/home/jovyan/work/ml_model.py",
            document_type="python",
            document_title="ML Model Training Script",
            participant_count=1,
            last_activity=datetime.utcnow() - timedelta(hours=2),
            is_active=False,
            is_owner=False,
        ),
    ]
    
    # Filter sessions based on user (in real implementation)
    return mock_sessions


def _mock_get_session_details(document_id: str, user_id: UUID) -> CollaborationSession:
    """Mock function to get detailed session information.
    
    In a real implementation, this would query the Y.js server or
    collaboration service for session details.
    """
    if document_id == "nonexistent":
        return None
    
    # Mock session data
    mock_owner = CollaborationUser(
        id=user_id,
        username="researcher",
        display_name="Research User",
        is_online=True,
        last_seen=datetime.utcnow(),
    )
    
    mock_participants = [
        mock_owner,
        CollaborationUser(
            id=uuid4(),
            username="collaborator1",
            display_name="Collaborator One",
            is_online=True,
            last_seen=datetime.utcnow() - timedelta(minutes=2),
        ),
        CollaborationUser(
            id=uuid4(),
            username="collaborator2",
            display_name="Collaborator Two",
            is_online=False,
            last_seen=datetime.utcnow() - timedelta(hours=1),
        ),
    ]
    
    return CollaborationSession(
        document_id=document_id,
        document_path=f"/home/jovyan/work/{document_id}.ipynb",
        document_type="notebook",
        created_at=datetime.utcnow() - timedelta(hours=3),
        last_activity=datetime.utcnow() - timedelta(minutes=5),
        participants=mock_participants,
        owner=mock_owner,
        is_active=True,
        document_title="Collaborative Analysis Notebook",
    )


def _mock_send_invite(
    document_id: str, 
    invite: InviteRequest, 
    inviter_id: UUID
) -> InviteResponse:
    """Mock function to send collaboration invite.
    
    In a real implementation, this would:
    1. Validate the invited user exists
    2. Create an invitation record
    3. Send notification to the invited user
    4. Return invitation details
    """
    # Mock validation - check if user exists (simplified)
    if invite.username == "nonexistent_user":
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"User '{invite.username}' not found"
        )
    
    # Mock invite creation
    invite_id = uuid4()
    expires_at = datetime.utcnow() + timedelta(hours=24)
    
    logger.info(
        f"Collaboration invite sent: document={document_id}, "
        f"inviter={inviter_id}, invited={invite.username}"
    )
    
    return InviteResponse(
        success=True,
        invite_id=invite_id,
        invited_user=invite.username,
        document_id=document_id,
        message=f"Successfully invited {invite.username} to collaborate on document",
        expires_at=expires_at,
    )


@router.get("/sessions", response_model=SessionListResponse)
def list_sessions(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> SessionListResponse:
    """List active collaboration sessions for current user.
    
    Returns both sessions owned by the user and sessions they've been
    invited to participate in.
    """
    logger.info(f"Listing collaboration sessions for user {current_user.username}")
    
    try:
        # Get active sessions for the user
        sessions = _mock_get_active_sessions(current_user.id)
        
        # Calculate counts
        total_count = len(sessions)
        active_count = sum(1 for session in sessions if session.is_active)
        
        logger.info(
            f"Found {total_count} sessions for user {current_user.username} "
            f"({active_count} active)"
        )
        
        return SessionListResponse(
            sessions=sessions,
            total_count=total_count,
            active_count=active_count,
        )
        
    except Exception as e:
        logger.error(f"Failed to list sessions for user {current_user.username}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to retrieve collaboration sessions"
        )


@router.get("/sessions/{document_id}", response_model=SessionDetailResponse)
def get_session(
    document_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> SessionDetailResponse:
    """Get collaboration session details for a specific document.
    
    Returns detailed information about the session including participants,
    activity, and connection information.
    """
    logger.info(
        f"Getting session details for document {document_id} "
        f"by user {current_user.username}"
    )
    
    try:
        # Get session details
        session = _mock_get_session_details(document_id, current_user.id)
        
        if not session:
            logger.warning(
                f"Session not found: document={document_id}, "
                f"user={current_user.username}"
            )
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Collaboration session for document '{document_id}' not found"
            )
        
        # Check if user has access to this session
        user_in_session = any(
            participant.id == current_user.id 
            for participant in session.participants
        )
        
        if not user_in_session:
            logger.warning(
                f"Access denied: user {current_user.username} "
                f"not participant in session {document_id}"
            )
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="You don't have access to this collaboration session"
            )
        
        # Determine user permissions (simplified)
        is_owner = session.owner.id == current_user.id
        user_permissions = "read-write" if is_owner else "read-write"
        
        # Build WebSocket connection URL
        yjs_url = _get_yjs_server_url()
        connection_url = f"{yjs_url}/doc/{document_id}"
        
        # Mock recent activity
        recent_activity = [
            f"{session.owner.username} created the session",
            "collaborator1 joined the session",
            "collaborator2 edited cell 3",
            "collaborator1 added a comment",
        ]
        
        logger.info(
            f"Session details retrieved for document {document_id} "
            f"with {len(session.participants)} participants"
        )
        
        return SessionDetailResponse(
            session=session,
            user_permissions=user_permissions,
            connection_url=connection_url,
            recent_activity=recent_activity,
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(
            f"Failed to get session details for document {document_id}: {e}"
        )
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to retrieve session details"
        )


@router.post("/sessions/{document_id}/invite", response_model=InviteResponse)
def invite_to_session(
    document_id: str,
    invite: InviteRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> InviteResponse:
    """Invite a user to join a collaboration session.
    
    Sends an invitation to the specified user to join the collaborative
    editing session for the given document.
    """
    logger.info(
        f"Inviting user {invite.username} to session {document_id} "
        f"by {current_user.username}"
    )
    
    try:
        # Check if session exists and user has permission to invite
        session = _mock_get_session_details(document_id, current_user.id)
        
        if not session:
            logger.warning(f"Session not found for invite: document={document_id}")
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Collaboration session for document '{document_id}' not found"
            )
        
        # Check if user has permission to invite (owner or admin)
        is_owner = session.owner.id == current_user.id
        if not is_owner:
            logger.warning(
                f"Invite permission denied: user {current_user.username} "
                f"not owner of session {document_id}"
            )
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Only session owners can invite users"
            )
        
        # Check if user is trying to invite themselves
        if invite.username == current_user.username:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Cannot invite yourself to a session"
            )
        
        # Check if user is already in the session
        already_invited = any(
            participant.username == invite.username 
            for participant in session.participants
        )
        
        if already_invited:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"User '{invite.username}' is already in this session"
            )
        
        # Send the invite
        invite_response = _mock_send_invite(document_id, invite, current_user.id)
        
        logger.info(
            f"Collaboration invite sent successfully: "
            f"document={document_id}, invited={invite.username}, "
            f"invite_id={invite_response.invite_id}"
        )
        
        return invite_response
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(
            f"Failed to send invite for document {document_id} "
            f"to user {invite.username}: {e}"
        )
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to send collaboration invite"
        )
