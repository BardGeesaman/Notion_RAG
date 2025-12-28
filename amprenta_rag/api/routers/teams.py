"""API router for team management functionality."""

from typing import List
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.api.schemas import (
    TeamCreate,
    TeamResponse,
    TeamMemberResponse,
    TeamListResponse,
)
from amprenta_rag.models.auth import User, Team, TeamMember

router = APIRouter(prefix="/teams", tags=["teams"])


@router.get("/my-teams")
def get_my_teams(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> TeamListResponse:
    """
    Get all teams the current user is a member of.
    
    Args:
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        List of user's teams with membership info
    """
    try:
        # Get teams where user is a member
        team_memberships = (
            db.query(TeamMember)
            .filter(TeamMember.user_id == current_user.id)
            .all()
        )
        
        teams = []
        for membership in team_memberships:
            team = db.query(Team).filter(Team.id == membership.team_id).first()
            if team:
                # Count members
                member_count = (
                    db.query(TeamMember)
                    .filter(TeamMember.team_id == team.id)
                    .count()
                )
                
                teams.append(TeamResponse(
                    id=team.id,
                    name=team.name,
                    description=team.description,
                    created_at=team.created_at,
                    member_count=member_count,
                    user_role=membership.role,
                ))
        
        return TeamListResponse(teams=teams)
        
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to fetch teams: {str(e)}"
        )


@router.post("/teams", status_code=status.HTTP_201_CREATED)
def create_team(
    team_data: TeamCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> TeamResponse:
    """
    Create a new team.
    
    Args:
        team_data: Team creation data
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        Created team details
    """
    try:
        # Check if team name already exists
        existing_team = db.query(Team).filter(Team.name == team_data.name).first()
        if existing_team:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Team name already exists"
            )
        
        # Create team
        team = Team(
            name=team_data.name,
            description=team_data.description,
        )
        
        db.add(team)
        db.flush()  # Get the team ID
        
        # Add creator as team owner
        membership = TeamMember(
            team_id=team.id,
            user_id=current_user.id,
            role="owner"
        )
        
        db.add(membership)
        db.commit()
        
        return TeamResponse(
            id=team.id,
            name=team.name,
            description=team.description,
            created_at=team.created_at,
            member_count=1,
            user_role="owner",
        )
        
    except HTTPException:
        raise
    except Exception as e:
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to create team: {str(e)}"
        )


@router.get("/teams/{team_id}/members")
def get_team_members(
    team_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> List[TeamMemberResponse]:
    """
    Get all members of a team.
    
    Args:
        team_id: UUID of the team
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        List of team members
    """
    try:
        # Check if team exists
        team = db.query(Team).filter(Team.id == team_id).first()
        if not team:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Team not found"
            )
        
        # Check if user is a member of the team
        user_membership = (
            db.query(TeamMember)
            .filter(
                TeamMember.team_id == team_id,
                TeamMember.user_id == current_user.id
            )
            .first()
        )
        
        if not user_membership:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Access denied. You are not a member of this team."
            )
        
        # Get all team members
        memberships = (
            db.query(TeamMember)
            .filter(TeamMember.team_id == team_id)
            .all()
        )
        
        members = []
        for membership in memberships:
            user = db.query(User).filter(User.id == membership.user_id).first()
            if user:
                members.append(TeamMemberResponse(
                    id=membership.id,
                    user_id=user.id,
                    username=user.username,
                    email=user.email,
                    role=membership.role,
                    joined_at=membership.joined_at,
                ))
        
        return members
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to fetch team members: {str(e)}"
        )


@router.get("/teams/{team_id}")
def get_team(
    team_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> TeamResponse:
    """
    Get team details.
    
    Args:
        team_id: UUID of the team
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        Team details
    """
    try:
        # Check if team exists
        team = db.query(Team).filter(Team.id == team_id).first()
        if not team:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Team not found"
            )
        
        # Check if user is a member of the team
        user_membership = (
            db.query(TeamMember)
            .filter(
                TeamMember.team_id == team_id,
                TeamMember.user_id == current_user.id
            )
            .first()
        )
        
        if not user_membership:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Access denied. You are not a member of this team."
            )
        
        # Count members
        member_count = (
            db.query(TeamMember)
            .filter(TeamMember.team_id == team_id)
            .count()
        )
        
        return TeamResponse(
            id=team.id,
            name=team.name,
            description=team.description,
            created_at=team.created_at,
            member_count=member_count,
            user_role=user_membership.role,
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to fetch team: {str(e)}"
        )
