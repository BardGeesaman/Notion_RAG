"""Permission checking utilities for teams and projects."""
from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from amprenta_rag.database.models import Team, TeamMember, Project
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_user_teams(user_id: str, db) -> List[Team]:
    """
    Get all teams a user is a member of.
    
    Args:
        user_id: UUID of the user
        db: Database session
        
    Returns:
        List of Team objects
    """
    try:
        memberships = db.query(TeamMember).filter(
            TeamMember.user_id == UUID(user_id)
        ).all()
        teams = [m.team for m in memberships]
        logger.debug("[PERMISSIONS] User %s is in %d teams", user_id[:8], len(teams))
        return teams
    except Exception as e:
        logger.error("[PERMISSIONS] Error getting user teams: %r", e)
        return []


def get_user_projects(user_id: str, db) -> List[Project]:
    """
    Get all projects a user can access (via team membership + public projects).
    
    Args:
        user_id: UUID of the user
        db: Database session
        
    Returns:
        List of Project objects
    """
    try:
        # Get projects via team membership
        memberships = db.query(TeamMember).filter(
            TeamMember.user_id == UUID(user_id)
        ).all()
        team_ids = [m.team_id for m in memberships]
        
        team_projects = []
        if team_ids:
            team_projects = db.query(Project).filter(
                Project.team_id.in_(team_ids)
            ).all()
        
        # Get public projects
        public_projects = db.query(Project).filter(
            Project.is_public == True
        ).all()
        
        # Combine and deduplicate
        all_projects = {p.id: p for p in team_projects + public_projects}.values()
        projects = list(all_projects)
        
        logger.debug("[PERMISSIONS] User %s can access %d projects", user_id[:8], len(projects))
        return projects
    except Exception as e:
        logger.error("[PERMISSIONS] Error getting user projects: %r", e)
        return []


def can_view_project(user_id: str, project_id: str, db) -> bool:
    """
    Check if user can view a project.
    
    Args:
        user_id: UUID of the user
        project_id: UUID of the project
        db: Database session
        
    Returns:
        True if user can view the project
    """
    try:
        project = db.query(Project).filter(Project.id == UUID(project_id)).first()
        if not project:
            return False
        
        # Public projects are viewable by everyone
        if project.is_public:
            return True
        
        # Check if user is a member of the project's team
        membership = db.query(TeamMember).filter(
            TeamMember.user_id == UUID(user_id),
            TeamMember.team_id == project.team_id
        ).first()
        
        can_view = membership is not None
        logger.debug("[PERMISSIONS] User %s can_view_project %s: %s", user_id[:8], project_id[:8], can_view)
        return can_view
    except Exception as e:
        logger.error("[PERMISSIONS] Error checking view permission: %r", e)
        return False


def can_edit_project(user_id: str, project_id: str, db) -> bool:
    """
    Check if user can edit a project (owner/admin only).
    
    Args:
        user_id: UUID of the user
        project_id: UUID of the project
        db: Database session
        
    Returns:
        True if user can edit the project
    """
    try:
        project = db.query(Project).filter(Project.id == UUID(project_id)).first()
        if not project:
            return False
        
        # Check if user is owner or admin of the project's team
        membership = db.query(TeamMember).filter(
            TeamMember.user_id == UUID(user_id),
            TeamMember.team_id == project.team_id,
            TeamMember.role.in_(["owner", "admin"])
        ).first()
        
        can_edit = membership is not None
        logger.debug("[PERMISSIONS] User %s can_edit_project %s: %s", user_id[:8], project_id[:8], can_edit)
        return can_edit
    except Exception as e:
        logger.error("[PERMISSIONS] Error checking edit permission: %r", e)
        return False


def get_team_role(user_id: str, team_id: str, db) -> Optional[str]:
    """
    Get user's role in a team.
    
    Args:
        user_id: UUID of the user
        team_id: UUID of the team
        db: Database session
        
    Returns:
        Role string (owner, admin, member, viewer) or None if not a member
    """
    try:
        membership = db.query(TeamMember).filter(
            TeamMember.user_id == UUID(user_id),
            TeamMember.team_id == UUID(team_id)
        ).first()
        
        role = membership.role if membership else None
        logger.debug("[PERMISSIONS] User %s role in team %s: %s", user_id[:8], team_id[:8], role)
        return role
    except Exception as e:
        logger.error("[PERMISSIONS] Error getting team role: %r", e)
        return None
