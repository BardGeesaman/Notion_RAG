"""Entity sharing service for collaboration features."""

import logging
from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session
from sqlalchemy import and_, or_

from amprenta_rag.models.auth import EntityShare, TeamMember
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import ActivityEventType
from amprenta_rag.services.activity import log_activity

logger = logging.getLogger(__name__)


def share_entity(
    entity_type: str,
    entity_id: UUID,
    user_id: Optional[UUID] = None,
    team_id: Optional[UUID] = None,
    permission: str = "view",
    shared_by: UUID = None,
    db: Optional[Session] = None
) -> Optional[EntityShare]:
    """
    Share an entity with a user or team.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity to share
        user_id: UUID of user to share with (exclusive with team_id)
        team_id: UUID of team to share with (exclusive with user_id)
        permission: Permission level (view, edit, admin)
        shared_by: UUID of user creating the share
        db: Database session (optional)
    
    Returns:
        EntityShare object if successful, None if failed
    """
    if not shared_by:
        logger.error("shared_by is required")
        return None
        
    if not user_id and not team_id:
        logger.error("Either user_id or team_id must be provided")
        return None
        
    if user_id and team_id:
        logger.error("Cannot share with both user and team")
        return None
        
    if permission not in ["view", "edit", "admin"]:
        logger.error(f"Invalid permission: {permission}")
        return None
    
    try:
        session_provided = db is not None
        if not session_provided:
            db = db_session()
            
        with db if session_provided else db_session() as session:
            # Check if share already exists
            existing_share = session.query(EntityShare).filter(
                and_(
                    EntityShare.entity_type == entity_type,
                    EntityShare.entity_id == entity_id,
                    or_(
                        EntityShare.shared_with_user_id == user_id,
                        EntityShare.shared_with_team_id == team_id
                    )
                )
            ).first()
            
            if existing_share:
                # Update existing share
                existing_share.permission = permission
                existing_share.shared_by_id = shared_by
                session.commit()
                session.expunge(existing_share)
                
                # Log activity event for share update
                try:
                    log_activity(
                        event_type=ActivityEventType.ENTITY_SHARED,
                        target_type=entity_type,
                        target_id=entity_id,
                        target_name=f"{entity_type}:{entity_id}",
                        actor_id=shared_by,
                        program_id=None,  # Entity sharing is not program-specific
                        metadata={
                            "entity_type": entity_type,
                            "entity_id": str(entity_id),
                            "shared_with_user_id": str(user_id) if user_id else None,
                            "shared_with_team_id": str(team_id) if team_id else None,
                            "permission": permission,
                            "action": "updated"
                        }
                    )
                except Exception as e:
                    logger.error(f"Failed to log ENTITY_SHARED activity (update): {e}")
                
                return existing_share
            
            # Create new share
            share = EntityShare(
                entity_type=entity_type,
                entity_id=entity_id,
                shared_with_user_id=user_id,
                shared_with_team_id=team_id,
                permission=permission,
                shared_by_id=shared_by
            )
            
            session.add(share)
            session.commit()
            session.expunge(share)
            
            logger.info(
                f"Entity shared: {entity_type}:{entity_id} with "
                f"{'user:' + str(user_id) if user_id else 'team:' + str(team_id)} "
                f"permission:{permission}"
            )
            
            # Log activity event
            try:
                log_activity(
                    event_type=ActivityEventType.ENTITY_SHARED,
                    target_type=entity_type,
                    target_id=entity_id,
                    target_name=f"{entity_type}:{entity_id}",
                    actor_id=shared_by,
                    program_id=None,  # Entity sharing is not program-specific
                    metadata={
                        "entity_type": entity_type,
                        "entity_id": str(entity_id),
                        "shared_with_user_id": str(user_id) if user_id else None,
                        "shared_with_team_id": str(team_id) if team_id else None,
                        "permission": permission
                    }
                )
            except Exception as e:
                logger.error(f"Failed to log ENTITY_SHARED activity: {e}")
            
            return share
            
    except Exception as e:
        logger.error(f"Failed to share entity {entity_type}:{entity_id}: {e}")
        return None


def unshare_entity(
    entity_type: str,
    entity_id: UUID,
    user_id: Optional[UUID] = None,
    team_id: Optional[UUID] = None,
    db: Optional[Session] = None
) -> bool:
    """
    Remove sharing for an entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        user_id: UUID of user to unshare from (exclusive with team_id)
        team_id: UUID of team to unshare from (exclusive with user_id)
        db: Database session (optional)
    
    Returns:
        True if successful, False if failed
    """
    if not user_id and not team_id:
        logger.error("Either user_id or team_id must be provided")
        return False
        
    if user_id and team_id:
        logger.error("Cannot specify both user_id and team_id")
        return False
    
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            share = session.query(EntityShare).filter(
                and_(
                    EntityShare.entity_type == entity_type,
                    EntityShare.entity_id == entity_id,
                    or_(
                        EntityShare.shared_with_user_id == user_id,
                        EntityShare.shared_with_team_id == team_id
                    )
                )
            ).first()
            
            if not share:
                logger.warning(f"No share found for {entity_type}:{entity_id}")
                return False
            
            # Store share info before deletion for activity logging
            share_info = {
                "entity_type": share.entity_type,
                "entity_id": str(share.entity_id),
                "shared_with_user_id": str(share.shared_with_user_id) if share.shared_with_user_id else None,
                "shared_with_team_id": str(share.shared_with_team_id) if share.shared_with_team_id else None,
                "permission": share.permission,
                "shared_by_id": share.shared_by_id
            }
            
            session.delete(share)
            session.commit()
            
            logger.info(
                f"Entity unshared: {entity_type}:{entity_id} from "
                f"{'user:' + str(user_id) if user_id else 'team:' + str(team_id)}"
            )
            
            # Log activity event
            try:
                log_activity(
                    event_type=ActivityEventType.ENTITY_UNSHARED,
                    target_type=entity_type,
                    target_id=entity_id,
                    target_name=f"{entity_type}:{entity_id}",
                    actor_id=share_info["shared_by_id"],  # Use original sharer as actor
                    program_id=None,  # Entity sharing is not program-specific
                    metadata=share_info
                )
            except Exception as e:
                logger.error(f"Failed to log ENTITY_UNSHARED activity: {e}")
            
            return True
            
    except Exception as e:
        logger.error(f"Failed to unshare entity {entity_type}:{entity_id}: {e}")
        return False


def list_shares_for_entity(
    entity_type: str,
    entity_id: UUID,
    db: Optional[Session] = None
) -> List[EntityShare]:
    """
    List all shares for a specific entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        db: Database session (optional)
    
    Returns:
        List of EntityShare objects
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            shares = session.query(EntityShare).filter(
                and_(
                    EntityShare.entity_type == entity_type,
                    EntityShare.entity_id == entity_id
                )
            ).all()
            
            # Detach from session
            for share in shares:
                session.expunge(share)
            
            return shares
            
    except Exception as e:
        logger.error(f"Failed to list shares for {entity_type}:{entity_id}: {e}")
        return []


def get_my_shares(
    user_id: UUID,
    entity_type: Optional[str] = None,
    db: Optional[Session] = None
) -> List[EntityShare]:
    """
    Get all entities shared with a user (directly or through teams).
    
    Args:
        user_id: UUID of the user
        entity_type: Optional filter by entity type
        db: Database session (optional)
    
    Returns:
        List of EntityShare objects
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            # Get user's team memberships
            team_ids = session.query(TeamMember.team_id).filter(
                TeamMember.user_id == user_id
            ).subquery()
            
            # Build query for shares
            query = session.query(EntityShare).filter(
                or_(
                    EntityShare.shared_with_user_id == user_id,
                    EntityShare.shared_with_team_id.in_(team_ids)
                )
            )
            
            if entity_type:
                query = query.filter(EntityShare.entity_type == entity_type)
            
            shares = query.all()
            
            # Detach from session
            for share in shares:
                session.expunge(share)
            
            return shares
            
    except Exception as e:
        logger.error(f"Failed to get shares for user {user_id}: {e}")
        return []


def check_share_permission(
    user_id: UUID,
    entity_type: str,
    entity_id: UUID,
    required_permission: str = "view",
    db: Optional[Session] = None
) -> bool:
    """
    Check if user has permission to access an entity through sharing.
    
    Args:
        user_id: UUID of the user
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        required_permission: Required permission level (view, edit, admin)
        db: Database session (optional)
    
    Returns:
        True if user has permission, False otherwise
    """
    # Permission hierarchy: admin > edit > view
    permission_levels = {"view": 1, "edit": 2, "admin": 3}
    
    if required_permission not in permission_levels:
        logger.error(f"Invalid required_permission: {required_permission}")
        return False
    
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            # Get user's team memberships
            team_ids = session.query(TeamMember.team_id).filter(
                TeamMember.user_id == user_id
            ).subquery()
            
            # Check for shares (direct user or team)
            shares = session.query(EntityShare).filter(
                and_(
                    EntityShare.entity_type == entity_type,
                    EntityShare.entity_id == entity_id,
                    or_(
                        EntityShare.shared_with_user_id == user_id,
                        EntityShare.shared_with_team_id.in_(team_ids)
                    )
                )
            ).all()
            
            # Check if any share meets the permission requirement
            required_level = permission_levels[required_permission]
            
            for share in shares:
                share_level = permission_levels.get(share.permission, 0)
                if share_level >= required_level:
                    return True
            
            return False
            
    except Exception as e:
        logger.error(
            f"Failed to check permission for user {user_id} "
            f"on {entity_type}:{entity_id}: {e}"
        )
        return False
