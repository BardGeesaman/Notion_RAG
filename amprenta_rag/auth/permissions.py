"""Entity-level permissions and authorization checks."""

import logging
from typing import Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.session import db_session
from amprenta_rag.models.auth import User
from amprenta_rag.services.sharing import check_share_permission

logger = logging.getLogger(__name__)


def can_view_entity(
    user_id: UUID,
    entity_type: str,
    entity_id: UUID,
    db: Optional[Session] = None
) -> bool:
    """
    Check if user can view an entity.
    
    This checks:
    1. Entity ownership (if user created/owns the entity)
    2. Sharing permissions (if entity is shared with user/team)
    3. Admin privileges (if user is admin)
    
    Args:
        user_id: UUID of the user
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        db: Database session (optional)
    
    Returns:
        True if user can view entity, False otherwise
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            # Check if user is admin
            user = session.query(User).filter(User.id == user_id).first()
            if user and user.role == "admin":
                return True
            
            # Check entity ownership
            if _is_entity_owner(user_id, entity_type, entity_id, session):
                return True
            
            # Check sharing permissions
            if check_share_permission(
                user_id=user_id,
                entity_type=entity_type,
                entity_id=entity_id,
                required_permission="view",
                db=session
            ):
                return True
            
            return False
            
    except Exception as e:
        logger.error(f"Failed to check view permission for {entity_type}:{entity_id}: {e}")
        return False


def can_edit_entity(
    user_id: UUID,
    entity_type: str,
    entity_id: UUID,
    db: Optional[Session] = None
) -> bool:
    """
    Check if user can edit an entity.
    
    This checks:
    1. Entity ownership (if user created/owns the entity)
    2. Sharing permissions with edit or admin level
    3. Admin privileges (if user is admin)
    
    Args:
        user_id: UUID of the user
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        db: Database session (optional)
    
    Returns:
        True if user can edit entity, False otherwise
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            # Check if user is admin
            user = session.query(User).filter(User.id == user_id).first()
            if user and user.role == "admin":
                return True
            
            # Check entity ownership
            if _is_entity_owner(user_id, entity_type, entity_id, session):
                return True
            
            # Check sharing permissions (edit or admin)
            if check_share_permission(
                user_id=user_id,
                entity_type=entity_type,
                entity_id=entity_id,
                required_permission="edit",
                db=session
            ):
                return True
            
            return False
            
    except Exception as e:
        logger.error(f"Failed to check edit permission for {entity_type}:{entity_id}: {e}")
        return False


def can_share_entity(
    user_id: UUID,
    entity_type: str,
    entity_id: UUID,
    db: Optional[Session] = None
) -> bool:
    """
    Check if user can share an entity.
    
    This checks:
    1. Entity ownership (if user created/owns the entity)
    2. Sharing permissions with admin level
    3. Admin privileges (if user is admin)
    
    Args:
        user_id: UUID of the user
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        db: Database session (optional)
    
    Returns:
        True if user can share entity, False otherwise
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            # Check if user is admin
            user = session.query(User).filter(User.id == user_id).first()
            if user and user.role == "admin":
                return True
            
            # Check entity ownership
            if _is_entity_owner(user_id, entity_type, entity_id, session):
                return True
            
            # Check sharing permissions (admin level required)
            if check_share_permission(
                user_id=user_id,
                entity_type=entity_type,
                entity_id=entity_id,
                required_permission="admin",
                db=session
            ):
                return True
            
            return False
            
    except Exception as e:
        logger.error(f"Failed to check share permission for {entity_type}:{entity_id}: {e}")
        return False


def _is_entity_owner(
    user_id: UUID,
    entity_type: str,
    entity_id: UUID,
    db: Session
) -> bool:
    """
    Check if user owns/created an entity.
    
    This is a helper function that checks entity ownership across different tables.
    
    Args:
        user_id: UUID of the user
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        db: Database session
    
    Returns:
        True if user owns the entity, False otherwise
    """
    try:
        # Import here to avoid circular imports
        from amprenta_rag.database.models import Dataset, Experiment, Compound, Signature
        
        if entity_type == "dataset":
            entity = db.query(Dataset).filter(Dataset.id == entity_id).first()
            return entity and entity.created_by_id == user_id
            
        elif entity_type == "experiment":
            entity = db.query(Experiment).filter(Experiment.id == entity_id).first()
            return entity and entity.created_by_id == user_id
            
        elif entity_type == "compound":
            entity = db.query(Compound).filter(Compound.id == entity_id).first()
            # Compounds might not have created_by_id, so check if user has edit access
            # For now, assume compounds are publicly viewable but editable by admins only
            return False
            
        elif entity_type == "signature":
            entity = db.query(Signature).filter(Signature.id == entity_id).first()
            return entity and entity.created_by_id == user_id
            
        else:
            logger.warning(f"Unknown entity_type: {entity_type}")
            return False
            
    except Exception as e:
        logger.error(f"Failed to check ownership for {entity_type}:{entity_id}: {e}")
        return False


def get_user_role(user_id: UUID, db: Optional[Session] = None) -> Optional[str]:
    """
    Get user's role.
    
    Args:
        user_id: UUID of the user
        db: Database session (optional)
    
    Returns:
        User role string or None if user not found
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            user = session.query(User).filter(User.id == user_id).first()
            return user.role if user else None
            
    except Exception as e:
        logger.error(f"Failed to get role for user {user_id}: {e}")
        return None


def is_admin(user_id: UUID, db: Optional[Session] = None) -> bool:
    """
    Check if user is an admin.
    
    Args:
        user_id: UUID of the user
        db: Database session (optional)
    
    Returns:
        True if user is admin, False otherwise
    """
    role = get_user_role(user_id, db)
    return role == "admin"