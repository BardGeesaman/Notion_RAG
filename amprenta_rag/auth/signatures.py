"""Electronic signature service for regulatory compliance."""

from __future__ import annotations

import hashlib
import hmac
import json
import os
from typing import List, Optional

from sqlalchemy.orm import Session

from amprenta_rag.auth.password import verify_password
from amprenta_rag.database.models import ElectronicSignature, User
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Secret key for HMAC (should be in env var in production)
SECRET_KEY = os.getenv("SIGNATURE_SECRET_KEY", "default-secret-key-change-in-production").encode()


def create_signature(
    user_id: str,
    action: str,
    entity_type: str,
    entity_id: str,
    password: str,
    meaning: str,
    db: Session,
    ip_address: Optional[str] = None,
) -> Optional[ElectronicSignature]:
    """
    Create electronic signature after password verification.

    Args:
        user_id: User UUID
        action: Action being signed (approve, reject, etc.)
        entity_type: Type of entity
        entity_id: Entity UUID
        password: User's password for verification
        meaning: Statement of signature meaning
        db: Database session
        ip_address: Optional IP address

    Returns:
        ElectronicSignature if successful, None if password wrong
    """
    # Verify user password
    user = db.query(User).filter(User.id == user_id).first()
    
    if not user:
        logger.warning("[SIGNATURE] User not found: %s", user_id)
        return None
    
    if not verify_password(password, user.password_hash):
        logger.warning("[SIGNATURE] Password verification failed for user: %s", user.username)
        return None
    
    # Generate signature hash (HMAC of signing data)
    signing_data = {
        "user_id": str(user_id),
        "action": action,
        "entity_type": entity_type,
        "entity_id": str(entity_id),
        "meaning": meaning,
    }
    
    data_str = json.dumps(signing_data, sort_keys=True)
    signature_hash = hmac.new(SECRET_KEY, data_str.encode(), hashlib.sha256).hexdigest()
    
    # Create signature record
    signature = ElectronicSignature(
        user_id=user_id,
        action=action,
        entity_type=entity_type,
        entity_id=entity_id,
        signature_hash=signature_hash,
        meaning=meaning,
        ip_address=ip_address,
    )
    
    db.add(signature)
    db.commit()
    db.refresh(signature)
    
    logger.info("[SIGNATURE] Created signature for %s by %s", entity_type, user.username)
    
    return signature


def verify_signature(signature_id: str, db: Session) -> bool:
    """
    Verify signature integrity.

    Args:
        signature_id: Signature UUID
        db: Database session

    Returns:
        True if signature is valid and untampered
    """
    signature = db.query(ElectronicSignature).filter(ElectronicSignature.id == signature_id).first()
    
    if not signature:
        return False
    
    # Regenerate expected hash
    signing_data = {
        "user_id": str(signature.user_id),
        "action": signature.action,
        "entity_type": signature.entity_type,
        "entity_id": str(signature.entity_id),
        "meaning": signature.meaning,
    }
    
    data_str = json.dumps(signing_data, sort_keys=True)
    expected_hash = hmac.new(SECRET_KEY, data_str.encode(), hashlib.sha256).hexdigest()
    
    is_valid = hmac.compare_digest(signature.signature_hash, expected_hash)
    
    if is_valid and not signature.verified_at:
        # Mark as verified
        signature.verified_at = db.func.now()
        db.commit()
    
    return is_valid


def get_signatures(entity_type: str, entity_id: str, db: Session) -> List[ElectronicSignature]:
    """
    Get all signatures for an entity.

    Args:
        entity_type: Type of entity
        entity_id: Entity UUID
        db: Database session

    Returns:
        List of ElectronicSignature records
    """
    return (
        db.query(ElectronicSignature)
        .filter(
            ElectronicSignature.entity_type == entity_type,
            ElectronicSignature.entity_id == entity_id,
        )
        .order_by(ElectronicSignature.timestamp.desc())
        .all()
    )

