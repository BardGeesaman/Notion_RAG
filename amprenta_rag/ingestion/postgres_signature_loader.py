"""
Postgres-based signature loading.

Loads signatures from Postgres database instead of Notion.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import (
    Signature as SignatureModel,
    SignatureComponent as SignatureComponentModel,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature, SignatureComponent

logger = get_logger(__name__)


def load_signature_from_postgres(signature_model: SignatureModel) -> Optional[Signature]:
    """
    Load a Signature object from a Postgres SignatureModel.

    Args:
        signature_model: SignatureModel instance from Postgres

    Returns:
        Signature object or None if loading failed
    """
    try:
        # Load signature components
        components: List[SignatureComponent] = []

        with db_session() as db:
            try:
                sig_components = (
                    db.query(SignatureComponentModel)
                    .filter(SignatureComponentModel.signature_id == signature_model.id)
                    .all()
                )

                for comp_model in sig_components:
                    component = SignatureComponent(
                        feature_name=comp_model.feature_name or "",
                        feature_type=comp_model.feature_type,
                        direction=comp_model.direction,
                        weight=float(comp_model.weight) if comp_model.weight is not None else 1.0,
                    )
                    components.append(component)
            except Exception as e:
                logger.error(
                    "[POSTGRES-SIGNATURE] Error loading components for signature %s: %r",
                    signature_model.id,
                    e,
                )
                return None

        # Build Signature object
        signature = Signature(
            name=signature_model.name or "",
            components=components,
            modalities=list(signature_model.modalities or []),
            description=signature_model.description,
        )

        logger.debug(
            "[POSTGRES-SIGNATURE] Loaded signature %s with %d components",
            signature_model.name,
            len(components),
        )

        return signature

    except Exception as e:
        logger.error(
            "[POSTGRES-SIGNATURE] Error loading signature %s from Postgres: %r",
            signature_model.id,
            e,
        )
        return None


def fetch_all_signatures_from_postgres() -> List[SignatureModel]:
    """
    Fetch all Signature models from Postgres.

    Returns:
        List of SignatureModel instances
    """
    with db_session() as db:
        try:
            signatures = db.query(SignatureModel).all()

            logger.debug(
                "[POSTGRES-SIGNATURE] Fetched %d signature(s) from Postgres",
                len(signatures),
            )

            return signatures

        except Exception as e:
            logger.error(
                "[POSTGRES-SIGNATURE] Error fetching signatures from Postgres: %r",
                e,
            )
            return []


def find_signature_by_id(signature_id: UUID) -> Optional[SignatureModel]:
    """
    Find a signature in Postgres by ID.

    Args:
        signature_id: Postgres UUID of the signature

    Returns:
        SignatureModel instance or None if not found
    """
    with db_session() as db:
        return (
            db.query(SignatureModel)
            .filter(SignatureModel.id == signature_id)
            .first()
        )


def find_signatures_by_name(name: str) -> List[SignatureModel]:
    """
    Find signatures in Postgres by name.

    Args:
        name: Signature name to search for

    Returns:
        List of SignatureModel instances matching the name
    """
    with db_session() as db:
        return (
            db.query(SignatureModel)
            .filter(SignatureModel.name.ilike(f"%{name}%"))
            .all()
        )

