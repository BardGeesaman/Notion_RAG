"""Email notification service."""
from __future__ import annotations

import os
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
from typing import List, Optional, Tuple
from uuid import UUID

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.notifications.email_templates import (
    get_experiment_summary_html,
    get_share_email_html,
)

logger = get_logger(__name__)

# Email configuration from environment variables
SMTP_HOST = os.getenv("SMTP_HOST", "smtp.gmail.com")
SMTP_PORT = int(os.getenv("SMTP_PORT", "587"))
SMTP_USER = os.getenv("SMTP_USER")
SMTP_PASSWORD = os.getenv("SMTP_PASSWORD")
FROM_EMAIL = os.getenv("FROM_EMAIL", "") or (SMTP_USER or "")


def is_email_configured() -> bool:
    """
    Check if email service is configured.

    Returns:
        True if SMTP_USER and SMTP_PASSWORD are set, False otherwise
    """
    return bool(SMTP_USER and SMTP_PASSWORD)


def send_email(
    to: str,
    subject: str,
    body: str,
    html_body: Optional[str] = None,
    attachments: Optional[List[Tuple[str, bytes, str]]] = None,
) -> bool:
    """
    Send an email using SMTP with TLS.

    Args:
        to: Recipient email address
        subject: Email subject line
        body: Plain text email body
        html_body: Optional HTML email body (if provided, email will be multipart/alternative)
        attachments: Optional list of attachments as tuples of (filename, content_bytes, content_type)

    Returns:
        True if email sent successfully, False otherwise

    Note:
        If SMTP is not configured, logs a warning and returns False without crashing.
    """
    if not is_email_configured():
        logger.warning("[EMAIL] Email service not configured. Set SMTP_USER and SMTP_PASSWORD environment variables.")
        return False

    if not to:
        logger.error("[EMAIL] No recipient email address provided")
        return False

    try:
        # Create message
        msg = MIMEMultipart("alternative" if html_body else "mixed")
        msg["From"] = (FROM_EMAIL or SMTP_USER or "")
        msg["To"] = to
        msg["Subject"] = subject

        # Add plain text body
        text_part = MIMEText(body, "plain")
        msg.attach(text_part)

        # Add HTML body if provided
        if html_body:
            html_part = MIMEText(html_body, "html")
            msg.attach(html_part)

        # Add attachments if provided
        if attachments:
            for filename, content_bytes, content_type in attachments:
                maintype, subtype = content_type.split("/", 1) if "/" in content_type else ("application", "octet-stream")
                attachment = MIMEBase(maintype, subtype)
                attachment.set_payload(content_bytes)
                encoders.encode_base64(attachment)
                attachment.add_header(
                    "Content-Disposition",
                    f'attachment; filename="{filename}"',
                )
                msg.attach(attachment)

        # Send email via SMTP with TLS
        with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
            server.starttls()
            server.login(SMTP_USER or "", SMTP_PASSWORD or "")
            server.send_message(msg)

        logger.info("[EMAIL] Email sent successfully to %s", to)
        return True

    except smtplib.SMTPException as e:
        logger.error("[EMAIL] SMTP error sending email to %s: %s", to, e)
        return False
    except Exception as e:
        logger.error("[EMAIL] Unexpected error sending email to %s: %s", to, e)
        return False


def send_experiment_summary(experiment_id: UUID, to: str, db) -> bool:
    """
    Send an experiment summary email.

    Args:
        experiment_id: UUID of the experiment
        to: Recipient email address
        db: Database session

    Returns:
        True if email sent successfully, False otherwise
    """
    try:
        from amprenta_rag.database.models import Experiment

        # Load experiment from DB
        experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()
        if not experiment:
            logger.error("[EMAIL] Experiment %s not found", experiment_id)
            return False

        # Get related datasets through the relationship
        datasets = list(experiment.datasets) if hasattr(experiment, "datasets") else []

        # Generate HTML
        html_body = get_experiment_summary_html(experiment, datasets)

        # Generate plain text fallback
        plain_body = f"""
Experiment Summary: {experiment.name}

Description: {experiment.description or 'No description available'}
Design Type: {experiment.design_type or 'Unknown'}
Organism: {experiment.organism if hasattr(experiment, 'organism') and experiment.organism else 'Unknown'}
Datasets: {len(datasets)} dataset(s) associated
        """.strip()

        # Send email
        subject = f"Experiment Summary: {experiment.name}"
        return send_email(to, subject, plain_body, html_body=html_body)

    except Exception as e:
        logger.error("[EMAIL] Error sending experiment summary: %s", e)
        return False


def send_share_notification(
    entity_type: str,
    entity_id: UUID,
    to: str,
    from_user: str,
    message: Optional[str],
    db,
) -> bool:
    """
    Send a share notification email.

    Args:
        entity_type: Type of entity being shared (e.g., "Experiment", "Compound")
        entity_id: UUID of the entity
        to: Recipient email address
        from_user: Username of the person sharing
        message: Optional message from the sharer
        db: Database session

    Returns:
        True if email sent successfully, False otherwise
    """
    try:
        from amprenta_rag.database.models import Experiment, Compound

        # Load entity from DB based on type
        entity = None
        entity_name = "Unknown"

        if entity_type.lower() == "experiment":
            entity = db.query(Experiment).filter(Experiment.id == entity_id).first()
            if entity:
                entity_name = entity.name
        elif entity_type.lower() == "compound":
            entity = db.query(Compound).filter(Compound.id == entity_id).first()
            if entity:
                entity_name = entity.compound_id or entity.name or "Unknown Compound"
        else:
            logger.warning("[EMAIL] Unknown entity type: %s", entity_type)
            entity_name = f"{entity_type} {entity_id}"

        # Generate HTML
        html_body = get_share_email_html(entity_type, entity_name, from_user, message)

        # Generate plain text fallback
        plain_body = f"""
{from_user} has shared a {entity_type} with you: {entity_name}
        """.strip()
        if message:
            plain_body += f"\n\nMessage: {message}"

        # Send email
        subject = f"Shared {entity_type}: {entity_name}"
        return send_email(to, subject, plain_body, html_body=html_body)

    except Exception as e:
        logger.error("[EMAIL] Error sending share notification: %s", e)
        return False

