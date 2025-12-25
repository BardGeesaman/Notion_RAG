"""Notifications for scheduled digests (email + Slack)."""

from __future__ import annotations

import os
import smtplib
from email.message import EmailMessage
from typing import Iterable, List, Optional

import httpx


def _split_recipients(recipients: Iterable[str]) -> List[str]:
    out: List[str] = []
    for r in recipients or []:
        s = str(r).strip()
        if not s:
            continue
        out.append(s)
    return out


def _send_email(to_emails: List[str], subject: str, body: str) -> None:
    host = os.environ.get("SMTP_HOST", "").strip()
    port = int(os.environ.get("SMTP_PORT", "587"))
    user = os.environ.get("SMTP_USER", "").strip()
    password = os.environ.get("SMTP_PASS", "").strip()
    from_email = os.environ.get("SMTP_FROM", user or "no-reply@amprenta.local").strip()

    if not host:
        raise RuntimeError("SMTP_HOST is not configured")
    if not to_emails:
        return

    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"] = from_email
    msg["To"] = ", ".join(to_emails)
    msg.set_content(body)

    with smtplib.SMTP(host, port, timeout=30) as s:
        s.starttls()
        if user and password:
            s.login(user, password)
        s.send_message(msg)


def _send_slack(message: str, *, webhook_url: str) -> None:
    if not webhook_url:
        raise RuntimeError("Slack webhook URL missing")
    with httpx.Client(timeout=30) as client:
        r = client.post(webhook_url, json={"text": message})
    r.raise_for_status()


def send_digest_notification(
    recipients: List[str],
    digest_url: str,
    program_name: str,
    *,
    slack_webhook_url: Optional[str] = None,
) -> dict:
    """Send a digest-ready notification.

    recipients: list of emails (and/or identifiers; non-emails are ignored for SMTP)
    """
    rs = _split_recipients(recipients)
    subject = f"Weekly {program_name} digest is ready"
    body = f"Weekly {program_name} digest is ready: {digest_url}"

    emailed: List[str] = [r for r in rs if "@" in r]
    slack_url = (slack_webhook_url or os.environ.get("SLACK_WEBHOOK_URL", "")).strip() or None

    out = {"emailed": emailed, "slack": bool(slack_url)}

    if emailed:
        _send_email(emailed, subject, body)
    if slack_url:
        _send_slack(body, webhook_url=slack_url)

    return out


__all__ = ["send_digest_notification"]


