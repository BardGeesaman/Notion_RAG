"""
Gmail API client for fetching emails directly from Gmail.

Replaces Zapier workflow by fetching emails directly from Gmail inbox
and ingesting them using Postgres-only ingestion.
"""

from __future__ import annotations

import base64
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Try to import Gmail API libraries
try:
    from google.auth.transport.requests import Request
    from google.oauth2.credentials import Credentials
    from google_auth_oauthlib.flow import InstalledAppFlow
    from googleapiclient.discovery import build
    from googleapiclient.errors import HttpError

    GMAIL_AVAILABLE = True
except ImportError:
    GMAIL_AVAILABLE = False
    logger.warning(
        "Gmail API libraries not installed. Install with: pip install google-auth google-auth-oauthlib google-auth-httplib2 google-api-python-client"
    )


# Gmail API scopes needed
SCOPES = ["https://www.googleapis.com/auth/gmail.readonly"]


class GmailClient:
    """Client for fetching emails from Gmail API."""

    def __init__(
        self,
        credentials_file: Optional[str] = None,
        token_file: Optional[str] = None,
    ):
        """
        Initialize Gmail client.
        
        Args:
            credentials_file: Path to OAuth2 credentials JSON file
            token_file: Path to store/load OAuth2 token JSON file
        """
        if not GMAIL_AVAILABLE:
            raise ImportError(
                "Gmail API libraries not installed. "
                "Install with: pip install google-auth google-auth-oauthlib google-auth-httplib2 google-api-python-client"
            )
        
        self.credentials_file = credentials_file or os.getenv(
            "GMAIL_CREDENTIALS_FILE", "credentials/gmail_credentials.json"
        )
        self.token_file = token_file or os.getenv(
            "GMAIL_TOKEN_FILE", "credentials/gmail_token.json"
        )
        
        self.service = None
        self._authenticate()

    def _authenticate(self) -> None:
        """Authenticate with Gmail API using OAuth2."""
        creds = None
        
        # Load existing token if available
        token_path = Path(self.token_file)
        if token_path.exists():
            try:
                creds = Credentials.from_authorized_user_file(str(token_path), SCOPES)
                logger.info("[GMAIL] Loaded existing credentials from token file")
            except Exception as e:
                logger.warning(
                    "[GMAIL] Could not load existing token file: %r", e
                )
        
        # If no valid credentials, authenticate
        if not creds or not creds.valid:
            if creds and creds.expired and creds.refresh_token:
                try:
                    creds.refresh(Request())
                    logger.info("[GMAIL] Refreshed expired credentials")
                except Exception as e:
                    logger.warning(
                        "[GMAIL] Could not refresh credentials: %r", e
                    )
                    creds = None
            
            if not creds:
                # Start OAuth flow
                creds_path = Path(self.credentials_file)
                if not creds_path.exists():
                    raise FileNotFoundError(
                        f"Gmail credentials file not found: {self.credentials_file}\n"
                        "Please download OAuth2 credentials from Google Cloud Console:\n"
                        "1. Go to https://console.cloud.google.com/\n"
                        "2. Create/select a project\n"
                        "3. Enable Gmail API\n"
                        "4. Create OAuth 2.0 credentials\n"
                        "5. Save as credentials/gmail_credentials.json"
                    )
                
                flow = InstalledAppFlow.from_client_secrets_file(
                    str(creds_path), SCOPES
                )
                creds = flow.run_local_server(port=0)
                logger.info("[GMAIL] Completed OAuth2 authentication")
            
            # Save token for future use
            token_path.parent.mkdir(parents=True, exist_ok=True)
            with open(token_path, "w") as token:
                token.write(creds.to_json())
            logger.info("[GMAIL] Saved credentials to token file")
        
        # Build Gmail API service
        try:
            self.service = build("gmail", "v1", credentials=creds)
            logger.info("[GMAIL] Gmail API service initialized")
        except Exception as e:
            logger.error("[GMAIL] Error building Gmail service: %r", e)
            raise

    def fetch_emails(
        self,
        query: str = "in:inbox",
        max_results: int = 100,
        since: Optional[datetime] = None,
    ) -> List[Dict[str, Any]]:
        """
        Fetch emails from Gmail.
        
        Args:
            query: Gmail search query (e.g., "in:inbox", "from:example@gmail.com")
            max_results: Maximum number of emails to fetch
            since: Only fetch emails since this datetime
            
        Returns:
            List of email dictionaries with parsed email data
        """
        if not self.service:
            raise RuntimeError("Gmail service not initialized. Call _authenticate() first.")
        
        # Build query with date filter if provided
        if since:
            since_str = since.strftime("%Y/%m/%d")
            query = f"{query} after:{since_str}"
        
        logger.info("[GMAIL] Fetching emails with query: %s", query)
        
        try:
            # List messages
            results = (
                self.service.users()
                .messages()
                .list(userId="me", q=query, maxResults=max_results)
                .execute()
            )
            
            messages = results.get("messages", [])
            logger.info("[GMAIL] Found %d message(s)", len(messages))
            
            if not messages:
                return []
            
            # Fetch full message details
            emails = []
            for msg in messages:
                try:
                    email_data = self._parse_email_message(msg["id"])
                    if email_data:
                        emails.append(email_data)
                except Exception as e:
                    logger.warning(
                        "[GMAIL] Error parsing message %s: %r", msg.get("id"), e
                    )
                    continue
            
            logger.info("[GMAIL] Successfully parsed %d email(s)", len(emails))
            return emails
            
        except HttpError as e:
            logger.error("[GMAIL] Gmail API error: %r", e)
            raise

    def _parse_email_message(self, message_id: str) -> Optional[Dict[str, Any]]:
        """
        Parse a single email message into structured data.
        
        Args:
            message_id: Gmail message ID
            
        Returns:
            Dictionary with email data or None if parsing fails
        """
        try:
            message = (
                self.service.users()
                .messages()
                .get(userId="me", id=message_id, format="full")
                .execute()
            )
            
            # Extract headers
            headers = message.get("payload", {}).get("headers", [])
            header_dict = {h["name"].lower(): h["value"] for h in headers}
            
            # Extract body
            body_text = self._extract_body(message.get("payload", {}))
            
            if not body_text:
                logger.debug("[GMAIL] No text body found for message %s", message_id)
                return None
            
            # Parse email data
            email_data = {
                "id": message_id,
                "thread_id": message.get("threadId"),
                "subject": header_dict.get("subject", "(no subject)"),
                "from": header_dict.get("from", ""),
                "to": header_dict.get("to", ""),
                "cc": header_dict.get("cc", ""),
                "date": header_dict.get("date", ""),
                "body": body_text,
                "snippet": message.get("snippet", ""),
                "labels": message.get("labelIds", []),
            }
            
            return email_data
            
        except Exception as e:
            logger.warning(
                "[GMAIL] Error parsing message %s: %r", message_id, e
            )
            return None

    def _extract_body(self, payload: Dict[str, Any]) -> str:
        """
        Extract text body from email payload.
        
        Args:
            payload: Email payload dictionary
            
        Returns:
            Email body text
        """
        body_text = ""
        
        # Check if this is a multipart message
        if payload.get("mimeType") == "multipart/alternative" or payload.get("mimeType") == "multipart/mixed":
            parts = payload.get("parts", [])
            
            # Prefer text/plain, fallback to text/html
            for part in parts:
                mime_type = part.get("mimeType", "")
                if mime_type == "text/plain":
                    body_text = self._decode_body(part)
                    break
                elif mime_type == "text/html" and not body_text:
                    body_text = self._decode_body(part)
            
            # If no text found, recursively check nested parts
            if not body_text:
                for part in parts:
                    nested_parts = part.get("parts", [])
                    if nested_parts:
                        nested_body = self._extract_body(part)
                        if nested_body:
                            body_text = nested_body
                            break
        else:
            # Single part message
            mime_type = payload.get("mimeType", "")
            if mime_type in ["text/plain", "text/html"]:
                body_text = self._decode_body(payload)
        
        return body_text.strip()

    def _decode_body(self, part: Dict[str, Any]) -> str:
        """
        Decode email body from base64.
        
        Args:
            part: Email part dictionary
            
        Returns:
            Decoded body text
        """
        body_data = part.get("body", {}).get("data")
        if not body_data:
            return ""
        
        try:
            # Gmail API uses URL-safe base64 encoding
            decoded = base64.urlsafe_b64decode(body_data)
            return decoded.decode("utf-8", errors="ignore")
        except Exception as e:
            logger.debug("[GMAIL] Error decoding body: %r", e)
            return ""

