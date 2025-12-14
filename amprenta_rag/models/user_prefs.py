"""User preference, bookmarking, and comment models."""

from __future__ import annotations

import uuid
from datetime import datetime

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, JSON, String, Text, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import relationship, backref

from amprenta_rag.database.base import Base


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class Feedback(Base):
    """User feedback and feature requests."""

    __tablename__ = "feedback"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    feedback_type = Column(String(20), nullable=False)  # bug, feature
    title = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    status = Column(String(20), default="new")  # new, reviewing, planned, done, wontfix
    priority = Column(String(20), default="medium")  # low, medium, high, critical
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())

    user = relationship("User", backref="feedback")


class UserFavorite(Base):
    """User favorite pages/bookmarks."""

    __tablename__ = "user_favorites"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    page_name = Column(String(100), nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", backref="favorites")


class Notification(Base):
    """User notifications."""

    __tablename__ = "notifications"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    title = Column(String(255), nullable=False)
    message = Column(Text, nullable=True)
    notification_type = Column(String(50), default="info")  # info, success, warning, discovery
    is_read = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", backref="notifications")


class EmailSubscription(Base):
    """Email subscription preferences for users."""

    __tablename__ = "email_subscriptions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False, index=True)
    subscription_type = Column(String(50), nullable=False)  # digest, alerts, shares
    frequency = Column(String(50), nullable=False)  # daily, weekly, immediate
    is_active = Column(Boolean, default=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", backref="email_subscriptions")


class Bookmark(Base):
    """User bookmarks for entities."""

    __tablename__ = "bookmarks"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    entity_type = Column(String(50), nullable=False)  # experiment, compound, signature
    entity_id = Column(UUID(as_uuid=True), nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", backref="bookmarks")


class Note(Base):
    """User notes attached to entities."""

    __tablename__ = "notes"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    entity_type = Column(String(50), nullable=False)  # experiment, compound, signature, dataset, etc.
    entity_id = Column(UUID(as_uuid=True), nullable=False)
    annotation_type = Column(String(50), nullable=True)  # Optional type/category for annotations
    content = Column(Text, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by = relationship("User", foreign_keys=[created_by_id])


class SavedFilter(Base):
    """Saved filter presets for entity lists."""

    __tablename__ = "saved_filters"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    entity_type = Column(String(50), nullable=False)  # experiment, compound, etc.
    filters = Column(JSON, nullable=False)  # Filter criteria as dict
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    user = relationship("User", foreign_keys=[user_id])


class Comment(Base):
    """Contextual comments on entities (experiments, datasets, signatures, compounds)."""

    __tablename__ = "comments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    entity_type = Column(String(50), nullable=False, index=True)  # "experiment", "dataset", "signature", "compound"
    entity_id = Column(UUID(as_uuid=True), nullable=False, index=True)
    parent_id = Column(UUID(as_uuid=True), ForeignKey("comments.id"), nullable=True, index=True)  # For replies
    content = Column(Text, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=True)

    replies = relationship("Comment", backref=backref("parent", remote_side=[id]))
    created_by = relationship("User", foreign_keys=[created_by_id])


__all__ = [
    "Feedback",
    "UserFavorite",
    "Notification",
    "EmailSubscription",
    "Bookmark",
    "Note",
    "SavedFilter",
    "Comment",
]

