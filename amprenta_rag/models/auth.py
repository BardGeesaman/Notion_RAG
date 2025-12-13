"""Authentication and authorization models."""

from __future__ import annotations

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, String, Text, UniqueConstraint, JSON
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import relationship
from sqlalchemy import func

from amprenta_rag.database.base import Base
import uuid


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class User(Base):
    """User account for authentication."""

    __tablename__ = "users"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    username = Column(String(100), unique=True, nullable=False, index=True)
    email = Column(String(255), unique=True, nullable=False)
    password_hash = Column(String(255), nullable=False)
    role = Column(String(50), default="researcher")  # admin, researcher, viewer
    is_active = Column(Boolean, default=True)
    created_at = Column(DateTime, default=func.now())
    last_login = Column(DateTime, nullable=True)


class Team(Base):
    """Team/organization for collaboration."""

    __tablename__ = "teams"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False, unique=True)
    description = Column(Text, nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    members = relationship("TeamMember", back_populates="team", cascade="all, delete-orphan")
    projects = relationship("Project", back_populates="team", cascade="all, delete-orphan")


class TeamMember(Base):
    """Team membership with role."""

    __tablename__ = "team_members"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    team_id = Column(UUID(as_uuid=True), ForeignKey("teams.id"), nullable=False)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    role = Column(String(20), default="member")  # owner, admin, member, viewer
    joined_at = Column(DateTime(timezone=True), server_default=func.now())

    team = relationship("Team", back_populates="members")
    user = relationship("User", backref="team_memberships")


class Project(Base):
    """Project within a team."""

    __tablename__ = "projects"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    team_id = Column(UUID(as_uuid=True), ForeignKey("teams.id"), nullable=False)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    is_public = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    team = relationship("Team", back_populates="projects")
    experiments = relationship("Experiment", backref="project")


class FeaturePermission(Base):
    """Feature visibility permissions by role."""

    __tablename__ = "feature_permissions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    role = Column(String(50), nullable=False)  # admin, researcher, viewer
    page_name = Column(String(100), nullable=False)  # page identifier
    is_visible = Column(Boolean, default=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    __table_args__ = (
        UniqueConstraint("role", "page_name", name="uq_feature_permissions_role_page"),
    )


class AuditLog(Base):
    """Audit trail for user actions."""

    __tablename__ = "audit_logs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    username = Column(String(100), nullable=True)  # Cached for when user deleted
    action = Column(String(100), nullable=False)  # login, logout, create, update, delete
    entity_type = Column(String(100), nullable=True)  # experiment, dataset, signature, etc.
    entity_id = Column(String(100), nullable=True)  # UUID or identifier
    details = Column(JSON, nullable=True)  # Additional context
    ip_address = Column(String(50), nullable=True)
    timestamp = Column(DateTime, default=func.now(), index=True)

    user = relationship("User", backref="audit_logs")


__all__ = [
    "User",
    "Team",
    "TeamMember",
    "Project",
    "FeaturePermission",
    "AuditLog",
]

