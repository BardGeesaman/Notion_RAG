"""Authentication and authorization models."""

from __future__ import annotations

import uuid
from typing import List, Optional, TYPE_CHECKING

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, String, Text, UniqueConstraint, JSON, func, Integer, Float
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import Experiment


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
    company_id = Column(UUID(as_uuid=True), ForeignKey("companies.id"), nullable=True, index=True)
    company_role = Column(String(20), nullable=True, default="member")  # member|admin|owner (tenant role)
    is_active = Column(Boolean, default=True)
    created_at = Column(DateTime, default=func.now())
    last_login = Column(DateTime, nullable=True)

    company: Mapped[Optional["Company"]] = relationship("Company", back_populates="users", foreign_keys=[company_id])


class Company(Base):
    """Tenant organization for multi-tenancy."""

    __tablename__ = "companies"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False, index=True)
    subdomain = Column(String(100), nullable=False, unique=True, index=True)
    logo_url = Column(String(500), nullable=True)
    primary_color = Column(String(50), nullable=True)

    max_users = Column(Integer, nullable=True)
    max_datasets = Column(Integer, nullable=True)
    storage_quota_gb = Column(Float, nullable=True)

    status = Column(String(50), nullable=True, default="active", index=True)
    trial_ends_at = Column(DateTime(timezone=True), nullable=True)

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)

    users: Mapped[List["User"]] = relationship("User", back_populates="company")


class Team(Base):
    """Team/organization for collaboration."""

    __tablename__ = "teams"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False, unique=True)
    description = Column(Text, nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    members: Mapped[List["TeamMember"]] = relationship(
        back_populates="team", cascade="all, delete-orphan"
    )
    projects: Mapped[List["Project"]] = relationship(
        back_populates="team", cascade="all, delete-orphan"
    )


class TeamMember(Base):
    """Team membership with role."""

    __tablename__ = "team_members"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    team_id = Column(UUID(as_uuid=True), ForeignKey("teams.id"), nullable=False)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    role = Column(String(20), default="member")  # owner, admin, member, viewer
    joined_at = Column(DateTime(timezone=True), server_default=func.now())

    team: Mapped["Team"] = relationship(back_populates="members")
    user: Mapped["User"] = relationship(backref="team_memberships")


class Project(Base):
    """Project within a team."""

    __tablename__ = "projects"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    team_id = Column(UUID(as_uuid=True), ForeignKey("teams.id"), nullable=False)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    is_public = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    team: Mapped["Team"] = relationship(back_populates="projects")
    experiments: Mapped[List["Experiment"]] = relationship(back_populates="project")


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

    user: Mapped[Optional["User"]] = relationship(backref="audit_logs")


__all__ = [
    "Company",
    "User",
    "Team",
    "TeamMember",
    "Project",
    "FeaturePermission",
    "AuditLog",
]

