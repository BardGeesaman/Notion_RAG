"""Company (tenant) management API endpoints."""

from __future__ import annotations

import secrets
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user_with_company, get_database_session
from amprenta_rag.auth.password import hash_password
from amprenta_rag.database.models import Company, User


router = APIRouter(prefix="/companies", tags=["Companies"])


def _require_superadmin(current_user: User) -> None:
    # Global user.role is used across the app for permissions (admin/researcher/viewer).
    if (getattr(current_user, "role", None) or "").lower() != "admin":
        raise HTTPException(status_code=403, detail="Superadmin only")


class CompanyCreateRequest(BaseModel):
    name: str
    subdomain: str
    logo_url: Optional[str] = None
    primary_color: Optional[str] = None
    max_users: Optional[int] = None
    max_datasets: Optional[int] = None
    storage_quota_gb: Optional[float] = None
    status: Optional[str] = "active"


class CompanyPatchRequest(BaseModel):
    name: Optional[str] = None
    logo_url: Optional[str] = None
    primary_color: Optional[str] = None
    max_users: Optional[int] = None
    max_datasets: Optional[int] = None
    storage_quota_gb: Optional[float] = None
    status: Optional[str] = None


class CompanyResponse(BaseModel):
    id: UUID
    name: str
    subdomain: str
    logo_url: Optional[str] = None
    primary_color: Optional[str] = None
    max_users: Optional[int] = None
    max_datasets: Optional[int] = None
    storage_quota_gb: Optional[float] = None
    status: Optional[str] = None

    class Config:
        from_attributes = True


class UserResponse(BaseModel):
    id: UUID
    username: str
    email: str
    role: Optional[str] = None
    company_id: Optional[UUID] = None
    company_role: Optional[str] = None
    is_active: Optional[bool] = None

    class Config:
        from_attributes = True


class InviteUserRequest(BaseModel):
    email: str
    role: str = Field("member", description="Tenant role: member|admin|owner (stored in company_role)")


class UpdateUserRoleRequest(BaseModel):
    role: str = Field(..., description="Tenant role: member|admin|owner (stored in company_role)")


@router.post("", response_model=CompanyResponse)
def create_company(
    payload: CompanyCreateRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> CompanyResponse:
    _require_superadmin(current_user)
    existing = db.query(Company).filter(Company.subdomain == payload.subdomain).first()
    if existing:
        raise HTTPException(status_code=400, detail="subdomain already exists")
    c = Company(
        name=payload.name,
        subdomain=payload.subdomain,
        logo_url=payload.logo_url,
        primary_color=payload.primary_color,
        max_users=payload.max_users,
        max_datasets=payload.max_datasets,
        storage_quota_gb=payload.storage_quota_gb,
        status=payload.status or "active",
    )
    db.add(c)
    db.commit()
    db.refresh(c)
    return CompanyResponse.model_validate(c)


@router.get("/{company_id}", response_model=CompanyResponse)
def get_company(
    company_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> CompanyResponse:
    _require_superadmin(current_user)
    c = db.query(Company).filter(Company.id == company_id).first()
    if not c:
        raise HTTPException(status_code=404, detail="Company not found")
    return CompanyResponse.model_validate(c)


@router.patch("/{company_id}", response_model=CompanyResponse)
def patch_company(
    company_id: UUID,
    payload: CompanyPatchRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> CompanyResponse:
    _require_superadmin(current_user)
    c = db.query(Company).filter(Company.id == company_id).first()
    if not c:
        raise HTTPException(status_code=404, detail="Company not found")
    for k, v in payload.model_dump(exclude_unset=True).items():
        setattr(c, k, v)
    db.add(c)
    db.commit()
    db.refresh(c)
    return CompanyResponse.model_validate(c)


@router.get("/{company_id}/users", response_model=List[UserResponse])
def list_company_users(
    company_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> List[UserResponse]:
    _require_superadmin(current_user)
    users = db.query(User).filter(User.company_id == company_id).order_by(User.username.asc()).all()
    return [UserResponse.model_validate(u) for u in users]


@router.post("/{company_id}/users/invite")
def invite_user(
    company_id: UUID,
    payload: InviteUserRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> dict:
    _require_superadmin(current_user)
    comp = db.query(Company).filter(Company.id == company_id).first()
    if not comp:
        raise HTTPException(status_code=404, detail="Company not found")

    email = payload.email.strip().lower()
    if not email or "@" not in email:
        raise HTTPException(status_code=400, detail="Invalid email")

    user = db.query(User).filter(User.email == email).first()
    temp_password = None
    if user is None:
        base = email.split("@", 1)[0]
        username = base
        i = 1
        while db.query(User).filter(User.username == username).first() is not None:
            i += 1
            username = f"{base}{i}"
        temp_password = secrets.token_urlsafe(16)
        user = User(
            username=username,
            email=email,
            password_hash=hash_password(temp_password),
            role="researcher",
            is_active=True,
        )
        db.add(user)
        db.flush()

    user.company_id = comp.id
    user.company_role = payload.role
    db.add(user)
    db.commit()
    db.refresh(user)

    out = {"user": UserResponse.model_validate(user).model_dump()}
    if temp_password:
        out["temp_password"] = temp_password
    return out


@router.delete("/{company_id}/users/{user_id}")
def remove_user_from_company(
    company_id: UUID,
    user_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> dict:
    _require_superadmin(current_user)
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    if user.company_id != company_id:
        raise HTTPException(status_code=400, detail="User not in this company")
    user.company_id = None
    user.company_role = None
    db.add(user)
    db.commit()
    return {"removed": True}


@router.patch("/{company_id}/users/{user_id}", response_model=UserResponse)
def update_company_user_role(
    company_id: UUID,
    user_id: UUID,
    payload: UpdateUserRoleRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> UserResponse:
    _require_superadmin(current_user)
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    if user.company_id != company_id:
        raise HTTPException(status_code=400, detail="User not in this company")
    user.company_role = payload.role
    db.add(user)
    db.commit()
    db.refresh(user)
    return UserResponse.model_validate(user)


__all__ = ["router"]


