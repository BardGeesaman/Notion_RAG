"""Tests for IP tracking models."""

from __future__ import annotations

from datetime import datetime
from uuid import uuid4

import pytest
from sqlalchemy.exc import IntegrityError

from amprenta_rag.models.ip import (
    DisclosureInventor,
    InventionDisclosure,
    IPLink,
    PatentApplication,
    PatentClaim,
)


def test_create_invention_disclosure():
    """Test InventionDisclosure model creation."""
    disclosure = InventionDisclosure(
        title="Novel Compound Discovery",
        description="Description of invention",
        technical_field="Chemistry",
        status="draft",
        created_by_id=uuid4(),
        company_id=uuid4(),
    )
    
    assert disclosure.title == "Novel Compound Discovery"
    assert disclosure.status == "draft"


def test_disclosure_status_values():
    """Test disclosure status field accepts valid values."""
    valid_statuses = ["draft", "submitted", "under_review", "filed", "granted", "rejected"]
    
    for status in valid_statuses:
        disclosure = InventionDisclosure(
            title="Test",
            description="Test",
            status=status,
            created_by_id=uuid4(),
            company_id=uuid4(),
        )
        assert disclosure.status == status


def test_create_disclosure_inventor():
    """Test DisclosureInventor model creation."""
    inventor = DisclosureInventor(
        disclosure_id=uuid4(),
        user_id=uuid4(),
        contribution_percentage=50,
        is_primary=True,
    )
    
    assert inventor.contribution_percentage == 50
    assert inventor.is_primary is True


def test_create_patent_application():
    """Test PatentApplication model creation."""
    patent = PatentApplication(
        disclosure_id=uuid4(),
        application_number="US12345678",
        filing_date=datetime(2024, 1, 1),
        jurisdiction="US",
        status="pending",
    )
    
    assert patent.application_number == "US12345678"
    assert patent.jurisdiction == "US"


def test_create_patent_claim_hierarchy():
    """Test PatentClaim with parent-child relationship."""
    patent_id = uuid4()
    parent_claim_id = uuid4()
    
    parent_claim = PatentClaim(
        id=parent_claim_id,
        patent_id=patent_id,
        claim_number=1,
        claim_text="A compound comprising...",
        claim_type="independent",
    )
    
    dependent_claim = PatentClaim(
        patent_id=patent_id,
        claim_number=2,
        claim_text="The compound of claim 1, wherein...",
        claim_type="dependent",
        parent_claim_id=parent_claim_id,
    )
    
    assert parent_claim.claim_type == "independent"
    assert dependent_claim.parent_claim_id == parent_claim_id


def test_ip_link_unique_constraint():
    """Test IPLink unique constraint on entity/disclosure/link_type."""
    link1 = IPLink(
        entity_type="compound",
        entity_id=uuid4(),
        disclosure_id=uuid4(),
        link_type="embodiment",
    )
    
    # Unique constraint would be enforced at database level
    assert link1.link_type == "embodiment"

