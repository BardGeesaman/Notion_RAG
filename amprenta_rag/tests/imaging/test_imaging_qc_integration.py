"""Integration tests for imaging QC persistence."""

import pytest
from uuid import uuid4
from unittest.mock import MagicMock, patch, AsyncMock
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from amprenta_rag.database.base import Base
from amprenta_rag.imaging.models import ImageQCRecord


@pytest.fixture
def db_session():
    """Create in-memory SQLite database for testing."""
    engine = create_engine("sqlite:///:memory:", echo=False)
    
    # Only create the ImageQCRecord table to avoid ARRAY type issues
    ImageQCRecord.__table__.create(engine, checkfirst=True)
    
    Session = sessionmaker(bind=engine)
    session = Session()
    
    yield session
    
    session.close()


class TestImageQCRecord:
    """Test ImageQCRecord model."""

    def test_qc_record_creation(self, db_session):
        """Test creating a QC record."""
        image_id = uuid4()
        
        record = ImageQCRecord(
            image_id=image_id,
            focus_score=0.85,
            focus_algorithm="laplacian",
            is_focused=True,
            saturation_percent=2.5,
            is_saturated=False,
            uniformity_score=0.92,
            vignetting_detected=False,
            overall_score=87.5,
            passed_qc=True,
            issues=[]
        )
        
        db_session.add(record)
        db_session.commit()
        
        assert record.id is not None
        assert record.focus_score == 0.85
        assert record.passed_qc is True

    def test_qc_record_with_issues(self, db_session):
        """Test QC record with issues list."""
        image_id = uuid4()
        
        record = ImageQCRecord(
            image_id=image_id,
            focus_score=0.3,
            is_focused=False,
            saturation_percent=15.0,
            is_saturated=True,
            overall_score=45.0,
            passed_qc=False,
            issues=["Low focus score", "High saturation detected"]
        )
        
        db_session.add(record)
        db_session.commit()
        db_session.refresh(record)
        
        assert record.passed_qc is False
        assert len(record.issues) == 2
        assert "Low focus score" in record.issues
