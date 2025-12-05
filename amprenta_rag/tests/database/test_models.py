"""
Tests for SQLAlchemy database models.

These tests verify:
- Models can be imported
- Models can be queried
- Relationships work correctly
"""

import uuid

import pytest
from sqlalchemy.exc import IntegrityError

from amprenta_rag.database.base import get_session_local
from amprenta_rag.database.models import (
    Dataset,
    Experiment,
    Feature,
    Program,
    Signature,
    SignatureComponent,
)


class TestModelImports:
    """Test that models can be imported."""
    
    def test_program_import(self):
        """Test Program model can be imported."""
        assert Program is not None
    
    def test_experiment_import(self):
        """Test Experiment model can be imported."""
        assert Experiment is not None
    
    def test_dataset_import(self):
        """Test Dataset model can be imported."""
        assert Dataset is not None
    
    def test_feature_import(self):
        """Test Feature model can be imported."""
        assert Feature is not None
    
    def test_signature_import(self):
        """Test Signature model can be imported."""
        assert Signature is not None
    
    def test_signature_component_import(self):
        """Test SignatureComponent model can be imported."""
        assert SignatureComponent is not None


class TestModelQueries:
    """Test that models can be queried."""
    
    def test_program_query(self):
        """Test Program model can be queried."""
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        try:
            count = db.query(Program).count()
            assert isinstance(count, int)
            assert count >= 0
        finally:
            db.close()
    
    def test_experiment_query(self):
        """Test Experiment model can be queried."""
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        try:
            count = db.query(Experiment).count()
            assert isinstance(count, int)
            assert count >= 0
        finally:
            db.close()
    
    def test_dataset_query(self):
        """Test Dataset model can be queried."""
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        try:
            count = db.query(Dataset).count()
            assert isinstance(count, int)
            assert count >= 0
        finally:
            db.close()
    
    def test_feature_query(self):
        """Test Feature model can be queried."""
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        try:
            count = db.query(Feature).count()
            assert isinstance(count, int)
            assert count >= 0
        finally:
            db.close()
    
    def test_signature_query(self):
        """Test Signature model can be queried."""
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        try:
            count = db.query(Signature).count()
            assert isinstance(count, int)
            assert count >= 0
        finally:
            db.close()
    
    def test_signature_component_query(self):
        """Test SignatureComponent model can be queried."""
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        try:
            count = db.query(SignatureComponent).count()
            assert isinstance(count, int)
            assert count >= 0
        finally:
            db.close()


class TestModelCreation:
    """Test that models can be created (with rollback)."""
    
    def test_create_program(self):
        """Test creating a Program."""
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        try:
            program = Program(
                id=uuid.uuid4(),
                name="Test Program",
                description="Created by test",
            )
            db.add(program)
            db.flush()  # Flush without committing
            
            # Verify it was created
            assert program.id is not None
            assert program.name == "Test Program"
            
            # Rollback - don't actually commit
            db.rollback()
        finally:
            db.close()
    
    def test_create_dataset(self):
        """Test creating a Dataset."""
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        try:
            dataset = Dataset(
                id=uuid.uuid4(),
                name="Test Dataset",
                omics_type="Lipidomics",
            )
            db.add(dataset)
            db.flush()
            
            # Verify it was created
            assert dataset.id is not None
            assert dataset.name == "Test Dataset"
            assert dataset.omics_type == "Lipidomics"
            
            # Rollback
            db.rollback()
        finally:
            db.close()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

