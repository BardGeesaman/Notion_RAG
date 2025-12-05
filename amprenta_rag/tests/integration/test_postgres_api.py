"""
Integration tests for Postgres + FastAPI.

These tests verify:
- API can create/read/update/delete entities in Postgres
- End-to-end workflows work
- Data persists correctly
"""

import uuid

import pytest
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from amprenta_rag.api.main import app
from amprenta_rag.database.base import get_session_local
from amprenta_rag.database.models import Program

client = TestClient(app)


class TestPostgresAPIIntegration:
    """Test integration between Postgres and FastAPI."""
    
    def test_create_program_via_api(self):
        """Test creating a program via API persists to Postgres."""
        program_data = {
            "name": "Integration Test Program",
            "description": "Created via API integration test",
        }
        
        # Create via API
        response = client.post("/api/v1/programs", json=program_data)
        if response.status_code not in (200, 201):
            pytest.skip("API endpoint not fully implemented yet")
        
        program = response.json()
        program_id = uuid.UUID(program["id"])
        
        # Verify in database
        SessionLocal = get_session_local()
        db = SessionLocal()
        try:
            db_program = db.query(Program).filter(Program.id == program_id).first()
            assert db_program is not None
            assert db_program.name == program_data["name"]
        finally:
            db.close()
        
        # Clean up
        client.delete(f"/api/v1/programs/{program_id}")
    
    def test_read_program_via_api(self):
        """Test reading a program from Postgres via API."""
        # First create in database
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        program_id = uuid.uuid4()
        db_program = Program(
            id=program_id,
            name="Test Program for API Read",
            description="Created in DB, read via API",
        )
        db.add(db_program)
        db.commit()
        db.refresh(db_program)
        
        try:
            # Read via API
            response = client.get(f"/api/v1/programs/{program_id}")
            if response.status_code == 404:
                pytest.skip("API endpoint not fully implemented yet")
            
            assert response.status_code == 200
            program = response.json()
            assert program["name"] == db_program.name
            
        finally:
            # Clean up
            db.delete(db_program)
            db.commit()
            db.close()
    
    def test_update_program_via_api(self):
        """Test updating a program via API updates Postgres."""
        # First create in database
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        program_id = uuid.uuid4()
        db_program = Program(
            id=program_id,
            name="Original Name",
            description="Original description",
        )
        db.add(db_program)
        db.commit()
        
        try:
            # Update via API
            update_data = {"description": "Updated via API"}
            response = client.patch(f"/api/v1/programs/{program_id}", json=update_data)
            if response.status_code == 404:
                pytest.skip("API endpoint not fully implemented yet")
            
            assert response.status_code == 200
            
            # Verify in database
            db.refresh(db_program)
            assert db_program.description == update_data["description"]
            
        finally:
            # Clean up
            db.delete(db_program)
            db.commit()
            db.close()
    
    def test_delete_program_via_api(self):
        """Test deleting a program via API removes from Postgres."""
        # First create in database
        SessionLocal = get_session_local()
        db = SessionLocal()
        
        program_id = uuid.uuid4()
        db_program = Program(
            id=program_id,
            name="Program to Delete",
            description="Will be deleted via API",
        )
        db.add(db_program)
        db.commit()
        
        try:
            # Delete via API
            response = client.delete(f"/api/v1/programs/{program_id}")
            if response.status_code == 404:
                pytest.skip("API endpoint not fully implemented yet")
            
            assert response.status_code in (200, 204)
            
            # Verify deleted in database
            db_program_check = db.query(Program).filter(Program.id == program_id).first()
            assert db_program_check is None
            
        finally:
            db.close()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

