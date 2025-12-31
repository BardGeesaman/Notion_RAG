"""Tests for InlineAnnotation model."""

from __future__ import annotations

import uuid
from datetime import datetime, timezone

from sqlalchemy import inspect

from amprenta_rag.database.base import get_engine
from amprenta_rag.database.session import db_session
from amprenta_rag.models.user_prefs import InlineAnnotation


def test_inline_annotation_create():
    """Test creating a basic inline annotation."""
    with db_session() as db:
        annotation = InlineAnnotation(
            entity_type="notebook",
            entity_id=uuid.uuid4(),
            position_type="cell",
            position_data={"cell_index": 5},
            content="This cell needs clarification",
            created_by_id=None
        )
        
        db.add(annotation)
        db.commit()
        db.refresh(annotation)
        
        assert annotation.id is not None
        assert annotation.entity_type == "notebook"
        assert annotation.position_type == "cell"
        assert annotation.position_data == {"cell_index": 5}
        assert annotation.content == "This cell needs clarification"
        assert annotation.status == "open"  # Default status
        assert annotation.created_at is not None
        assert annotation.parent_id is None
        assert annotation.resolved_by_id is None
        assert annotation.resolved_at is None


def test_inline_annotation_position_types():
    """Test different position types and their data structures."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        created_by_id = None
        
        # Test cell position
        cell_annotation = InlineAnnotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 3},
            content="Cell comment",
            created_by_id=created_by_id
        )
        
        # Test column position
        column_annotation = InlineAnnotation(
            entity_type="dataset",
            entity_id=entity_id,
            position_type="column",
            position_data={"column": "gene_name"},
            content="Column comment",
            created_by_id=created_by_id
        )
        
        # Test row position
        row_annotation = InlineAnnotation(
            entity_type="dataset",
            entity_id=entity_id,
            position_type="row",
            position_data={"row_index": 42},
            content="Row comment",
            created_by_id=created_by_id
        )
        
        # Test range position
        range_annotation = InlineAnnotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="range",
            position_data={"start_cell": 5, "end_cell": 8},
            content="Range comment",
            created_by_id=created_by_id
        )
        
        db.add_all([cell_annotation, column_annotation, row_annotation, range_annotation])
        db.commit()
        
        # Verify all were created with correct position data
        annotations = db.query(InlineAnnotation).filter(InlineAnnotation.entity_id == entity_id).all()
        assert len(annotations) == 4
        
        position_types = {ann.position_type: ann.position_data for ann in annotations}
        assert position_types["cell"] == {"cell_index": 3}
        assert position_types["column"] == {"column": "gene_name"}
        assert position_types["row"] == {"row_index": 42}
        assert position_types["range"] == {"start_cell": 5, "end_cell": 8}


def test_inline_annotation_status_default():
    """Test that default status is 'open' and can be changed to 'resolved'."""
    with db_session() as db:
        annotation = InlineAnnotation(
            entity_type="experiment",
            entity_id=uuid.uuid4(),
            position_type="field",
            position_data={"field": "protocol_version"},
            content="This field needs updating",
            created_by_id=None
        )
        
        db.add(annotation)
        db.commit()
        db.refresh(annotation)
        
        # Check default status
        assert annotation.status == "open"
        
        # Update status to resolved
        annotation.status = "resolved"
        annotation.resolved_by_id = None  # Use None since no real user
        annotation.resolved_at = datetime.now(timezone.utc)
        
        db.commit()
        db.refresh(annotation)
        
        assert annotation.status == "resolved"
        assert annotation.resolved_by_id is None
        assert annotation.resolved_at is not None


def test_inline_annotation_threading():
    """Test parent-child relationships for threaded replies."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        user1_id = None
        user2_id = None
        
        # Create parent annotation
        parent_annotation = InlineAnnotation(
            entity_type="dataset",
            entity_id=entity_id,
            position_type="column",
            position_data={"column": "expression_level"},
            content="This column has outliers",
            created_by_id=user1_id
        )
        
        db.add(parent_annotation)
        db.commit()
        db.refresh(parent_annotation)
        
        # Create reply annotation
        reply_annotation = InlineAnnotation(
            entity_type="dataset",
            entity_id=entity_id,
            position_type="column",
            position_data={"column": "expression_level"},
            content="I've verified the outliers are real data points",
            created_by_id=user2_id,
            parent_id=parent_annotation.id
        )
        
        db.add(reply_annotation)
        db.commit()
        db.refresh(reply_annotation)
        
        # Test relationships
        assert reply_annotation.parent_id == parent_annotation.id
        assert reply_annotation.parent is not None
        assert reply_annotation.parent.id == parent_annotation.id
        
        # Test that parent has replies
        db.refresh(parent_annotation)
        assert len(parent_annotation.replies) == 1
        assert parent_annotation.replies[0].id == reply_annotation.id
        assert parent_annotation.replies[0].content == "I've verified the outliers are real data points"


def test_inline_annotation_indexes_exist():
    """Test that the expected database indexes exist."""
    engine = get_engine()
    inspector = inspect(engine)
    
    # Get all indexes for the inline_annotations table
    indexes = inspector.get_indexes("inline_annotations")
    index_names = {idx["name"] for idx in indexes}
    
    # Check that our custom indexes exist
    expected_indexes = {
        "ix_inline_annotations_entity",
        "ix_inline_annotations_status", 
        "ix_inline_annotations_position_type",
        "ix_inline_annotations_entity_type",  # Auto-generated from index=True
        "ix_inline_annotations_entity_id"    # Auto-generated from index=True
    }
    
    for expected_index in expected_indexes:
        assert expected_index in index_names, f"Index {expected_index} not found. Available: {index_names}"


def test_inline_annotation_unique_constraint_none():
    """Test that multiple annotations can exist at the same position (no unique constraint)."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        user1_id = None
        user2_id = None
        
        # Create two annotations at the same position
        annotation1 = InlineAnnotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 10},
            content="First comment on this cell",
            created_by_id=user1_id
        )
        
        annotation2 = InlineAnnotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 10},
            content="Second comment on the same cell",
            created_by_id=user2_id
        )
        
        db.add_all([annotation1, annotation2])
        db.commit()  # Should succeed - no unique constraint violation
        
        # Verify both annotations exist
        annotations = db.query(InlineAnnotation).filter(
            InlineAnnotation.entity_id == entity_id,
            InlineAnnotation.position_type == "cell"
        ).all()
        
        assert len(annotations) == 2
        contents = {ann.content for ann in annotations}
        assert "First comment on this cell" in contents
        assert "Second comment on the same cell" in contents
