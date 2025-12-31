"""Tests for inline annotations service layer."""

import uuid

from amprenta_rag.database.session import db_session
from amprenta_rag.services.inline_annotations import (
    create_annotation,
    delete_annotation,
    get_annotation,
    get_annotation_count,
    get_annotations,
    reopen_annotation,
    reply_to_annotation,
    resolve_annotation,
    validate_position_data,
)


def test_validate_position_data_valid():
    """Test position data validation for valid schemas."""
    # Cell position
    assert validate_position_data("cell", {"cell_index": 5}) is True
    
    # Cell range position
    assert validate_position_data("cell_range", {"start_cell": 1, "end_cell": 5}) is True
    assert validate_position_data("range", {"start_cell": 1, "end_cell": 5}) is True  # Alternative name
    
    # Column position
    assert validate_position_data("column", {"column": "gene_name"}) is True
    
    # Row position
    assert validate_position_data("row", {"row_index": 42}) is True
    
    # Field position
    assert validate_position_data("field", {"field": "protocol_version"}) is True


def test_validate_position_data_invalid():
    """Test position data validation for invalid schemas."""
    # Unknown position type
    assert validate_position_data("unknown_type", {"some": "data"}) is False
    
    # Missing required fields
    assert validate_position_data("cell", {}) is False
    assert validate_position_data("cell", {"wrong_field": 5}) is False
    assert validate_position_data("cell_range", {"start_cell": 1}) is False  # Missing end_cell
    assert validate_position_data("column", {"row_index": 42}) is False  # Wrong field
    assert validate_position_data("row", {"column": "gene_name"}) is False  # Wrong field
    assert validate_position_data("field", {"cell_index": 5}) is False  # Wrong field


def test_create_annotation_success():
    """Test successful annotation creation."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        user_id = None  # System annotation
        
        annotation = create_annotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 3},
            content="This cell needs review",
            user_id=user_id,
            db=db,
        )
        
        assert annotation.id is not None
        assert annotation.entity_type == "notebook"
        assert annotation.entity_id == entity_id
        assert annotation.position_type == "cell"
        assert annotation.position_data == {"cell_index": 3}
        assert annotation.content == "This cell needs review"
        assert annotation.created_by_id is None
        assert annotation.status == "open"
        assert annotation.parent_id is None
        assert annotation.created_at is not None


def test_create_annotation_invalid_position():
    """Test annotation creation with invalid position data."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        try:
            create_annotation(
                entity_type="notebook",
                entity_id=entity_id,
                position_type="cell",
                position_data={"invalid_field": 3},  # Wrong field for cell type
                content="This should fail",
                user_id=None,
                db=db,
            )
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "Invalid position_data" in str(e)


def test_get_annotations_all():
    """Test getting all annotations for an entity."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create multiple annotations
        annotation1 = create_annotation(
            entity_type="dataset",
            entity_id=entity_id,
            position_type="column",
            position_data={"column": "gene_name"},
            content="First annotation",
            user_id=None,
            db=db,
        )
        
        annotation2 = create_annotation(
            entity_type="dataset",
            entity_id=entity_id,
            position_type="row",
            position_data={"row_index": 10},
            content="Second annotation",
            user_id=None,
            db=db,
        )
        
        # Get all annotations
        annotations = get_annotations("dataset", entity_id, db)
        
        assert len(annotations) == 2
        annotation_ids = {ann.id for ann in annotations}
        assert annotation1.id in annotation_ids
        assert annotation2.id in annotation_ids


def test_get_annotations_filtered_by_status():
    """Test getting annotations filtered by status."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create annotations with different statuses
        open_annotation = create_annotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 1},
            content="Open annotation",
            user_id=None,
            db=db,
        )
        
        resolved_annotation = create_annotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 2},
            content="Resolved annotation",
            user_id=None,
            db=db,
        )
        
        # Resolve the second annotation
        resolve_annotation(resolved_annotation.id, None, db)
        
        # Get only open annotations
        open_annotations = get_annotations("notebook", entity_id, db, status="open")
        assert len(open_annotations) == 1
        assert open_annotations[0].id == open_annotation.id
        
        # Get only resolved annotations
        resolved_annotations = get_annotations("notebook", entity_id, db, status="resolved")
        assert len(resolved_annotations) == 1
        assert resolved_annotations[0].id == resolved_annotation.id


def test_get_annotations_filtered_by_position_type():
    """Test getting annotations filtered by position type."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create annotations with different position types
        cell_annotation = create_annotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 1},
            content="Cell annotation",
            user_id=None,
            db=db,
        )
        
        range_annotation = create_annotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="range",
            position_data={"start_cell": 1, "end_cell": 5},
            content="Range annotation",
            user_id=None,
            db=db,
        )
        
        # Get only cell annotations
        cell_annotations = get_annotations("notebook", entity_id, db, position_type="cell")
        assert len(cell_annotations) == 1
        assert cell_annotations[0].id == cell_annotation.id
        
        # Get only range annotations
        range_annotations = get_annotations("notebook", entity_id, db, position_type="range")
        assert len(range_annotations) == 1
        assert range_annotations[0].id == range_annotation.id


def test_resolve_annotation():
    """Test resolving an annotation."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create annotation
        annotation = create_annotation(
            entity_type="experiment",
            entity_id=entity_id,
            position_type="field",
            position_data={"field": "protocol_version"},
            content="This needs updating",
            user_id=None,
            db=db,
        )
        
        assert annotation.status == "open"
        assert annotation.resolved_by_id is None
        assert annotation.resolved_at is None
        
        # Resolve annotation
        resolved_annotation = resolve_annotation(annotation.id, None, db)
        
        assert resolved_annotation is not None
        assert resolved_annotation.id == annotation.id
        assert resolved_annotation.status == "resolved"
        assert resolved_annotation.resolved_by_id is None  # System resolution
        assert resolved_annotation.resolved_at is not None


def test_reopen_annotation():
    """Test reopening a resolved annotation."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create and resolve annotation
        annotation = create_annotation(
            entity_type="dataset",
            entity_id=entity_id,
            position_type="column",
            position_data={"column": "expression_level"},
            content="Has outliers",
            user_id=None,
            db=db,
        )
        
        resolve_annotation(annotation.id, None, db)
        
        # Reopen annotation
        reopened_annotation = reopen_annotation(annotation.id, None, db)
        
        assert reopened_annotation is not None
        assert reopened_annotation.id == annotation.id
        assert reopened_annotation.status == "open"
        assert reopened_annotation.resolved_by_id is None
        assert reopened_annotation.resolved_at is None


def test_reply_to_annotation():
    """Test adding a reply to an annotation."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create parent annotation
        parent_annotation = create_annotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 5},
            content="This cell has an error",
            user_id=None,
            db=db,
        )
        
        # Add reply
        reply = reply_to_annotation(
            parent_annotation.id,
            "I've fixed the error",
            None,
            db,
        )
        
        assert reply is not None
        assert reply.parent_id == parent_annotation.id
        assert reply.entity_type == parent_annotation.entity_type
        assert reply.entity_id == parent_annotation.entity_id
        assert reply.position_type == parent_annotation.position_type
        assert reply.position_data == parent_annotation.position_data
        assert reply.content == "I've fixed the error"
        assert reply.created_by_id is None


def test_delete_annotation_author():
    """Test deleting annotation by author."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create annotation
        annotation = create_annotation(
            entity_type="dataset",
            entity_id=entity_id,
            position_type="row",
            position_data={"row_index": 15},
            content="This row is suspicious",
            user_id=None,
            db=db,
        )
        
        # Delete annotation (system can delete any annotation)
        deleted = delete_annotation(annotation.id, None, db)
        assert deleted is True
        
        # Verify annotation is gone
        retrieved = get_annotation(annotation.id, db)
        assert retrieved is None


def test_delete_annotation_not_author():
    """Test deleting annotation by non-author (should fail)."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create annotation by system (None)
        annotation = create_annotation(
            entity_type="notebook",
            entity_id=entity_id,
            position_type="cell",
            position_data={"cell_index": 8},
            content="System annotation",
            user_id=None,  # Created by system
            db=db,
        )
        
        # Try to delete as regular user (should fail - only system can delete system annotations)
        user_id = uuid.uuid4()
        deleted = delete_annotation(annotation.id, user_id, db)
        assert deleted is False
        
        # Verify annotation still exists
        retrieved = get_annotation(annotation.id, db)
        assert retrieved is not None
        assert retrieved.id == annotation.id


def test_get_annotation_count():
    """Test getting annotation counts by status."""
    with db_session() as db:
        entity_id = uuid.uuid4()
        
        # Create annotations with different statuses
        create_annotation(
            entity_type="experiment",
            entity_id=entity_id,
            position_type="field",
            position_data={"field": "temperature"},
            content="Check temperature range",
            user_id=None,
            db=db,
        )
        
        create_annotation(
            entity_type="experiment",
            entity_id=entity_id,
            position_type="field",
            position_data={"field": "duration"},
            content="Duration seems long",
            user_id=None,
            db=db,
        )
        
        resolved_annotation = create_annotation(
            entity_type="experiment",
            entity_id=entity_id,
            position_type="field",
            position_data={"field": "protocol"},
            content="Protocol needs update",
            user_id=None,
            db=db,
        )
        
        # Resolve one annotation
        resolve_annotation(resolved_annotation.id, None, db)
        
        # Get counts
        counts = get_annotation_count("experiment", entity_id, db)
        
        assert counts["open"] == 2
        assert counts["resolved"] == 1
        assert counts["total"] == 3
