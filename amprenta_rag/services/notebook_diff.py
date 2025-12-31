"""Service layer for notebook diff computation and snapshot management."""

from __future__ import annotations

import difflib
import json
from pathlib import Path
from typing import Dict, List, Optional, Any
from uuid import UUID

import nbformat
from sqlalchemy.orm import Session

from amprenta_rag.database.models import NotebookSnapshot
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.notebooks.registry import load_registry, resolve_notebook_path

logger = get_logger(__name__)


def capture_snapshot(db: Session, review_id: UUID, notebook_path: str) -> NotebookSnapshot:
    """
    Capture a snapshot of a notebook at review time.
    
    Args:
        db: Database session
        review_id: ID of the notebook review
        notebook_path: Path to the notebook file
    
    Returns:
        Created NotebookSnapshot instance
    
    Raises:
        FileNotFoundError: If notebook file not found
        ValueError: If notebook cannot be parsed
    """
    logger.info(f"Capturing snapshot for review {review_id}, notebook: {notebook_path}")
    
    try:
        # Try to resolve from registry first
        registry = load_registry()
        notebook_file_path = None
        
        for template in registry:
            if template.get("notebook_path") == notebook_path:
                notebook_file_path = resolve_notebook_path(template)
                break
        
        # If not in registry, treat as direct file path
        if not notebook_file_path:
            notebook_file_path = Path(notebook_path)
        
        if not notebook_file_path.exists():
            raise FileNotFoundError(f"Notebook not found: {notebook_path}")
        
        # Load and parse notebook
        with open(notebook_file_path, 'r', encoding='utf-8') as f:
            notebook_content = nbformat.read(f, as_version=nbformat.NO_CONVERT)
        
        # Convert to dict for JSON storage
        content_dict = nbformat.v4.writes(notebook_content)
        content_json = json.loads(content_dict)
        
        # Count cells
        cell_count = len(content_json.get('cells', []))
        
        # Create snapshot
        snapshot = NotebookSnapshot(
            review_id=review_id,
            content_json=content_json,
            cell_count=cell_count
        )
        
        db.add(snapshot)
        db.commit()
        db.refresh(snapshot)
        
        logger.info(f"Created snapshot {snapshot.id} for review {review_id} ({cell_count} cells)")
        return snapshot
        
    except Exception as e:
        logger.error(f"Failed to capture snapshot for review {review_id}: {e}")
        raise


def get_snapshot(db: Session, review_id: UUID) -> Optional[NotebookSnapshot]:
    """
    Get the snapshot for a review.
    
    Args:
        db: Database session
        review_id: ID of the notebook review
    
    Returns:
        NotebookSnapshot instance or None if not found
    """
    logger.debug(f"Fetching snapshot for review {review_id}")
    
    snapshot = db.query(NotebookSnapshot).filter(
        NotebookSnapshot.review_id == review_id
    ).first()
    
    if snapshot:
        logger.debug(f"Found snapshot {snapshot.id} for review {review_id}")
    else:
        logger.debug(f"No snapshot found for review {review_id}")
    
    return snapshot


def compute_cell_diff(old_cells: List[Dict], new_cells: List[Dict]) -> Dict[str, List[Dict]]:
    """
    Compute differences between two sets of notebook cells.
    
    Args:
        old_cells: Original cells from snapshot
        new_cells: Current cells from notebook
    
    Returns:
        Dict with keys: added, removed, modified, unchanged
        Each containing list of cell info dicts with: index, cell_type, source_preview
    """
    logger.debug(f"Computing diff between {len(old_cells)} old and {len(new_cells)} new cells")
    
    result = {
        "added": [],
        "removed": [],
        "modified": [],
        "unchanged": []
    }
    
    # Process cells for comparison
    
    # Track which cells we've processed
    processed_old = set()
    processed_new = set()
    
    # First pass: find exact matches and modifications
    for new_idx, new_cell in enumerate(new_cells):
        best_match_idx = None
        best_match_score = 0
        
        for old_idx, old_cell in enumerate(old_cells):
            if old_idx in processed_old:
                continue
                
            # Check if cells are similar enough to be considered the same
            old_source = old_cell.get('source', '')
            new_source = new_cell.get('source', '')
            
            if old_cell.get('cell_type') == new_cell.get('cell_type'):
                if old_source == new_source:
                    # Exact match
                    best_match_idx = old_idx
                    best_match_score = 1.0
                    break
                else:
                    # Check similarity
                    similarity = difflib.SequenceMatcher(None, old_source, new_source).ratio()
                    if similarity > best_match_score and similarity > 0.6:
                        best_match_idx = old_idx
                        best_match_score = similarity
        
        if best_match_idx is not None:
            old_cell = old_cells[best_match_idx]
            processed_old.add(best_match_idx)
            processed_new.add(new_idx)
            
            cell_info = {
                "index": new_idx,
                "cell_type": new_cell.get('cell_type', 'unknown'),
                "source_preview": _get_source_preview(new_cell.get('source', ''))
            }
            
            if best_match_score == 1.0:
                result["unchanged"].append(cell_info)
            else:
                result["modified"].append({
                    **cell_info,
                    "old_index": best_match_idx,
                    "similarity": best_match_score
                })
    
    # Second pass: identify added cells
    for new_idx, new_cell in enumerate(new_cells):
        if new_idx not in processed_new:
            result["added"].append({
                "index": new_idx,
                "cell_type": new_cell.get('cell_type', 'unknown'),
                "source_preview": _get_source_preview(new_cell.get('source', ''))
            })
    
    # Third pass: identify removed cells
    for old_idx, old_cell in enumerate(old_cells):
        if old_idx not in processed_old:
            result["removed"].append({
                "index": old_idx,
                "cell_type": old_cell.get('cell_type', 'unknown'),
                "source_preview": _get_source_preview(old_cell.get('source', ''))
            })
    
    logger.debug(f"Diff computed: {len(result['added'])} added, {len(result['removed'])} removed, "
                f"{len(result['modified'])} modified, {len(result['unchanged'])} unchanged")
    
    return result


def get_review_diff(db: Session, review_id: UUID, current_notebook_path: str) -> Dict[str, Any]:
    """
    Get diff between stored snapshot and current notebook content.
    
    Args:
        db: Database session
        review_id: ID of the notebook review
        current_notebook_path: Path to current notebook file
    
    Returns:
        Dict with diff summary and cell-level changes
    
    Raises:
        ValueError: If no snapshot found or notebook cannot be loaded
    """
    logger.info(f"Computing review diff for review {review_id}")
    
    # Get stored snapshot
    snapshot = get_snapshot(db, review_id)
    if not snapshot:
        raise ValueError(f"No snapshot found for review {review_id}")
    
    # Load current notebook
    try:
        # Try to resolve from registry first
        registry = load_registry()
        notebook_file_path = None
        
        for template in registry:
            if template.get("notebook_path") == current_notebook_path:
                notebook_file_path = resolve_notebook_path(template)
                break
        
        # If not in registry, treat as direct file path
        if not notebook_file_path:
            notebook_file_path = Path(current_notebook_path)
        
        if not notebook_file_path.exists():
            raise FileNotFoundError(f"Current notebook not found: {current_notebook_path}")
        
        with open(notebook_file_path, 'r', encoding='utf-8') as f:
            current_notebook = nbformat.read(f, as_version=nbformat.NO_CONVERT)
        
        current_content = json.loads(nbformat.v4.writes(current_notebook))
        
    except Exception as e:
        logger.error(f"Failed to load current notebook {current_notebook_path}: {e}")
        raise ValueError(f"Cannot load current notebook: {e}")
    
    # Extract cells
    old_cells = snapshot.content_json.get('cells', [])
    new_cells = current_content.get('cells', [])
    
    # Compute diff
    cell_diff = compute_cell_diff(old_cells, new_cells)
    
    # Create summary
    summary = {
        "snapshot_id": str(snapshot.id),
        "snapshot_created_at": snapshot.created_at.isoformat(),
        "old_cell_count": len(old_cells),
        "new_cell_count": len(new_cells),
        "changes": {
            "added": len(cell_diff["added"]),
            "removed": len(cell_diff["removed"]),
            "modified": len(cell_diff["modified"]),
            "unchanged": len(cell_diff["unchanged"])
        },
        "cell_diff": cell_diff
    }
    
    logger.info(f"Review diff computed for review {review_id}: "
               f"{summary['changes']['added']} added, {summary['changes']['removed']} removed, "
               f"{summary['changes']['modified']} modified cells")
    
    return summary


def _get_source_preview(source: str, max_length: int = 100) -> str:
    """
    Get a preview of cell source code.
    
    Args:
        source: Cell source content
        max_length: Maximum length of preview
    
    Returns:
        Truncated source preview
    """
    if not source:
        return ""
    
    # Handle list of strings (common in notebooks)
    if isinstance(source, list):
        source = "".join(source)
    
    # Clean up and truncate
    preview = source.strip().replace('\n', ' ').replace('\r', '')
    if len(preview) > max_length:
        preview = preview[:max_length-3] + "..."
    
    return preview
