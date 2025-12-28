"""Metadata enrichment service for datasets."""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset
from amprenta_rag.database.session import db_session
from amprenta_rag.ingestion.metadata.llm_semantic_extraction import extract_semantic_metadata_with_llm

logger = logging.getLogger(__name__)


@dataclass
class EnrichmentResult:
    """Result of metadata enrichment process."""
    
    dataset_id: UUID
    success: bool
    enriched_fields: List[str]
    extracted_metadata: Dict[str, Any]
    error_message: Optional[str] = None
    processing_time_seconds: Optional[float] = None


def enrich_dataset_metadata(dataset_id: UUID) -> EnrichmentResult:
    """
    Enrich dataset metadata using LLM semantic extraction.
    
    Fetches dataset from database, extracts metadata from description/abstract
    using LLM, merges with existing metadata (LLM takes precedence), and
    updates the dataset record.
    
    Args:
        dataset_id: UUID of the dataset to enrich
        
    Returns:
        EnrichmentResult with success status and enrichment details
    """
    import time
    start_time = time.time()
    
    try:
        with db_session() as db:
            # Fetch dataset
            dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if not dataset:
                return EnrichmentResult(
                    dataset_id=dataset_id,
                    success=False,
                    enriched_fields=[],
                    extracted_metadata={},
                    error_message="Dataset not found",
                    processing_time_seconds=time.time() - start_time,
                )
            
            # Prepare text for extraction
            text_content = _prepare_text_content(dataset)
            if not text_content:
                return EnrichmentResult(
                    dataset_id=dataset_id,
                    success=False,
                    enriched_fields=[],
                    extracted_metadata={},
                    error_message="No text content available for extraction",
                    processing_time_seconds=time.time() - start_time,
                )
            
            logger.info(f"Enriching metadata for dataset {dataset_id} with {len(text_content)} characters of text")
            
            # Extract metadata using LLM
            extracted_metadata = extract_semantic_metadata_with_llm(
                text=text_content,
                source_type="dataset",
                max_length=8000,
            )
            
            # Merge with existing metadata
            enriched_fields = _merge_metadata_with_dataset(db, dataset, extracted_metadata)
            
            # Commit changes
            db.commit()
            
            processing_time = time.time() - start_time
            logger.info(
                f"Successfully enriched dataset {dataset_id} metadata: "
                f"{len(enriched_fields)} fields updated in {processing_time:.2f}s"
            )
            
            return EnrichmentResult(
                dataset_id=dataset_id,
                success=True,
                enriched_fields=enriched_fields,
                extracted_metadata=extracted_metadata,
                processing_time_seconds=processing_time,
            )
            
    except Exception as e:
        processing_time = time.time() - start_time
        error_msg = f"Failed to enrich dataset metadata: {str(e)}"
        logger.error(f"Error enriching dataset {dataset_id}: {e}", exc_info=True)
        
        return EnrichmentResult(
            dataset_id=dataset_id,
            success=False,
            enriched_fields=[],
            extracted_metadata={},
            error_message=error_msg,
            processing_time_seconds=processing_time,
        )


def _prepare_text_content(dataset: Dataset) -> str:
    """Prepare text content from dataset for LLM extraction."""
    text_parts = []
    
    # Add title
    if dataset.title:
        text_parts.append(f"Title: {dataset.title}")
    
    # Add description
    if dataset.description:
        text_parts.append(f"Description: {dataset.description}")
    
    # Add abstract if available (assuming it might be in metadata)
    if hasattr(dataset, 'metadata') and dataset.metadata:
        if isinstance(dataset.metadata, dict):
            abstract = dataset.metadata.get('abstract') or dataset.metadata.get('summary')
            if abstract:
                text_parts.append(f"Abstract: {abstract}")
    
    # Add any existing disease/target info for context
    if hasattr(dataset, 'disease') and dataset.disease:
        if isinstance(dataset.disease, list):
            text_parts.append(f"Known diseases: {', '.join(dataset.disease)}")
        else:
            text_parts.append(f"Known diseases: {dataset.disease}")
    
    return "\n\n".join(text_parts)


def _merge_metadata_with_dataset(db: Session, dataset: Dataset, extracted_metadata: Dict[str, Any]) -> List[str]:
    """
    Merge extracted metadata with dataset record.
    
    LLM-extracted metadata takes precedence over existing values.
    
    Args:
        db: Database session
        dataset: Dataset object to update
        extracted_metadata: Metadata extracted by LLM
        
    Returns:
        List of field names that were updated
    """
    enriched_fields = []
    
    # Initialize metadata dict if it doesn't exist
    if not hasattr(dataset, 'metadata') or dataset.metadata is None:
        dataset.metadata = {}
    elif not isinstance(dataset.metadata, dict):
        dataset.metadata = {}
    
    # Map extracted fields to dataset attributes/metadata
    field_mappings = {
        # Direct attribute mappings (if they exist on the Dataset model)
        "diseases": ("disease", "list"),
        "species": ("organism", "list"),
        "sample_count": ("sample_count", "int"),
        
        # Metadata field mappings
        "targets": ("metadata.targets", "list"),
        "signatures": ("metadata.signatures", "list"),
        "phenotype_axes": ("metadata.phenotype_axes", "list"),
        "biomarker_roles": ("metadata.biomarker_roles", "list"),
        "study_design": ("metadata.study_design", "list"),
        "timepoints": ("metadata.timepoints", "list"),
        "cell_lines": ("metadata.cell_lines", "list"),
        "tissue_types": ("metadata.tissue_types", "list"),
        "assay_types": ("metadata.assay_types", "list"),
        "therapeutic_area": ("metadata.therapeutic_area", "list"),
    }
    
    for extracted_field, value in extracted_metadata.items():
        if not value:  # Skip empty values
            continue
            
        if extracted_field in field_mappings:
            field_path, value_type = field_mappings[extracted_field]
            
            if field_path.startswith("metadata."):
                # Update metadata field
                metadata_key = field_path.split(".", 1)[1]
                old_value = dataset.metadata.get(metadata_key)
                
                if value_type == "list":
                    # Merge lists, LLM values take precedence
                    existing_values = set(old_value) if isinstance(old_value, list) else set()
                    new_values = set(value) if isinstance(value, list) else {value}
                    merged_values = sorted(list(existing_values | new_values))
                    
                    if merged_values != old_value:
                        dataset.metadata[metadata_key] = merged_values
                        enriched_fields.append(metadata_key)
                        logger.debug(f"Updated metadata.{metadata_key}: {old_value} -> {merged_values}")
                
                elif value_type == "int":
                    if value != old_value:
                        dataset.metadata[metadata_key] = value
                        enriched_fields.append(metadata_key)
                        logger.debug(f"Updated metadata.{metadata_key}: {old_value} -> {value}")
            
            else:
                # Update direct attribute
                attr_name = field_path
                if hasattr(dataset, attr_name):
                    old_value = getattr(dataset, attr_name)
                    
                    if value_type == "list":
                        # Merge lists
                        existing_values = set(old_value) if isinstance(old_value, list) else set()
                        new_values = set(value) if isinstance(value, list) else {value}
                        merged_values = sorted(list(existing_values | new_values))
                        
                        if merged_values != old_value:
                            setattr(dataset, attr_name, merged_values)
                            enriched_fields.append(attr_name)
                            logger.debug(f"Updated {attr_name}: {old_value} -> {merged_values}")
                    
                    elif value_type == "int":
                        if value != old_value:
                            setattr(dataset, attr_name, value)
                            enriched_fields.append(attr_name)
                            logger.debug(f"Updated {attr_name}: {old_value} -> {value}")
    
    # Mark the dataset as having been enriched
    dataset.metadata["enrichment_timestamp"] = time.time()
    dataset.metadata["enrichment_version"] = "1.0"
    
    return enriched_fields


def batch_enrich_datasets(dataset_ids: List[UUID], max_concurrent: int = 5) -> List[EnrichmentResult]:
    """
    Enrich multiple datasets in batch.
    
    Args:
        dataset_ids: List of dataset UUIDs to enrich
        max_concurrent: Maximum number of concurrent enrichment operations
        
    Returns:
        List of EnrichmentResult objects
    """
    results = []
    
    # For now, process sequentially to avoid overwhelming the LLM API
    # In the future, this could be made concurrent with proper rate limiting
    for dataset_id in dataset_ids:
        try:
            result = enrich_dataset_metadata(dataset_id)
            results.append(result)
            
            # Brief pause between requests to respect API rate limits
            import time
            time.sleep(1)
            
        except Exception as e:
            logger.error(f"Error in batch enrichment for dataset {dataset_id}: {e}")
            results.append(EnrichmentResult(
                dataset_id=dataset_id,
                success=False,
                enriched_fields=[],
                extracted_metadata={},
                error_message=f"Batch processing error: {str(e)}",
            ))
    
    logger.info(f"Batch enrichment completed: {len(results)} datasets processed")
    return results


def get_enrichment_status(dataset_id: UUID) -> Dict[str, Any]:
    """
    Get enrichment status for a dataset.
    
    Args:
        dataset_id: UUID of the dataset
        
    Returns:
        Dictionary with enrichment status information
    """
    try:
        with db_session() as db:
            dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if not dataset:
                return {"error": "Dataset not found"}
            
            if not hasattr(dataset, 'metadata') or not dataset.metadata:
                return {
                    "enriched": False,
                    "enrichment_timestamp": None,
                    "enrichment_version": None,
                }
            
            metadata = dataset.metadata if isinstance(dataset.metadata, dict) else {}
            
            return {
                "enriched": "enrichment_timestamp" in metadata,
                "enrichment_timestamp": metadata.get("enrichment_timestamp"),
                "enrichment_version": metadata.get("enrichment_version"),
                "available_fields": list(metadata.keys()),
            }
            
    except Exception as e:
        logger.error(f"Error getting enrichment status for dataset {dataset_id}: {e}")
        return {"error": str(e)}


__all__ = [
    "EnrichmentResult",
    "enrich_dataset_metadata",
    "batch_enrich_datasets",
    "get_enrichment_status",
]
