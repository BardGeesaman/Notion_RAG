"""
Publication data extractor for scientific papers.

Extracts structured experiment data from publication PDFs using:
- Section detection (Abstract, Methods, Results, Discussion)
- LLM-based structured extraction with Pydantic validation
- Token limits and caching for efficiency
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from typing import Any, Dict, List, Optional
from uuid import UUID

from pydantic import BaseModel, Field

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.ingestion.text_extraction import extract_text_from_pdf_bytes
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Simple in-memory cache for extraction results
_EXTRACTION_CACHE: Dict[str, Any] = {}


class ExtractedExperiment(BaseModel):
    """Structured experiment data extracted from publication."""

    experiment_type: str = Field(..., description="Type of experiment: RNA-seq, Western blot, ELISA, qPCR, etc.")
    cell_line: Optional[str] = Field(None, description="Cell line or tissue type used")
    treatment: Optional[str] = Field(None, description="Treatment or intervention applied")
    concentration: Optional[str] = Field(None, description="Concentration of treatment")
    timepoint: Optional[str] = Field(None, description="Time point of measurement")
    replicate_count: Optional[int] = Field(None, description="Number of biological replicates")
    measured_entities: List[str] = Field(default_factory=list, description="Genes, proteins, metabolites measured")
    key_findings: Optional[str] = Field(None, description="Brief summary of key findings")


class PublicationExtractionResult(BaseModel):
    """Result of publication data extraction."""

    experiments: List[ExtractedExperiment] = Field(default_factory=list)
    methods_text: Optional[str] = Field(None, description="Extracted methods section text")
    confidence: float = Field(default=0.0, ge=0.0, le=1.0, description="Confidence score for extraction")
    cached: bool = Field(default=False, description="Whether result was retrieved from cache")


def detect_sections(text: str) -> Dict[str, str]:
    """
    Detect major sections in publication text.

    Uses regex patterns to identify section headers and extract content.

    Args:
        text: Full publication text

    Returns:
        Dictionary mapping section names to section content
    """
    sections = {}
    
    # Common section header patterns
    section_patterns = [
        (r"(?i)^abstract\s*$", "Abstract"),
        (r"(?i)^introduction\s*$", "Introduction"),
        (r"(?i)^methods?\s*$", "Methods"),
        (r"(?i)^materials?\s+and\s+methods?\s*$", "Methods"),
        (r"(?i)^experimental\s+procedures?\s*$", "Methods"),
        (r"(?i)^results?\s*$", "Results"),
        (r"(?i)^discussion\s*$", "Discussion"),
        (r"(?i)^conclusion\s*$", "Conclusion"),
        (r"(?i)^references?\s*$", "References"),
    ]
    
    # Split by lines and find section boundaries
    lines = text.split("\n")
    current_section = None
    current_content = []
    
    for line in lines:
        line_stripped = line.strip()
        
        # Check if this line is a section header
        matched_section = None
        for pattern, section_name in section_patterns:
            if re.match(pattern, line_stripped):
                matched_section = section_name
                break
        
        if matched_section:
            # Save previous section
            if current_section and current_content:
                sections[current_section] = "\n".join(current_content).strip()
            
            # Start new section
            current_section = matched_section
            current_content = []
        elif current_section:
            # Add to current section
            current_content.append(line)
    
    # Save last section
    if current_section and current_content:
        sections[current_section] = "\n".join(current_content).strip()
    
    logger.debug("[PUBLICATION_EXTRACTOR] Detected %d sections", len(sections))
    return sections


def extract_methods_section(text: str) -> Optional[str]:
    """
    Extract the Methods/Materials and Methods section from publication text.

    Args:
        text: Full publication text

    Returns:
        Methods section text, or None if not found
    """
    sections = detect_sections(text)
    
    methods_text = sections.get("Methods")
    
    if methods_text:
        logger.info("[PUBLICATION_EXTRACTOR] Extracted Methods section (%d chars)", len(methods_text))
    else:
        logger.warning("[PUBLICATION_EXTRACTOR] No Methods section found")
    
    return methods_text


def extract_experiments_from_text(text: str, max_tokens: int = 8000) -> PublicationExtractionResult:
    """
    Extract structured experiment data from publication text using LLM.

    Uses OpenAI API with structured output to extract experiment details.
    Implements token limiting and caching for efficiency.

    Args:
        text: Publication text (full text or methods section)
        max_tokens: Maximum tokens to send to LLM (default 8000)

    Returns:
        PublicationExtractionResult with extracted experiments
    """
    # Truncate text to approximate token limit (rough: 1 token ≈ 4 chars)
    max_chars = max_tokens * 4
    if len(text) > max_chars:
        text = text[:max_chars]
        logger.info("[PUBLICATION_EXTRACTOR] Truncated text to %d chars (~%d tokens)", max_chars, max_tokens)
    
    # Check cache
    cache_key = hashlib.md5(text.encode()).hexdigest()
    if cache_key in _EXTRACTION_CACHE:
        logger.info("[PUBLICATION_EXTRACTOR] Cache hit for text")
        result = _EXTRACTION_CACHE[cache_key]
        result.cached = True
        return result
    
    # Get LLM model
    cfg_model, _ = get_default_models()
    model = os.getenv("AMPRENTA_PUBLICATION_EXTRACTOR_MODEL", cfg_model or "gpt-4o-mini")
    
    # Build prompt
    schema = {
        "type": "object",
        "properties": {
            "experiments": {
                "type": "array",
                "items": ExtractedExperiment.model_json_schema()
            }
        }
    }
    
    system_prompt = (
        "You are a scientific data extractor specialized in biomedical research.\n"
        "Extract structured experiment data from the Methods section of scientific publications.\n"
        "Return ONLY valid JSON that matches the provided schema.\n"
        "Be conservative: only extract experiments explicitly described.\n"
    )
    
    user_prompt = (
        f"Extract all experiments from this publication text:\n\n{text}\n\n"
        f"Return JSON matching this schema:\n{json.dumps(schema, indent=2)}"
    )
    
    try:
        client = get_openai_client()
        
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0.0,
            max_tokens=4000,
        )
        
        raw_response = (response.choices[0].message.content or "").strip()
        
        # Parse JSON response
        data = json.loads(raw_response)
        
        # Validate and create result
        experiments = [ExtractedExperiment.model_validate(exp) for exp in data.get("experiments", [])]
        
        result = PublicationExtractionResult(
            experiments=experiments,
            methods_text=text[:1000] if text else None,  # Store truncated methods
            confidence=0.8 if experiments else 0.3,  # Simple confidence heuristic
            cached=False,
        )
        
        # Cache result
        _EXTRACTION_CACHE[cache_key] = result
        
        logger.info("[PUBLICATION_EXTRACTOR] Extracted %d experiments", len(experiments))
        return result
    
    except json.JSONDecodeError as e:
        logger.error("[PUBLICATION_EXTRACTOR] Failed to parse LLM response as JSON: %r", e)
        return PublicationExtractionResult(experiments=[], confidence=0.0)
    
    except Exception as e:
        logger.error("[PUBLICATION_EXTRACTOR] Extraction failed: %r", e)
        return PublicationExtractionResult(experiments=[], confidence=0.0)


def extract_from_pdf_bytes(pdf_bytes: bytes, literature_id: Optional[UUID] = None) -> PublicationExtractionResult:
    """
    Extract structured experiment data from PDF bytes.

    Full pipeline: PDF → text → section detection → LLM extraction.

    Args:
        pdf_bytes: PDF file content as bytes
        literature_id: Optional Literature record ID for tracking

    Returns:
        PublicationExtractionResult with extracted experiments
    """
    try:
        # Extract text from PDF
        text = extract_text_from_pdf_bytes(pdf_bytes)
        
        if not text or len(text) < 100:
            logger.warning("[PUBLICATION_EXTRACTOR] PDF text too short (%d chars)", len(text))
            return PublicationExtractionResult(experiments=[], confidence=0.0)
        
        logger.info(
            "[PUBLICATION_EXTRACTOR] Extracted %d chars from PDF (literature_id=%s)",
            len(text),
            literature_id,
        )
        
        # Try to extract Methods section specifically
        methods_text = extract_methods_section(text)
        
        # Use Methods section if available, otherwise use full text (truncated)
        extraction_text = methods_text if methods_text else text
        
        # Extract experiments
        result = extract_experiments_from_text(extraction_text, max_tokens=8000)
        
        return result
    
    except Exception as e:
        logger.error("[PUBLICATION_EXTRACTOR] PDF extraction failed: %r", e)
        return PublicationExtractionResult(experiments=[], confidence=0.0)


def clear_cache() -> None:
    """Clear the extraction cache."""
    _EXTRACTION_CACHE.clear()
    logger.info("[PUBLICATION_EXTRACTOR] Cache cleared")

