"""
LLM-based semantic metadata extraction.

Uses OpenAI API to extract structured metadata (diseases, targets, signatures, etc.)
from unstructured text content. Provides more accurate extraction than pattern matching.
"""

from __future__ import annotations

import json
from typing import Any, Dict

from amprenta_rag.clients.openai_client import get_openai_client
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


EXTRACTION_PROMPT_TEMPLATE = """Extract structured metadata from the following text content.

Text content:
{text_content}

Extract the following information:
1. **Diseases**: List any disease names mentioned (e.g., "Alzheimer's disease", "ALS", "Fragile X syndrome")
2. **Targets**: List any target molecules, proteins, or genes mentioned (e.g., "TP53", "TNF-alpha", "CDK4")
3. **Signatures**: List any biomarker signatures, lipid signatures, or signature identifiers mentioned
4. **Phenotype Axes**: List any phenotype axes, disease axes, or clinical axes mentioned
5. **Biomarker Roles**: List any biomarker roles mentioned (e.g., "diagnostic", "prognostic", "predictive")
6. **Sample Count**: Extract the number of samples/subjects if mentioned (return as integer, null if not found)
7. **Study Design**: List study design types mentioned (e.g., "case-control", "longitudinal", "cross-sectional", "randomized")
8. **Timepoints**: List any timepoints mentioned (e.g., "baseline", "6 months", "1 year", "pre-treatment", "post-treatment")
9. **Cell Lines**: List any cell line names mentioned (e.g., "HeLa", "MCF-7", "HEK293")
10. **Tissue Types**: List any tissue/organ types mentioned (e.g., "brain", "liver", "blood", "muscle", "adipose")
11. **Species**: List any species/organisms mentioned (e.g., "human", "mouse", "rat", "zebrafish")
12. **Assay Types**: List any experimental assay types mentioned (e.g., "RNA-seq", "proteomics", "metabolomics", "microarray", "ChIP-seq")
13. **Therapeutic Area**: List any therapeutic areas mentioned (e.g., "oncology", "neurology", "cardiology", "immunology")

Return your response as a JSON object with the following structure:
{{
  "diseases": ["disease1", "disease2"],
  "targets": ["target1", "target2"],
  "signatures": ["signature1", "signature2"],
  "phenotype_axes": ["axis1", "axis2"],
  "biomarker_roles": ["role1", "role2"],
  "sample_count": 123,
  "study_design": ["design1", "design2"],
  "timepoints": ["timepoint1", "timepoint2"],
  "cell_lines": ["cell_line1", "cell_line2"],
  "tissue_types": ["tissue1", "tissue2"],
  "species": ["species1", "species2"],
  "assay_types": ["assay1", "assay2"],
  "therapeutic_area": ["area1", "area2"]
}}

If a category has no matches, return an empty array for that field. For sample_count, return null if no count is mentioned.
Only extract information that is explicitly mentioned in the text. Do not infer or make assumptions.
"""


def extract_semantic_metadata_with_llm(
    text: str,
    source_type: str = "generic",
    max_length: int = 8000,
) -> Dict[str, Any]:
    """
    Extract semantic metadata from text using OpenAI LLM.

    Uses GPT to extract structured metadata (diseases, targets, signatures, etc.)
    from unstructured text content. More accurate than pattern matching.

    Args:
        text: Text content to extract metadata from
        source_type: Type of source (dataset, experiment, email, literature) - for logging
        max_length: Maximum text length to send to LLM (will truncate if longer)

    Returns:
        Dictionary with extracted metadata:
        - diseases: List[str]
        - targets: List[str]
        - signatures: List[str]
        - phenotype_axes: List[str]
        - biomarker_roles: List[str]
        - sample_count: int | None
        - study_design: List[str]
        - timepoints: List[str]
        - cell_lines: List[str]
        - tissue_types: List[str]
        - species: List[str]
        - assay_types: List[str]
        - therapeutic_area: List[str]
    """
    if not text or len(text.strip()) < 50:
        logger.debug("[LLM-SEMANTIC] Text too short for LLM extraction")
        return _empty_metadata()

    cfg = get_config()

    # Check if LLM extraction is enabled
    if not cfg.pipeline.enable_llm_semantic_extraction:
        logger.debug("[LLM-SEMANTIC] LLM semantic extraction is disabled")
        return _empty_metadata()

    # Truncate text if too long
    if len(text) > max_length:
        text = text[:max_length] + "... [truncated]"
        logger.debug(
            "[LLM-SEMANTIC] Truncated text to %d characters for LLM extraction",
            max_length,
        )

    try:
        client = get_openai_client()

        prompt = EXTRACTION_PROMPT_TEMPLATE.format(text_content=text)

        logger.debug(
            "[LLM-SEMANTIC] Calling OpenAI API for semantic extraction (%s, %d chars)",
            source_type,
            len(text),
        )

        response = client.chat.completions.create(
            model=cfg.openai.model,
            messages=[
                {
                    "role": "system",
                    "content": "You are a metadata extraction assistant. Extract structured information from text and return it as JSON.",
                },
                {"role": "user", "content": prompt},
            ],
            temperature=0.0,  # Deterministic extraction
            max_tokens=500,
            response_format={"type": "json_object"},
        )

        response_text = response.choices[0].message.content
        if not response_text:
            logger.warning("[LLM-SEMANTIC] Empty response from OpenAI")
            return _empty_metadata()

        # Parse JSON response
        try:
            extracted = json.loads(response_text)

            # Validate and normalize structure
            result = {
                "diseases": _normalize_list(extracted.get("diseases", [])),
                "targets": _normalize_list(extracted.get("targets", [])),
                "signatures": _normalize_list(extracted.get("signatures", [])),
                "phenotype_axes": _normalize_list(extracted.get("phenotype_axes", [])),
                "biomarker_roles": _normalize_list(extracted.get("biomarker_roles", [])),
                "sample_count": _normalize_int(extracted.get("sample_count")),
                "study_design": _normalize_list(extracted.get("study_design", [])),
                "timepoints": _normalize_list(extracted.get("timepoints", [])),
                "cell_lines": _normalize_list(extracted.get("cell_lines", [])),
                "tissue_types": _normalize_list(extracted.get("tissue_types", [])),
                "species": _normalize_list(extracted.get("species", [])),
                "assay_types": _normalize_list(extracted.get("assay_types", [])),
                "therapeutic_area": _normalize_list(extracted.get("therapeutic_area", [])),
            }

            logger.info(
                "[LLM-SEMANTIC] Extracted metadata: %d diseases, %d targets, %d signatures, %d axes, %d roles, %s samples, %d designs, %d timepoints, %d cell_lines, %d tissues, %d species, %d assays, %d therapeutic_areas",
                len(result["diseases"]),
                len(result["targets"]),
                len(result["signatures"]),
                len(result["phenotype_axes"]),
                len(result["biomarker_roles"]),
                result["sample_count"] or "no",
                len(result["study_design"]),
                len(result["timepoints"]),
                len(result["cell_lines"]),
                len(result["tissue_types"]),
                len(result["species"]),
                len(result["assay_types"]),
                len(result["therapeutic_area"]),
            )

            return result

        except json.JSONDecodeError as e:
            logger.error(
                "[LLM-SEMANTIC] Failed to parse JSON response from OpenAI: %r. Response: %s",
                e,
                response_text[:200],
            )
            return _empty_metadata()

    except Exception as e:
        logger.warning(
            "[LLM-SEMANTIC] Error calling OpenAI API for semantic extraction: %r",
            e,
        )
        return _empty_metadata()


def enhance_metadata_with_llm(
    metadata: Dict[str, Any],
    text: str,
    source_type: str = "generic",
) -> Dict[str, Any]:
    """
    Enhance existing metadata by extracting additional information with LLM.

    Extracts semantic metadata from text and merges it with existing metadata,
    avoiding duplicates.

    Args:
        metadata: Existing metadata dictionary
        text: Text content to extract from
        source_type: Type of source (for logging)

    Returns:
        Enhanced metadata dictionary with merged results
    """
    extracted = extract_semantic_metadata_with_llm(text, source_type)

    # Merge with existing metadata
    result = metadata.copy()

    # Merge diseases
    existing_diseases = set(result.get("diseases", []))
    extracted_diseases = set(extracted.get("diseases", []))
    result["diseases"] = sorted(list(existing_diseases | extracted_diseases))

    # Merge targets
    existing_targets = set(result.get("targets", []))
    extracted_targets = set(extracted.get("targets", []))
    result["targets"] = sorted(list(existing_targets | extracted_targets))

    # Merge signatures
    existing_sigs = set(result.get("lipid_signatures", []))
    extracted_sigs = set(extracted.get("signatures", []))
    result["lipid_signatures"] = sorted(list(existing_sigs | extracted_sigs))

    # Merge phenotype axes
    existing_axes = set(result.get("phenotype_axes", []))
    extracted_axes = set(extracted.get("phenotype_axes", []))
    result["phenotype_axes"] = sorted(list(existing_axes | extracted_axes))

    # Merge biomarker roles
    existing_roles = set(result.get("biomarker_role", []))
    extracted_roles = set(extracted.get("biomarker_roles", []))
    result["biomarker_role"] = sorted(list(existing_roles | extracted_roles))

    return result


def _normalize_list(value: Any) -> list[str]:
    """Normalize a value to a list of strings."""
    if not value:
        return []
    if isinstance(value, list):
        return [str(item).strip() for item in value if item and str(item).strip()]
    return [str(value).strip()] if str(value).strip() else []


def _normalize_int(value: Any) -> int | None:
    """Normalize a value to an integer or None."""
    if value is None:
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def _empty_metadata() -> Dict[str, Any]:
    """Return empty metadata structure."""
    return {
        "diseases": [],
        "targets": [],
        "signatures": [],
        "phenotype_axes": [],
        "biomarker_roles": [],
        "sample_count": None,
        "study_design": [],
        "timepoints": [],
        "cell_lines": [],
        "tissue_types": [],
        "species": [],
        "assay_types": [],
        "therapeutic_area": [],
    }

