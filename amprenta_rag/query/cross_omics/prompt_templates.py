"""
Prompt templates for cross-omics LLM synthesis.

Contains the system prompt and helper functions for building user prompts
for cross-omics summary generation with enhanced context awareness.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

# Maximum context chunks to send to LLM
MAX_CONTEXT_CHUNKS = 50
MAX_CHUNK_LENGTH = 2000


def get_cross_omics_system_prompt(
    include_comparative: bool = False,
) -> str:
    """
    Get the system prompt for cross-omics summary synthesis.

    Args:
        include_comparative: Whether to include comparative analysis instructions

    Returns:
        System prompt string for OpenAI API
    """
    base_prompt = (
        "You are an expert in multi-omics analysis (lipidomics, metabolomics, proteomics, transcriptomics). "
        "You are given snippets of information about multi-omics evidence.\n\n"
        "Summarize the cross-omics evidence with the following structure:\n\n"
        "1. High-level context\n"
        "   - Disease context (if available): What diseases or conditions are being studied?\n"
        "   - Model system context: Are these in vitro, in vivo, or patient-derived samples?\n"
        "   - Matrix context: What sample types (CSF, plasma, serum, tissue, etc.)?\n"
        "2. Per-omics findings:\n"
        "   - Lipidomics: Key lipid species changes, classes affected\n"
        "   - Metabolomics: Key metabolite changes, pathways affected\n"
        "   - Proteomics: Key protein changes, functional categories\n"
        "   - Transcriptomics: Key gene expression changes, pathways\n"
        "3. Cross-omics convergence:\n"
        "   - Features (genes/proteins/metabolites/lipids) that consistently change across multiple omics.\n"
        "   - Pathways or biological processes showing convergent signals.\n"
        "4. Cross-omics divergence:\n"
        "   - Conflicting signals or modality-specific changes.\n"
        "   - Omics-specific patterns that don't align with other modalities.\n"
    )

    if include_comparative:
        base_prompt += (
            "5. Comparative analysis:\n"
            "   - Disease vs control comparisons (if applicable).\n"
            "   - Direction of changes (increased/decreased) across omics.\n"
            "   - Magnitude and consistency of changes.\n"
        )
    else:
        base_prompt += "5. Disease, model system, and matrix context (detailed).\n"

    base_prompt += (
        "6. Key open questions and next experimental steps.\n\n"
        "Only use information provided in the context. Do not hallucinate external facts. "
        "Label modality-specific findings clearly. Be concise but comprehensive. "
        "When disease, model system, or matrix context is provided, integrate it throughout the summary."
    )

    return base_prompt


def build_enhanced_prompt(
    entity_name: str,
    entity_type: str,
    context_info: Optional[Dict[str, Any]] = None,
    omics_counts: Optional[Dict[str, int]] = None,
    additional_info: Optional[str] = None,
) -> str:
    """
    Build an enhanced user prompt with context information.

    Args:
        entity_name: Name of the entity being summarized
        entity_type: Type of entity ("program", "dataset", "signature", "feature")
        context_info: Dictionary with disease, matrix, model_systems context
        omics_counts: Dictionary with counts per omics type
        additional_info: Additional information to include

    Returns:
        Formatted user prompt string
    """
    prompt_parts: List[str] = [f"Generate a cross-omics summary for the {entity_type}: {entity_name}"]

    # Add context information
    if context_info:
        context_str = format_context_for_prompt(context_info)
        if context_str:
            prompt_parts.append(f"\nContext:\n{context_str}")

    # Add omics counts
    if omics_counts:
        omics_lines = [
            f"- {omics}: {count} chunks"
            for omics, count in omics_counts.items()
            if count > 0
        ]
        if omics_lines:
            prompt_parts.append(f"\nChunks retrieved:\n" + "\n".join(omics_lines))

    # Add additional information
    if additional_info:
        prompt_parts.append(f"\n{additional_info}")

    return "\n".join(prompt_parts)


def format_context_for_prompt(context: Dict[str, Any]) -> str:
    """
    Format context dictionary for inclusion in prompts.

    Args:
        context: Dictionary with diseases, matrix, model_systems

    Returns:
        Formatted string
    """
    parts: List[str] = []

    if context.get("diseases"):
        diseases_str = ", ".join(context["diseases"])
        parts.append(f"Disease context: {diseases_str}")

    if context.get("matrix"):
        matrix_str = ", ".join(context["matrix"])
        parts.append(f"Matrix types: {matrix_str}")

    if context.get("model_systems"):
        model_str = ", ".join(context["model_systems"])
        parts.append(f"Model systems: {model_str}")

    return "\n".join(parts) if parts else "No context information available."


def prepare_context_chunks(
    context_chunks: list[str],
    max_chunks: int = MAX_CONTEXT_CHUNKS,
) -> str:
    """
    Prepare context chunks for LLM synthesis.

    Args:
        context_chunks: List of chunk text strings
        max_chunks: Maximum number of chunks to include

    Returns:
        Formatted context string
    """
    # Truncate chunks to fit token budget
    truncated_chunks = context_chunks[:max_chunks]

    # Limit each chunk length
    limited_chunks = [
        chunk[:MAX_CHUNK_LENGTH] + ("..." if len(chunk) > MAX_CHUNK_LENGTH else "")
        for chunk in truncated_chunks
    ]

    return "\n\n---\n\n".join(limited_chunks)

