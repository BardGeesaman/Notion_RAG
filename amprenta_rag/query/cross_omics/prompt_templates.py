"""
Prompt templates for cross-omics LLM synthesis.

Contains the system prompt and helper functions for building user prompts
for cross-omics summary generation.
"""

from __future__ import annotations

# Maximum context chunks to send to LLM
MAX_CONTEXT_CHUNKS = 50
MAX_CHUNK_LENGTH = 2000


def get_cross_omics_system_prompt() -> str:
    """
    Get the system prompt for cross-omics summary synthesis.
    
    Returns:
        System prompt string for OpenAI API
    """
    return (
        "You are an expert in multi-omics analysis (lipidomics, metabolomics, proteomics, transcriptomics). "
        "You are given snippets of information about multi-omics evidence.\n\n"
        "Summarize the cross-omics evidence with the following structure:\n\n"
        "1. High-level context\n"
        "2. Per-omics findings:\n"
        "   - Lipidomics\n"
        "   - Metabolomics\n"
        "   - Proteomics\n"
        "   - Transcriptomics\n"
        "3. Cross-omics convergence:\n"
        "   - Features (genes/proteins/metabolites/lipids) that consistently change across multiple omics.\n"
        "4. Cross-omics divergence:\n"
        "   - Conflicting signals or modality-specific changes.\n"
        "5. Disease, model system, and matrix context.\n"
        "6. Key open questions and next experimental steps.\n\n"
        "Only use information provided in the context. Do not hallucinate external facts. "
        "Label modality-specific findings clearly. Be concise but comprehensive."
    )


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

