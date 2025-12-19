"""Hallucination detection and groundedness checking for RAG answers."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class GroundednessResult:
    """Result of groundedness check."""
    score: float  # 0.0 to 1.0, where 1.0 is fully grounded
    unsupported_claims: List[str]
    supported_claims: List[str]


def check_groundedness(answer: str, chunks: List[str]) -> GroundednessResult:
    """
    Check if an answer is grounded in the provided source chunks.

    Args:
        answer: The generated answer to check
        chunks: List of source chunks that should support the answer

    Returns:
        GroundednessResult with score and lists of supported/unsupported claims
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    context = "\n\n---\n\n".join(chunks[:8])  # Cap context

    prompt = (
        "Analyze the following answer and determine which claims are supported by the provided sources.\n\n"
        f"Sources:\n{context}\n\n"
        f"Answer to check:\n{answer}\n\n"
        "For each factual claim in the answer, determine if it is:\n"
        "1. SUPPORTED - directly stated or clearly inferable from the sources\n"
        "2. UNSUPPORTED - not found in the sources or contradicts them\n\n"
        "Respond in this exact format:\n"
        "SUPPORTED:\n- [claim 1]\n- [claim 2]\n\n"
        "UNSUPPORTED:\n- [claim 1]\n- [claim 2]\n\n"
        "If a claim is partially supported but has unsupported details, list it under UNSUPPORTED."
    )

    try:
        logger.info("[HALLUCINATION] Checking groundedness for answer (length=%d)", len(answer))
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,  # Low temperature for consistent parsing
        )
        response_text = resp.choices[0].message.content.strip()  # type: ignore[union-attr]

        # Parse response
        supported_claims: List[str] = []
        unsupported_claims: List[str] = []

        current_section = None
        for line in response_text.split("\n"):
            line = line.strip()
            if not line:
                continue
            if line.upper().startswith("SUPPORTED:"):
                current_section = "supported"
            elif line.upper().startswith("UNSUPPORTED:"):
                current_section = "unsupported"
            elif line.startswith("-"):
                claim = line[1:].strip()
                if current_section == "supported":
                    supported_claims.append(claim)
                elif current_section == "unsupported":
                    unsupported_claims.append(claim)

        # Calculate score: supported / total claims
        total_claims = len(supported_claims) + len(unsupported_claims)
        if total_claims == 0:
            score = 1.0  # No claims to check, assume grounded
        else:
            score = len(supported_claims) / total_claims

        logger.info(
            "[HALLUCINATION] Score=%.2f (%d supported, %d unsupported)",
            score,
            len(supported_claims),
            len(unsupported_claims),
        )

        return GroundednessResult(
            score=score,
            unsupported_claims=unsupported_claims,
            supported_claims=supported_claims,
        )
    except Exception as e:
        logger.error("[HALLUCINATION] Error checking groundedness: %r", e)
        # Return neutral result on error
        return GroundednessResult(
            score=0.5,
            unsupported_claims=[],
            supported_claims=[],
        )
