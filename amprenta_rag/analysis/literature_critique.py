"""Literature critical analysis service."""
from __future__ import annotations

import json
from typing import List, Dict, Any

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def generate_critique(text: str) -> Dict[str, Any]:
    """
    Generate a critical analysis of scientific text.

    Args:
        text: Scientific text to analyze

    Returns:
        Dictionary with strengths, weaknesses, limitations, and methodology_score
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    prompt = f"""Analyze the following scientific text with rigorous scientific evaluation.
Provide a critical analysis including strengths, weaknesses, limitations, and a methodology score (0-100).

Text to analyze:
{text}

Respond with a JSON object containing:
- "strengths": array of strings describing methodological or analytical strengths
- "weaknesses": array of strings describing methodological or analytical weaknesses
- "limitations": array of strings describing study limitations
- "methodology_score": integer from 0-100 rating the overall methodology quality

Be thorough and critical. Focus on scientific rigor, experimental design, statistical methods, and reproducibility."""

    try:
        response = client.chat.completions.create(
            model=chat_model,
            messages=[
                {"role": "system", "content": "You are a rigorous scientific reviewer. Provide detailed, critical analysis in JSON format."},
                {"role": "user", "content": prompt}
            ],
            response_format={"type": "json_object"},
            temperature=0.3,
        )

        result = json.loads(response.choices[0].message.content)

        # Ensure all required fields exist
        return {
            "strengths": result.get("strengths", []),
            "weaknesses": result.get("weaknesses", []),
            "limitations": result.get("limitations", []),
            "methodology_score": result.get("methodology_score", 50),
        }
    except Exception as e:
        logger.error("[LITERATURE] Failed to generate critique: %s", e)
        return {
            "strengths": [],
            "weaknesses": [],
            "limitations": [],
            "methodology_score": 50,
        }


def extract_unanswered_questions(text: str) -> List[str]:
    """
    Extract unanswered questions and gaps from scientific text.

    Args:
        text: Scientific text to analyze

    Returns:
        List of question strings identifying gaps and open questions
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    prompt = f"""Analyze the following scientific text and identify gaps, unanswered questions, and areas requiring further research.

Text to analyze:
{text}

Respond with a JSON object containing:
- "questions": array of strings, each representing an unanswered question or research gap

Focus on:
- Questions that the study raises but doesn't answer
- Gaps in the methodology or analysis
- Areas where further research is needed
- Unresolved contradictions or inconsistencies"""

    try:
        response = client.chat.completions.create(
            model=chat_model,
            messages=[
                {"role": "system", "content": "You are a scientific reviewer identifying research gaps. Provide questions in JSON format."},
                {"role": "user", "content": prompt}
            ],
            response_format={"type": "json_object"},
            temperature=0.3,
        )

        result = json.loads(response.choices[0].message.content)
        return result.get("questions", [])
    except Exception as e:
        logger.error("[LITERATURE] Failed to extract questions: %s", e)
        return []


def detect_contradictions(text1: str, text2: str) -> List[Dict[str, str]]:
    """
    Detect contradictions between two scientific texts.

    Args:
        text1: First scientific text
        text2: Second scientific text

    Returns:
        List of dictionaries with claim_a, claim_b, topic, and severity
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    prompt = f"""Compare the following two scientific texts and identify any contradictory claims or conflicting findings.

Text 1:
{text1}

Text 2:
{text2}

Respond with a JSON object containing:
- "contradictions": array of objects, each with:
  - "claim_a": string describing the claim from text 1
  - "claim_b": string describing the contradictory claim from text 2
  - "topic": string describing the topic/domain of the contradiction
  - "severity": string ("low", "medium", or "high") indicating how significant the contradiction is

Focus on factual contradictions, conflicting experimental results, or opposing conclusions."""

    try:
        response = client.chat.completions.create(
            model=chat_model,
            messages=[
                {"role": "system", "content": "You are a scientific reviewer comparing texts for contradictions. Provide analysis in JSON format."},
                {"role": "user", "content": prompt}
            ],
            response_format={"type": "json_object"},
            temperature=0.3,
        )

        result = json.loads(response.choices[0].message.content)
        contradictions = result.get("contradictions", [])

        # Validate structure
        validated = []
        for c in contradictions:
            if isinstance(c, dict) and all(k in c for k in ["claim_a", "claim_b", "topic", "severity"]):
                validated.append({
                    "claim_a": str(c["claim_a"]),
                    "claim_b": str(c["claim_b"]),
                    "topic": str(c["topic"]),
                    "severity": str(c["severity"]),
                })

        return validated
    except Exception as e:
        logger.error("[LITERATURE] Failed to detect contradictions: %s", e)
        return []
