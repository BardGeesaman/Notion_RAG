"""RAG evaluation metrics for answer quality assessment."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class EvalResult:
    """RAG evaluation results."""
    faithfulness: float  # 0.0 to 1.0
    relevance: float  # 0.0 to 1.0
    context_precision: float  # 0.0 to 1.0
    overall: float  # Average of the three metrics


def evaluate_faithfulness(answer: str, chunks: List[str]) -> float:
    """
    Evaluate if the answer's claims are supported by the context chunks.

    Args:
        answer: The generated answer
        chunks: List of context chunks used to generate the answer

    Returns:
        Score from 0.0 to 1.0 (1.0 = all claims supported)
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    context = "\n\n---\n\n".join(chunks[:8])

    prompt = (
        "Evaluate if the following answer's claims are supported by the provided context.\n\n"
        f"Context:\n{context}\n\n"
        f"Answer:\n{answer}\n\n"
        "Rate the faithfulness on a scale of 0.0 to 1.0:\n"
        "- 1.0: All claims are directly supported by the context\n"
        "- 0.5: Some claims are supported, some are not\n"
        "- 0.0: No claims are supported or answer contradicts context\n\n"
        "Respond with ONLY a number between 0.0 and 1.0 (e.g., 0.85):"
    )

    try:
        logger.debug("[EVAL] Evaluating faithfulness")
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
        )
        score_text = resp.choices[0].message.content.strip()  # type: ignore[union-attr]
        score = float(score_text)
        score = max(0.0, min(1.0, score))  # Clamp to [0, 1]
        logger.debug("[EVAL] Faithfulness score: %.2f", score)
        return score
    except Exception as e:
        logger.error("[EVAL] Error evaluating faithfulness: %r", e)
        return 0.5  # Default neutral score on error


def evaluate_relevance(question: str, answer: str) -> float:
    """
    Evaluate if the answer addresses the question.

    Args:
        question: The original question
        answer: The generated answer

    Returns:
        Score from 0.0 to 1.0 (1.0 = fully relevant)
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    prompt = (
        "Evaluate if the following answer addresses the question.\n\n"
        f"Question:\n{question}\n\n"
        f"Answer:\n{answer}\n\n"
        "Rate the relevance on a scale of 0.0 to 1.0:\n"
        "- 1.0: Answer fully addresses the question\n"
        "- 0.5: Answer partially addresses the question\n"
        "- 0.0: Answer does not address the question\n\n"
        "Respond with ONLY a number between 0.0 and 1.0 (e.g., 0.90):"
    )

    try:
        logger.debug("[EVAL] Evaluating relevance")
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
        )
        score_text = resp.choices[0].message.content.strip()  # type: ignore[union-attr]
        score = float(score_text)
        score = max(0.0, min(1.0, score))  # Clamp to [0, 1]
        logger.debug("[EVAL] Relevance score: %.2f", score)
        return score
    except Exception as e:
        logger.error("[EVAL] Error evaluating relevance: %r", e)
        return 0.5  # Default neutral score on error


def evaluate_context_precision(question: str, chunks: List[str]) -> float:
    """
    Evaluate if the context chunks are relevant to the question.

    Args:
        question: The original question
        chunks: List of context chunks

    Returns:
        Score from 0.0 to 1.0 (1.0 = all chunks relevant)
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    context = "\n\n---\n\n".join(chunks[:8])

    prompt = (
        "Evaluate if the following context chunks are relevant to the question.\n\n"
        f"Question:\n{question}\n\n"
        f"Context Chunks:\n{context}\n\n"
        "Rate the context precision on a scale of 0.0 to 1.0:\n"
        "- 1.0: All chunks are highly relevant\n"
        "- 0.5: Some chunks are relevant, some are not\n"
        "- 0.0: No chunks are relevant\n\n"
        "Respond with ONLY a number between 0.0 and 1.0 (e.g., 0.75):"
    )

    try:
        logger.debug("[EVAL] Evaluating context precision")
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
        )
        score_text = resp.choices[0].message.content.strip()  # type: ignore[union-attr]
        score = float(score_text)
        score = max(0.0, min(1.0, score))  # Clamp to [0, 1]
        logger.debug("[EVAL] Context precision score: %.2f", score)
        return score
    except Exception as e:
        logger.error("[EVAL] Error evaluating context precision: %r", e)
        return 0.5  # Default neutral score on error


def evaluate_rag_response(question: str, answer: str, chunks: List[str]) -> EvalResult:
    """
    Evaluate a RAG response using multiple metrics.

    Args:
        question: The original question
        answer: The generated answer
        chunks: List of context chunks used

    Returns:
        EvalResult with all metrics and overall score
    """
    logger.info("[EVAL] Evaluating RAG response for question: %s", question[:50])

    faithfulness = evaluate_faithfulness(answer, chunks)
    relevance = evaluate_relevance(question, answer)
    context_precision = evaluate_context_precision(question, chunks)

    overall = (faithfulness + relevance + context_precision) / 3.0

    logger.info(
        "[EVAL] Scores - Faithfulness: %.2f, Relevance: %.2f, Context Precision: %.2f, Overall: %.2f",
        faithfulness,
        relevance,
        context_precision,
        overall,
    )

    return EvalResult(
        faithfulness=faithfulness,
        relevance=relevance,
        context_precision=context_precision,
        overall=overall,
    )
