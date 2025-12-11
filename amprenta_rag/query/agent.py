"""Agentic RAG for complex multi-step queries."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Dict, Any

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class AgenticResult:
    """Result from agentic RAG query."""
    answer: str
    steps: List[Dict[str, Any]] = field(default_factory=list)
    total_chunks: int = 0


def analyze_question(question: str) -> Dict[str, Any]:
    """
    Analyze a question to determine if it's complex and needs decomposition.
    
    Args:
        question: The user's question
        
    Returns:
        Dict with is_complex (bool) and sub_questions (List[str] if complex)
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()
    
    prompt = (
        "Analyze this question and determine if it needs to be broken down into sub-questions.\n\n"
        f"Question: {question}\n\n"
        "A complex question requires multiple pieces of information or multiple steps to answer.\n"
        "A simple question can be answered with a single search.\n\n"
        "Respond in JSON format:\n"
        '{"is_complex": true/false, "sub_questions": ["sub-q1", "sub-q2"] if complex, else []}'
    )
    
    try:
        logger.debug("[AGENT] Analyzing question: %s", question[:50])
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
            response_format={"type": "json_object"},
        )
        import json
        result = json.loads(resp.choices[0].message.content.strip())  # type: ignore[union-attr]
        logger.info("[AGENT] Question analysis: complex=%s, sub_questions=%d", 
                   result.get("is_complex", False), len(result.get("sub_questions", [])))
        return result
    except Exception as e:
        logger.error("[AGENT] Error analyzing question: %r", e)
        return {"is_complex": False, "sub_questions": []}


def should_continue(question: str, current_context: str) -> bool:
    """
    Determine if we have enough information to answer the question.
    
    Args:
        question: The original question
        current_context: Accumulated context from RAG searches
        
    Returns:
        True if we need more information, False if we have enough
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()
    
    prompt = (
        "Evaluate if we have enough information to answer this question.\n\n"
        f"Question: {question}\n\n"
        f"Current Context:\n{current_context[:2000]}\n\n"
        "Respond with JSON:\n"
        '{"has_enough": true/false, "reason": "brief explanation"}'
    )
    
    try:
        logger.debug("[AGENT] Checking if should continue")
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
            response_format={"type": "json_object"},
        )
        import json
        result = json.loads(resp.choices[0].message.content.strip())  # type: ignore[union-attr]
        has_enough = result.get("has_enough", True)
        logger.info("[AGENT] Should continue: %s (%s)", not has_enough, result.get("reason", ""))
        return not has_enough  # Return True if we need more (not enough)
    except Exception as e:
        logger.error("[AGENT] Error checking should_continue: %r", e)
        return False  # Default to having enough


def refine_query(original: str, context: str) -> str:
    """
    Suggest a refined query based on information gaps.
    
    Args:
        original: The original query
        context: Current context gathered
        
    Returns:
        Refined query string
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()
    
    prompt = (
        "Based on the original question and current context, suggest a refined query to fill information gaps.\n\n"
        f"Original Question: {original}\n\n"
        f"Current Context:\n{context[:1500]}\n\n"
        "What additional information do we need? Provide a refined query that would help find it.\n"
        "Respond with ONLY the refined query text, no explanation."
    )
    
    try:
        logger.debug("[AGENT] Refining query")
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3,
        )
        refined = resp.choices[0].message.content.strip()  # type: ignore[union-attr]
        logger.info("[AGENT] Refined query: %s", refined[:50])
        return refined
    except Exception as e:
        logger.error("[AGENT] Error refining query: %r", e)
        return original


def agentic_rag(question: str, max_steps: int = 3) -> AgenticResult:
    """
    Perform agentic RAG with multi-step reasoning and query refinement.
    
    Args:
        question: The user's question
        max_steps: Maximum number of refinement steps
        
    Returns:
        AgenticResult with answer, steps, and total chunks used
    """
    # Lazy import to avoid circular dependency
    from amprenta_rag.query.rag.query import query_rag
    
    logger.info("[AGENT] Starting agentic RAG for: %s", question[:50])
    
    # Step 1: Analyze question
    analysis = analyze_question(question)
    is_complex = analysis.get("is_complex", False)
    sub_questions = analysis.get("sub_questions", [])
    
    steps = []
    all_chunks = []
    total_chunks = 0
    
    if not is_complex:
        # Simple question: single RAG call
        logger.info("[AGENT] Simple question, performing single RAG call")
        result = query_rag(question, generate_answer=True)
        steps.append({
            "step": 1,
            "type": "simple_rag",
            "query": question,
            "chunks_used": len(result.context_chunks),
        })
        all_chunks.extend(result.context_chunks)
        total_chunks += len(result.context_chunks)
        
        # Final synthesis
        from amprenta_rag.query.rag.synthesis import synthesize_answer
        answer = synthesize_answer(question, all_chunks)
        
        return AgenticResult(
            answer=answer,
            steps=steps,
            total_chunks=total_chunks,
        )
    
    # Complex question: decompose and search
    logger.info("[AGENT] Complex question, decomposing into %d sub-questions", len(sub_questions))
    
    all_context = []
    step_count = 1
    
    for sub_q in sub_questions:
        logger.info("[AGENT] Processing sub-question %d/%d: %s", step_count, len(sub_questions), sub_q[:50])
        
        # Search for this sub-question
        result = query_rag(sub_q, generate_answer=False)
        chunks = result.context_chunks
        all_chunks.extend(chunks)
        total_chunks += len(chunks)
        all_context.append(f"Sub-question: {sub_q}\nContext: {' '.join(chunks[:3])}")
        
        steps.append({
            "step": step_count,
            "type": "sub_question",
            "query": sub_q,
            "chunks_used": len(chunks),
        })
        step_count += 1
        
        # Check if we need more info
        current_context_str = "\n\n".join(all_context)
        if should_continue(sub_q, current_context_str) and step_count <= max_steps:
            refined = refine_query(sub_q, current_context_str)
            logger.info("[AGENT] Refining and searching again: %s", refined[:50])
            
            # Search with refined query
            refined_result = query_rag(refined, generate_answer=False)
            refined_chunks = refined_result.context_chunks
            all_chunks.extend(refined_chunks)
            total_chunks += len(refined_chunks)
            all_context.append(f"Refined query: {refined}\nContext: {' '.join(refined_chunks[:3])}")
            
            steps.append({
                "step": step_count,
                "type": "refined_search",
                "query": refined,
                "chunks_used": len(refined_chunks),
            })
            step_count += 1
    
    # Final synthesis across all gathered context
    logger.info("[AGENT] Synthesizing final answer from %d total chunks", total_chunks)
    from amprenta_rag.query.rag.synthesis import synthesize_answer
    final_context = "\n\n".join(all_context)
    answer = synthesize_answer(question, all_chunks)
    
    steps.append({
        "step": step_count,
        "type": "final_synthesis",
        "query": question,
        "chunks_used": total_chunks,
    })
    
    logger.info("[AGENT] Agentic RAG complete: %d steps, %d chunks", len(steps), total_chunks)
    
    return AgenticResult(
        answer=answer,
        steps=steps,
        total_chunks=total_chunks,
    )
