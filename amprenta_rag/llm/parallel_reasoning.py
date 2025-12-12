"""Parallel reasoning across multiple LLM models."""
from __future__ import annotations

from typing import List, Dict, Any, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

from amprenta_rag.llm.model_registry import call_model, get_available_models, AVAILABLE_MODELS
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def run_parallel_models(prompt: str, models: List[str]) -> List[Dict[str, Any]]:
    """
    Run multiple models in parallel with the same prompt.
    
    Args:
        prompt: The prompt to send to all models
        models: List of model names to use
        
    Returns:
        List of dicts with keys: "model", "response", "success", "error"
    """
    results: List[Dict[str, Any]] = []
    
    def call_single_model(model_name: str) -> Dict[str, Any]:
        """Call a single model and return result dict."""
        try:
            messages = [{"role": "user", "content": prompt}]
            response = call_model(model_name, messages, temperature=0.2)
            return {
                "model": model_name,
                "response": response,
                "success": True,
                "error": None,
            }
        except Exception as e:
            logger.error("[PARALLEL] Error calling model %s: %s", model_name, e)
            return {
                "model": model_name,
                "response": None,
                "success": False,
                "error": str(e),
            }
    
    # Run models in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=len(models)) as executor:
        future_to_model = {executor.submit(call_single_model, model): model for model in models}
        
        for future in as_completed(future_to_model):
            result = future.result()
            results.append(result)
    
    # Sort results to match input order
    model_order = {model: idx for idx, model in enumerate(models)}
    results.sort(key=lambda x: model_order.get(x["model"], 999))
    
    return results


def synthesize_responses(question: str, responses: List[Dict[str, Any]]) -> str:
    """
    Synthesize multiple model responses into a single answer.
    
    Uses GPT-4o to analyze agreements/disagreements and create a combined answer.
    
    Args:
        question: The original question
        responses: List of response dicts from run_parallel_models
        
    Returns:
        Synthesized answer text
    """
    # Filter to successful responses only
    successful_responses = [r for r in responses if r.get("success") and r.get("response")]
    
    if not successful_responses:
        return "No successful model responses to synthesize."
    
    if len(successful_responses) == 1:
        return successful_responses[0]["response"]
    
    # Build synthesis prompt
    responses_text = "\n\n".join(
        f"**Model: {r['model']}**\n{r['response']}"
        for r in successful_responses
    )
    
    synthesis_prompt = (
        f"Question: {question}\n\n"
        f"Multiple AI models have provided answers to this question:\n\n"
        f"{responses_text}\n\n"
        f"Please synthesize these responses into a single, comprehensive answer.\n"
        f"- Note areas where models agree\n"
        f"- Note areas where models disagree or provide different perspectives\n"
        f"- Combine the best insights from each response\n"
        f"- If there are contradictions, acknowledge them and explain the different viewpoints\n"
        f"- Provide a clear, well-structured final answer\n\n"
        f"Synthesized Answer:"
    )
    
    try:
        messages = [{"role": "user", "content": synthesis_prompt}]
        synthesized = call_model("gpt-4o", messages, temperature=0.2)
        return synthesized
    except Exception as e:
        logger.error("[PARALLEL] Error synthesizing responses: %s", e)
        # Fallback: return the first successful response
        return successful_responses[0]["response"]


def parallel_query(question: str, models: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Run a query across multiple models in parallel and synthesize the results.
    
    Args:
        question: The question to ask all models
        models: Optional list of model names. If None, uses all available models.
        
    Returns:
        Dict with keys:
            - "individual_responses": List of response dicts from each model
            - "synthesis": Synthesized combined answer
    """
    if models is None:
        # Use all available models
        available = get_available_models()
        models = [m["name"] for m in available]
    
    if not models:
        return {
            "individual_responses": [],
            "synthesis": "No models available.",
        }
    
    logger.info("[PARALLEL] Running parallel query across %d models: %s", len(models), models)
    
    # Run models in parallel
    individual_responses = run_parallel_models(question, models)
    
    # Synthesize responses
    synthesis = synthesize_responses(question, individual_responses)
    
    return {
        "individual_responses": individual_responses,
        "synthesis": synthesis,
    }

