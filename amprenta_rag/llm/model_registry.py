"""Model registry for multi-model intelligence."""
from __future__ import annotations

import os
from typing import List, Dict, Any

from amprenta_rag.clients.openai_client import get_openai_client
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

AVAILABLE_MODELS = {
    "gpt-4o": {
        "provider": "openai",
        "name": "gpt-4o",
        "description": "OpenAI GPT-4o",
    },
    "gpt-4o-mini": {
        "provider": "openai",
        "name": "gpt-4o-mini",
        "description": "OpenAI GPT-4o Mini (faster)",
    },
    "gpt-3.5-turbo": {
        "provider": "openai",
        "name": "gpt-3.5-turbo",
        "description": "OpenAI GPT-3.5",
    },
    "claude-3-5-sonnet": {
        "provider": "anthropic",
        "name": "claude-3-5-sonnet-20241022",
        "description": "Anthropic Claude 3.5 Sonnet",
    },
    "gemini-1.5-pro": {
        "provider": "google",
        "name": "gemini-1.5-pro",
        "description": "Google Gemini 1.5 Pro",
    },
}


def get_available_models() -> List[Dict[str, Any]]:
    """
    Get list of available models.
    
    Returns:
        List of model info dictionaries
    """
    return list(AVAILABLE_MODELS.values())


def get_model_client(model_name: str):
    """
    Get client for a specific model.
    
    Args:
        model_name: Name of the model (e.g., "gpt-4o")
        
    Returns:
        Model client (OpenAI, Anthropic, or Google client)
        
    Raises:
        ValueError: If model not found or API key missing
    """
    if model_name not in AVAILABLE_MODELS:
        raise ValueError(f"Model '{model_name}' not found. Available models: {list(AVAILABLE_MODELS.keys())}")
    
    model_info = AVAILABLE_MODELS[model_name]
    provider = model_info["provider"]
    
    if provider == "openai":
        return get_openai_client()
    elif provider == "anthropic":
        api_key = os.getenv("ANTHROPIC_API_KEY")
        if not api_key:
            raise ValueError("ANTHROPIC_API_KEY environment variable not set")
        try:
            import anthropic
            return anthropic.Anthropic(api_key=api_key)
        except ImportError:
            raise ValueError("anthropic package not installed. Install with: pip install anthropic")
    elif provider == "google":
        api_key = os.getenv("GOOGLE_API_KEY")
        if not api_key:
            raise ValueError("GOOGLE_API_KEY environment variable not set")
        try:
            import google.generativeai as genai
            genai.configure(api_key=api_key)
            return genai
        except ImportError:
            raise ValueError("google-generativeai package not installed. Install with: pip install google-generativeai")
    else:
        raise ValueError(f"Unsupported provider: {provider}")


def call_model(model_name: str, messages: List[Dict[str, str]], temperature: float = 0.2) -> str:
    """
    Call a model with messages and return response text.
    
    Args:
        model_name: Name of the model to use
        messages: List of message dicts with "role" and "content"
        temperature: Temperature for generation (default: 0.2)
        
    Returns:
        Response text from the model
        
    Raises:
        ValueError: If model not found or API key missing
        Exception: If API call fails
    """
    if model_name not in AVAILABLE_MODELS:
        raise ValueError(f"Model '{model_name}' not found. Available models: {list(AVAILABLE_MODELS.keys())}")
    
    model_info = AVAILABLE_MODELS[model_name]
    provider = model_info["provider"]
    actual_model_name = model_info["name"]
    
    if provider == "openai":
        client = get_openai_client()
        try:
            response = client.chat.completions.create(
                model=actual_model_name,
                messages=messages,
                temperature=temperature,
            )
            return response.choices[0].message.content.strip()  # type: ignore[union-attr]
        except Exception as e:
            logger.error("[LLM] Error calling model %s: %s", model_name, e)
            raise
    elif provider == "anthropic":
        api_key = os.getenv("ANTHROPIC_API_KEY")
        if not api_key:
            raise ValueError("ANTHROPIC_API_KEY environment variable not set. Please set it to use Claude models.")
        
        try:
            import anthropic
        except ImportError:
            raise ValueError("anthropic package not installed. Install with: pip install anthropic")
        
        client = anthropic.Anthropic(api_key=api_key)
        try:
            # Convert messages format for Anthropic (they use similar format)
            anthropic_messages = []
            system_message = None
            for msg in messages:
                role = msg.get("role", "user")
                content = msg.get("content", "")
                if role == "system":
                    system_message = content
                else:
                    # Anthropic uses "user" and "assistant" roles
                    anthropic_role = "user" if role == "user" else "assistant"
                    anthropic_messages.append({"role": anthropic_role, "content": content})
            
            response = client.messages.create(
                model=actual_model_name,
                messages=anthropic_messages,
                max_tokens=2000,
                temperature=temperature,
                system=system_message if system_message else "",
            )
            # Anthropic returns text in content[0].text
            return response.content[0].text.strip()  # type: ignore[union-attr]
        except Exception as e:
            logger.error("[LLM] Error calling Claude model %s: %s", model_name, e)
            raise
    elif provider == "google":
        api_key = os.getenv("GOOGLE_API_KEY")
        if not api_key:
            raise ValueError("GOOGLE_API_KEY environment variable not set. Please set it to use Gemini models.")
        
        try:
            import google.generativeai as genai
        except ImportError:
            raise ValueError("google-generativeai package not installed. Install with: pip install google-generativeai")
        
        genai.configure(api_key=api_key)
        model = genai.GenerativeModel(actual_model_name)
        
        try:
            # Convert messages to a single prompt for Gemini
            # Gemini uses a simpler prompt format, so we'll combine messages
            prompt_parts = []
            for msg in messages:
                role = msg.get("role", "user")
                content = msg.get("content", "")
                if role == "system":
                    prompt_parts.append(f"System: {content}")
                elif role == "user":
                    prompt_parts.append(f"User: {content}")
                elif role == "assistant":
                    prompt_parts.append(f"Assistant: {content}")
            
            prompt = "\n\n".join(prompt_parts)
            
            # Configure generation parameters
            generation_config = genai.types.GenerationConfig(
                temperature=temperature,
            )
            
            response = model.generate_content(
                prompt,
                generation_config=generation_config,
            )
            return response.text.strip()
        except Exception as e:
            logger.error("[LLM] Error calling Gemini model %s: %s", model_name, e)
            raise
    else:
        raise ValueError(f"Unsupported provider: {provider}")
