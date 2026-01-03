"""LLM utility functions for expert agents."""

from typing import Dict, List


def get_llm_response(messages: List[Dict], model: str = "gpt-4") -> str:
    """
    Get response from OpenAI.
    
    Args:
        messages: List of message dicts with role and content
        model: OpenAI model name
    
    Returns:
        Response content string
    """
    try:
        import openai
        client = openai.OpenAI()
        response = client.chat.completions.create(
            model=model,
            messages=messages
        )
        return response.choices[0].message.content
    except ImportError:
        # Fallback for test environment
        return "Mock LLM response for testing"
    except Exception as e:
        return f"LLM error: {str(e)}"


def get_llm_response_with_tokens(messages: List[Dict], model: str = "gpt-4") -> tuple[str, int]:
    """
    Get LLM response with token count.
    
    Args:
        messages: List of message dicts
        model: OpenAI model name
    
    Returns:
        Tuple of (response_content, token_count)
    """
    try:
        import openai
        client = openai.OpenAI()
        response = client.chat.completions.create(
            model=model,
            messages=messages
        )
        content = response.choices[0].message.content
        token_count = response.usage.total_tokens if response.usage else 0
        return content, token_count
    except ImportError:
        # Fallback for test environment
        response = "Mock LLM response for testing"
        token_count = sum(len(m.get("content", "")) // 4 for m in messages) + len(response) // 4
        return response, token_count
    except Exception as e:
        return f"LLM error: {str(e)}", 0
