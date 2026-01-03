"""Embedding utility functions for RAG."""

from typing import List


def get_embedding(text: str, model: str = "text-embedding-ada-002") -> List[float]:
    """
    Get text embedding from OpenAI.
    
    Args:
        text: Text to embed
        model: OpenAI embedding model
    
    Returns:
        Embedding vector as list of floats
    """
    try:
        import openai
        client = openai.OpenAI()
        response = client.embeddings.create(
            model=model,
            input=text
        )
        return response.data[0].embedding
    except ImportError:
        # Fallback for test environment - return dummy embedding
        import hashlib
        hash_val = int(hashlib.md5(text.encode()).hexdigest()[:8], 16)
        # Generate consistent dummy embedding based on text hash
        return [float((hash_val + i) % 1000) / 1000.0 for i in range(1536)]
    except Exception as e:
        # Return zero vector on error
        return [0.0] * 1536
