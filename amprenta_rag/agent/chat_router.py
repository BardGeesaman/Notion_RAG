from typing import Literal

ChatIntent = Literal[
    "freeform_rag",
    "dataset_summary",
    "program_summary",
    "signature_summary",
    "feature_summary",
    "similar_datasets",
    "help",
]


def route_intent(user_text: str) -> ChatIntent:
    text = user_text.lower()
    if "similar dataset" in text or "similar study" in text:
        return "similar_datasets"
    if "summarize dataset" in text or "tell me about dataset" in text:
        return "dataset_summary"
    if "program" in text and ("summary" in text or "overview" in text):
        return "program_summary"
    if "signature" in text:
        return "signature_summary"
    if "help" in text or "what can you do" in text:
        return "help"
    return "freeform_rag"
