from uuid import UUID
import re
from datetime import datetime
from typing import Tuple

from amprenta_rag.query.cross_omics_reasoning import (
    cross_omics_dataset_summary_postgres,
    cross_omics_feature_summary_postgres,
    cross_omics_program_summary_postgres,
    cross_omics_signature_summary_postgres,
)
from amprenta_rag.query.rag_engine import query_rag

from .chat_router import route_intent
from .chat_types import ChatMessage, ChatSessionState, ChatTurn


def run_chat_turn(session: ChatSessionState, user_text: str) -> Tuple[ChatSessionState, str]:
    intent = route_intent(user_text)
    # For v1: require explicit IDs in question
    try:
        if intent == "freeform_rag":
            rag_result = query_rag(user_query=user_text, top_k=10, use_postgres=True)
            answer = rag_result.answer or "I couldn't find enough evidence to answer that."
        elif intent == "dataset_summary":
            # Extract dataset_id explicitly (v1: expects UUID after 'dataset:')
            m = re.search(r"dataset:\s*([0-9a-f-]+)", user_text, re.I)
            if m:
                dataset_uuid = UUID(m.group(1))
                summary = cross_omics_dataset_summary_postgres(dataset_uuid)
                answer = summary
            else:
                answer = "Please specify a dataset ID as 'dataset: <UUID>'."
        elif intent == "program_summary":
            m = re.search(r"program:\s*([-\w]+)", user_text, re.I)
            if m:
                program_uuid = UUID(m.group(1))
                answer = cross_omics_program_summary_postgres(program_uuid)
            else:
                answer = "Please specify a program ID as 'program: <ID>'."
        elif intent == "signature_summary":
            m = re.search(r"signature:\s*([-\w]+)", user_text, re.I)
            if m:
                sig_uuid = UUID(m.group(1))
                answer = cross_omics_signature_summary_postgres(sig_uuid)
            else:
                answer = "Please specify a signature ID as 'signature: <ID>'."
        elif intent == "feature_summary":
            m = re.search(r"feature:\s*([\w]+)", user_text, re.I)
            if m:
                feat_uuid = UUID(m.group(1))
                answer = cross_omics_feature_summary_postgres(feat_uuid, "gene")  # v1: default to gene
            else:
                answer = "Please specify a feature as 'feature: <NAME>'."
        elif intent == "similar_datasets":
            answer = "Similarity search is coming soon."
        elif intent == "help":
            answer = (
                "I can help you with:\n"
                "- Freeform questions about your data (RAG search)\n"
                "- Summaries of programs, datasets, features, and signatures\n"
                "- Finding similar datasets\n"
                "- Generating evidence reports\n"
                "Try e.g. 'Summarize dataset <ID>' or 'Find datasets similar to <ID>'."
            )
        else:
            answer = "I'm not sure how to handle that request yet. Try rephrasing or ask for 'help'."
    except Exception as e:
        answer = f"Error: {e}"
    session.turns.append(
        ChatTurn(
            messages=[
                ChatMessage(role="user", content=user_text, timestamp=datetime.utcnow()),
                ChatMessage(role="assistant", content=answer, timestamp=datetime.utcnow()),
            ]
        )
    )
    return session, answer
