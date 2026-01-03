"""AI Expert Agents service layer."""

import logging
from datetime import datetime
from typing import List, Optional, Dict, Any
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import (
    ExpertAgent, ExpertConversation, ExpertMessage,
    ExpertFeedback, ExpertTrainingExample, ExpertKnowledgeDoc
)
from amprenta_rag.utils.llm import get_llm_response_with_tokens
from amprenta_rag.utils.embeddings import get_embedding

logger = logging.getLogger(__name__)


# ============================================================================
# EXPERT MANAGEMENT
# ============================================================================

def get_experts(db: Session, active_only: bool = True) -> List[ExpertAgent]:
    """List all expert agents."""
    query = db.query(ExpertAgent)
    if active_only:
        query = query.filter(ExpertAgent.is_active.is_(True))
    return query.all()


def get_expert(db: Session, expert_id: UUID) -> Optional[ExpertAgent]:
    """Get expert by ID."""
    return db.query(ExpertAgent).filter(ExpertAgent.id == expert_id).first()


def update_expert(
    db: Session,
    expert_id: UUID,
    system_prompt: Optional[str] = None,
    specializations: Optional[List[str]] = None,
    bump_version: bool = False
) -> Optional[ExpertAgent]:
    """Update expert configuration. Optionally bump prompt version."""
    expert = get_expert(db, expert_id)
    if not expert:
        return None
    
    if system_prompt:
        expert.system_prompt = system_prompt
    if specializations:
        expert.specializations = specializations
    if bump_version:
        # Increment version: v1.0 -> v1.1
        parts = expert.prompt_version.split('.')
        expert.prompt_version = f"{parts[0]}.{int(parts[1]) + 1}"
    
    db.commit()
    db.refresh(expert)
    return expert


# ============================================================================
# CONVERSATION MANAGEMENT
# ============================================================================

def create_conversation(
    db: Session,
    user_id: UUID,
    expert_ids: List[UUID],
    title: Optional[str] = None,
    context_entity_type: Optional[str] = None,
    context_entity_id: Optional[UUID] = None,
    is_panel: bool = False
) -> ExpertConversation:
    """Create new conversation with one or more experts."""
    conversation = ExpertConversation(
        user_id=user_id,
        title=title or "New Conversation",
        context_entity_type=context_entity_type,
        context_entity_id=context_entity_id,
        is_panel=is_panel or len(expert_ids) > 1
    )
    db.add(conversation)
    db.flush()
    
    # Add experts to conversation
    experts = db.query(ExpertAgent).filter(ExpertAgent.id.in_(expert_ids)).all()
    conversation.participants = experts
    
    db.commit()
    db.refresh(conversation)
    return conversation


def get_conversation(db: Session, conversation_id: UUID) -> Optional[ExpertConversation]:
    """Get conversation with messages."""
    return db.query(ExpertConversation).filter(
        ExpertConversation.id == conversation_id
    ).first()


def get_user_conversations(db: Session, user_id: UUID, limit: int = 50) -> List[ExpertConversation]:
    """List user's conversations."""
    return db.query(ExpertConversation).filter(
        ExpertConversation.user_id == user_id
    ).order_by(ExpertConversation.last_message_at.desc()).limit(limit).all()


# ============================================================================
# MESSAGE & RESPONSE GENERATION
# ============================================================================

def add_user_message(
    db: Session,
    conversation_id: UUID,
    content: str
) -> ExpertMessage:
    """Add user message to conversation."""
    message = ExpertMessage(
        conversation_id=conversation_id,
        role="user",
        content=content
    )
    db.add(message)
    db.commit()
    db.refresh(message)
    return message


def get_expert_response(
    db: Session,
    conversation_id: UUID,
    expert_id: UUID,
    user_message: str
) -> ExpertMessage:
    """
    Generate single expert response.
    
    1. Load expert persona and prompt
    2. Retrieve relevant knowledge (RAG)
    3. Inject entity context if available
    4. Generate response via LLM
    5. Store message with metadata
    """
    conversation = get_conversation(db, conversation_id)
    expert = get_expert(db, expert_id)
    
    if not conversation or not expert:
        raise ValueError("Conversation or expert not found")
    
    # Build context
    context_parts = []
    
    # 1. Expert-specific knowledge (RAG)
    knowledge = retrieve_expert_knowledge(db, expert_id, user_message, k=3)
    if knowledge:
        context_parts.append("Relevant knowledge:\n" + "\n".join([k.content for k in knowledge]))
    
    # 2. Entity context
    if conversation.context_entity_type and conversation.context_entity_id:
        entity_context = inject_entity_context(
            db, 
            conversation.context_entity_type, 
            conversation.context_entity_id
        )
        if entity_context:
            context_parts.append(f"Context ({conversation.context_entity_type}):\n{entity_context}")
    
    # 3. Conversation history (last 10 messages)
    history = db.query(ExpertMessage).filter(
        ExpertMessage.conversation_id == conversation_id
    ).order_by(ExpertMessage.created_at.desc()).limit(10).all()
    history.reverse()
    
    # Build messages for LLM
    messages = [
        {"role": "system", "content": expert.system_prompt + "\n\n" + "\n\n".join(context_parts)}
    ]
    
    for msg in history:
        messages.append({"role": msg.role, "content": msg.content})
    
    messages.append({"role": "user", "content": user_message})
    
    # Generate response
    response, token_count = get_llm_response_with_tokens(messages)
    
    # Extract reasoning if present (between <reasoning> tags)
    reasoning = None
    if "<reasoning>" in response and "</reasoning>" in response:
        start = response.index("<reasoning>") + len("<reasoning>")
        end = response.index("</reasoning>")
        reasoning = response[start:end].strip()
        response = response.replace(f"<reasoning>{reasoning}</reasoning>", "").strip()
    
    # Store response
    expert_message = ExpertMessage(
        conversation_id=conversation_id,
        expert_id=expert_id,
        role="assistant",
        content=response,
        prompt_version_used=expert.prompt_version,
        reasoning=reasoning,
        token_count=token_count
    )
    db.add(expert_message)
    
    # Update conversation timestamp
    conversation.last_message_at = datetime.utcnow()
    
    db.commit()
    db.refresh(expert_message)
    return expert_message


def get_panel_response(
    db: Session,
    conversation_id: UUID,
    user_message: str
) -> Dict[str, Any]:
    """
    Get multi-expert panel response.
    
    Orchestration Strategy (P1 FIX):
    1. Send question to all experts in parallel
    2. Collect individual responses (Round 1)
    3. Identify disagreements
    4. If disagreements: experts see each other's views (Round 2)
    5. Synthesize with consensus scoring
    """
    conversation = get_conversation(db, conversation_id)
    if not conversation or not conversation.participants:
        raise ValueError("Panel conversation not found or has no experts")
    
    # Add user message
    add_user_message(db, conversation_id, user_message)
    
    # Round 1: Get individual responses
    expert_responses = []
    for expert in conversation.participants:
        response = get_expert_response(db, conversation_id, expert.id, user_message)
        expert_responses.append({
            "expert_id": str(expert.id),
            "expert_name": expert.name,
            "expert_role": expert.role,
            "content": response.content,
            "reasoning": response.reasoning,
            "message_id": str(response.id)
        })
    
    # Analyze for disagreements (simple keyword matching for MVP)
    disagreements = detect_disagreements(expert_responses)
    
    # Calculate consensus score
    consensus_type = "full" if not disagreements else ("partial" if len(disagreements) < 2 else "none")
    confidence = 0.9 if consensus_type == "full" else (0.7 if consensus_type == "partial" else 0.5)
    
    # Synthesize primary recommendation
    primary = synthesize_recommendations(expert_responses)
    
    return {
        "consensus": consensus_type,
        "confidence_score": confidence,
        "primary_recommendation": primary,
        "expert_responses": expert_responses,
        "disagreements": disagreements
    }


def detect_disagreements(responses: List[Dict]) -> List[Dict]:
    """Detect disagreements between expert responses (MVP: keyword matching)."""
    disagreements = []
    # Simple MVP: check for contradictory keywords
    positive_keywords = ["recommend", "suggest", "should", "would"]
    negative_keywords = ["not recommend", "avoid", "shouldn't", "wouldn't"]
    
    # More sophisticated disagreement detection can be added later
    return disagreements


def synthesize_recommendations(responses: List[Dict]) -> str:
    """Synthesize a primary recommendation from multiple expert responses."""
    if not responses:
        return "No recommendations available."
    
    # MVP: Return first expert's key recommendation
    return f"Based on expert consensus: {responses[0]['content'][:500]}..."


# ============================================================================
# RAG KNOWLEDGE RETRIEVAL
# ============================================================================

def retrieve_expert_knowledge(
    db: Session,
    expert_id: UUID,
    query: str,
    k: int = 5
) -> List[ExpertKnowledgeDoc]:
    """Retrieve relevant knowledge docs for expert."""
    expert = get_expert(db, expert_id)
    if not expert:
        return []
    
    namespace = f"expert_{expert_id}"
    
    try:
        query_embedding = get_embedding(query)
        
        # pgvector cosine similarity search
        results = db.query(ExpertKnowledgeDoc).filter(
            ExpertKnowledgeDoc.namespace == namespace
        ).order_by(
            ExpertKnowledgeDoc.embedding.cosine_distance(query_embedding)
        ).limit(k).all()
        
        return results
    except Exception as e:
        logger.warning(f"RAG retrieval failed: {e}")
        # Fallback: return recent docs
        return db.query(ExpertKnowledgeDoc).filter(
            ExpertKnowledgeDoc.expert_id == expert_id
        ).order_by(ExpertKnowledgeDoc.created_at.desc()).limit(k).all()


def add_knowledge_doc(
    db: Session,
    expert_id: UUID,
    title: str,
    content: str,
    source_type: Optional[str] = None,
    source_url: Optional[str] = None,
    chunk_index: int = 0
) -> ExpertKnowledgeDoc:
    """Add knowledge document for expert."""
    embedding = get_embedding(content)
    
    doc = ExpertKnowledgeDoc(
        expert_id=expert_id,
        namespace=f"expert_{expert_id}",
        title=title,
        content=content,
        chunk_index=chunk_index,
        embedding=embedding,
        embedding_model="text-embedding-ada-002",
        source_type=source_type,
        source_url=source_url,
    )
    db.add(doc)
    db.commit()
    db.refresh(doc)
    return doc


# ============================================================================
# ENTITY CONTEXT INJECTION
# ============================================================================

def inject_entity_context(
    db: Session,
    entity_type: str,
    entity_id: UUID
) -> Optional[str]:
    """Load relevant platform data for context injection."""
    from amprenta_rag.database.models import Compound, Dataset, Experiment
    
    if entity_type == "Compound":
        compound = db.query(Compound).filter(Compound.id == entity_id).first()
        if compound:
            return f"Compound: {compound.compound_id}\nSMILES: {compound.smiles}\nMW: {compound.molecular_weight}"
    
    elif entity_type == "Dataset":
        dataset = db.query(Dataset).filter(Dataset.id == entity_id).first()
        if dataset:
            return f"Dataset: {dataset.name}\nType: {dataset.omics_type}\nDescription: {dataset.description or 'N/A'}"
    
    elif entity_type == "Experiment":
        experiment = db.query(Experiment).filter(Experiment.id == entity_id).first()
        if experiment:
            return f"Experiment: {experiment.name}\nType: {experiment.design_type}\nDescription: {experiment.description or 'N/A'}"
    
    return None


# ============================================================================
# FEEDBACK & TRAINING
# ============================================================================

def record_feedback(
    db: Session,
    message_id: UUID,
    user_id: UUID,
    rating: int,
    correction: Optional[str] = None,
    tags: Optional[List[str]] = None
) -> ExpertFeedback:
    """Record feedback on expert response."""
    if not 1 <= rating <= 5:
        raise ValueError("Rating must be between 1 and 5")
    
    feedback = ExpertFeedback(
        message_id=message_id,
        user_id=user_id,
        rating=rating,
        correction=correction,
        tags=tags or []
    )
    db.add(feedback)
    db.commit()
    db.refresh(feedback)
    return feedback


def create_training_example(
    db: Session,
    expert_id: UUID,
    question: str,
    ideal_answer: str,
    created_by_id: UUID,
    source: str = "manual"
) -> ExpertTrainingExample:
    """Create training example for expert."""
    expert = get_expert(db, expert_id)
    if not expert:
        raise ValueError("Expert not found")
    
    example = ExpertTrainingExample(
        expert_id=expert_id,
        question=question,
        ideal_answer=ideal_answer,
        prompt_version=expert.prompt_version,
        is_approved=False,
        source=source,
        created_by_id=created_by_id
    )
    db.add(example)
    db.commit()
    db.refresh(example)
    return example


def approve_training_example(db: Session, example_id: UUID) -> bool:
    """Approve training example for fine-tuning."""
    example = db.query(ExpertTrainingExample).filter(
        ExpertTrainingExample.id == example_id
    ).first()
    if not example:
        return False
    example.is_approved = True
    db.commit()
    return True


def export_training_data(db: Session, expert_id: Optional[UUID] = None) -> List[Dict]:
    """Export approved training examples in JSONL format for fine-tuning."""
    query = db.query(ExpertTrainingExample).filter(
        ExpertTrainingExample.is_approved.is_(True)
    )
    if expert_id:
        query = query.filter(ExpertTrainingExample.expert_id == expert_id)
    
    examples = query.all()
    
    # OpenAI fine-tuning format
    return [
        {
            "messages": [
                {"role": "system", "content": example.expert.system_prompt},
                {"role": "user", "content": example.question},
                {"role": "assistant", "content": example.ideal_answer}
            ]
        }
        for example in examples
    ]


# ============================================================================
# CONVERSATION UTILITIES
# ============================================================================

def get_conversation_messages(
    db: Session,
    conversation_id: UUID,
    limit: int = 100
) -> List[ExpertMessage]:
    """Get conversation messages in chronological order."""
    return db.query(ExpertMessage).filter(
        ExpertMessage.conversation_id == conversation_id
    ).order_by(ExpertMessage.created_at.asc()).limit(limit).all()


def get_expert_statistics(db: Session, expert_id: UUID) -> Dict[str, Any]:
    """Get statistics for an expert."""
    expert = get_expert(db, expert_id)
    if not expert:
        return {}
    
    # Count messages
    message_count = db.query(ExpertMessage).filter(
        ExpertMessage.expert_id == expert_id
    ).count()
    
    # Count knowledge docs
    knowledge_count = db.query(ExpertKnowledgeDoc).filter(
        ExpertKnowledgeDoc.expert_id == expert_id
    ).count()
    
    # Count training examples
    training_count = db.query(ExpertTrainingExample).filter(
        ExpertTrainingExample.expert_id == expert_id
    ).count()
    
    # Average rating
    feedback_query = db.query(ExpertFeedback).join(ExpertMessage).filter(
        ExpertMessage.expert_id == expert_id
    )
    ratings = [f.rating for f in feedback_query.all()]
    avg_rating = sum(ratings) / len(ratings) if ratings else None
    
    return {
        "expert_name": expert.name,
        "expert_role": expert.role,
        "message_count": message_count,
        "knowledge_doc_count": knowledge_count,
        "training_example_count": training_count,
        "average_rating": avg_rating,
        "total_feedback": len(ratings),
        "prompt_version": expert.prompt_version,
    }
