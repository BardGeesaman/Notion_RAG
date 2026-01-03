"""AI Expert Agents service layer with CrewAI integration."""

import logging
from datetime import datetime
from typing import List, Optional, Dict, Any
from uuid import UUID

from crewai import Agent, Task, Crew
from crewai.tools import BaseTool
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from amprenta_rag.database.models import (
    ExpertAgent, ExpertConversation, ExpertMessage,
    ExpertFeedback, ExpertTrainingExample, ExpertKnowledgeDoc
)

logger = logging.getLogger(__name__)


# ============================================================================
# CREWAI TOOLS
# ============================================================================

class CompoundLookupTool(BaseTool):
    """Tool to look up compound information."""
    name: str = "compound_lookup"
    description: str = "Look up compound by ID or SMILES to get properties, structure, and activity data"
    
    def _run(self, compound_id: str = None, smiles: str = None) -> str:
        """Execute compound lookup."""
        from amprenta_rag.database.session import db_session
        from amprenta_rag.models.chemistry import Compound
        
        with db_session() as db:
            if compound_id:
                compound = db.query(Compound).filter(Compound.compound_id == compound_id).first()
            elif smiles:
                compound = db.query(Compound).filter(Compound.smiles == smiles).first()
            else:
                return "Error: Must provide compound_id or smiles"
            
            if not compound:
                return "Compound not found"
            
            return f"""Compound: {compound.compound_id}
SMILES: {compound.smiles}
Molecular Weight: {compound.molecular_weight}
LogP: {compound.logp}
HBD Count: {compound.hbd_count}
HBA Count: {compound.hba_count}"""


class RAGSearchTool(BaseTool):
    """Tool to search expert knowledge base."""
    name: str = "rag_search"
    description: str = "Search expert-specific knowledge documents for relevant information"
    
    def __init__(self, expert_id: UUID):
        super().__init__()
        self.expert_id = expert_id
    
    def _run(self, query: str) -> str:
        """Search expert knowledge base."""
        from amprenta_rag.database.session import db_session
        from amprenta_rag.utils.embeddings import get_embedding
        
        with db_session() as db:
            try:
                query_embedding = get_embedding(query)
                namespace = f"expert_{self.expert_id}"
                
                docs = db.query(ExpertKnowledgeDoc).filter(
                    ExpertKnowledgeDoc.namespace == namespace
                ).order_by(
                    ExpertKnowledgeDoc.embedding.cosine_distance(query_embedding)
                ).limit(3).all()
                
                if not docs:
                    return "No relevant knowledge found"
                
                return "\n\n".join([f"Source: {doc.title}\n{doc.content}" for doc in docs])
            except Exception as e:
                return f"Knowledge search failed: {e}"


class DatasetInfoTool(BaseTool):
    """Tool to get dataset information."""
    name: str = "dataset_info"
    description: str = "Get detailed information about a dataset including samples, features, and metadata"
    
    def _run(self, dataset_id: str) -> str:
        """Get dataset information."""
        from amprenta_rag.database.session import db_session
        from amprenta_rag.database.models import Dataset, Feature
        
        with db_session() as db:
            try:
                dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                if not dataset:
                    return f"Dataset {dataset_id} not found"
                
                feature_count = db.query(Feature).filter(Feature.dataset_id == dataset_id).count()
                
                return f"""Dataset: {dataset.name}
Type: {dataset.omics_type}
Description: {dataset.description or 'N/A'}
Sample Count: {dataset.sample_count}
Feature Count: {feature_count}
Source: {dataset.source or 'N/A'}"""
            except Exception as e:
                return f"Dataset lookup failed: {e}"


class ExperimentInfoTool(BaseTool):
    """Tool to get experiment information."""
    name: str = "experiment_info"
    description: str = "Get detailed information about an experiment including design, datasets, and results"
    
    def _run(self, experiment_id: str) -> str:
        """Get experiment information."""
        from amprenta_rag.database.session import db_session
        from amprenta_rag.database.models import Experiment, Dataset
        
        with db_session() as db:
            try:
                experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()
                if not experiment:
                    return f"Experiment {experiment_id} not found"
                
                dataset_count = db.query(Dataset).filter(Dataset.experiment_id == experiment_id).count()
                
                return f"""Experiment: {experiment.name}
Design Type: {experiment.design_type}
Description: {experiment.description or 'N/A'}
Status: {experiment.status}
Dataset Count: {dataset_count}
Created: {experiment.created_at.strftime('%Y-%m-%d') if experiment.created_at else 'N/A'}"""
            except Exception as e:
                return f"Experiment lookup failed: {e}"


# ============================================================================
# CREWAI AGENT CREATION
# ============================================================================

def create_expert_agent(expert: ExpertAgent) -> Agent:
    """Build CrewAI agent from database ExpertAgent model."""
    tools = [
        CompoundLookupTool(),
        RAGSearchTool(expert.id),
        DatasetInfoTool(),
        ExperimentInfoTool(),
    ]
    
    agent = Agent(
        role=expert.role,
        goal=f"Provide expert {expert.role.lower()} advice and analysis",
        backstory=expert.system_prompt,
        tools=tools,
        verbose=True,
        memory=True,  # Enable CrewAI memory for session context
    )
    
    return agent


def create_expert_crew(experts: List[ExpertAgent], process: str = "sequential") -> Crew:
    """Build CrewAI crew for panel discussions."""
    agents = [create_expert_agent(expert) for expert in experts]
    
    crew = Crew(
        agents=agents,
        process=process,
        verbose=True,
        memory=True,
    )
    
    return crew


# ============================================================================
# RESPONSE GENERATION
# ============================================================================

def get_expert_response(
    db: Session,
    conversation_id: UUID,
    expert_id: UUID,
    user_message: str
) -> ExpertMessage:
    """
    Generate single expert response using CrewAI.
    
    Memory Strategy:
    - CrewAI memory for in-session context
    - Load recent DB messages as task context
    - Persist response to ExpertMessage table
    """
    # Get expert and conversation
    expert = db.query(ExpertAgent).filter(ExpertAgent.id == expert_id).first()
    conversation = db.query(ExpertConversation).filter(ExpertConversation.id == conversation_id).first()
    
    if not expert or not conversation:
        raise ValueError("Expert or conversation not found")
    
    # Load recent conversation history
    recent_messages = db.query(ExpertMessage).filter(
        ExpertMessage.conversation_id == conversation_id
    ).order_by(ExpertMessage.created_at.desc()).limit(10).all()
    recent_messages.reverse()
    
    # Build context
    context_parts = []
    
    # Add entity context if available
    if conversation.context_entity_type and conversation.context_entity_id:
        entity_context = inject_entity_context(
            db, conversation.context_entity_type, conversation.context_entity_id
        )
        if entity_context:
            context_parts.append(f"Context ({conversation.context_entity_type}):\n{entity_context}")
    
    # Add recent conversation history
    if recent_messages:
        history = "\n".join([f"{msg.role}: {msg.content}" for msg in recent_messages])
        context_parts.append(f"Recent conversation:\n{history}")
    
    # Create CrewAI agent and task
    agent = create_expert_agent(expert)
    
    task_description = user_message
    if context_parts:
        task_description = "\n\n".join(context_parts) + f"\n\nUser Question: {user_message}"
    
    task = Task(
        description=task_description,
        agent=agent,
        expected_output="A detailed expert response addressing the user's question"
    )
    
    # Execute task
    crew = Crew(agents=[agent], tasks=[task], verbose=False)
    result = crew.kickoff()
    
    # Extract response content
    response_content = str(result)
    
    # Estimate token count (rough approximation)
    token_count = len(response_content) // 4
    
    # Store response in database
    expert_message = ExpertMessage(
        conversation_id=conversation_id,
        expert_id=expert_id,
        role="assistant",
        content=response_content,
        prompt_version_used=expert.prompt_version,
        token_count=token_count
    )
    
    db.add(expert_message)
    
    # Update conversation timestamp
    conversation.last_message_at = datetime.utcnow()
    
    db.commit()
    db.refresh(expert_message)
    
    logger.info(f"Generated response from {expert.name} for conversation {conversation_id}")
    return expert_message


def get_panel_response(
    db: Session,
    conversation_id: UUID,
    expert_ids: List[UUID],
    user_message: str
) -> Dict[str, Any]:
    """
    Generate multi-expert panel response using CrewAI.
    
    CrewAI handles the orchestration between experts.
    """
    # Get experts and conversation
    experts = db.query(ExpertAgent).filter(ExpertAgent.id.in_(expert_ids)).all()
    conversation = db.query(ExpertConversation).filter(ExpertConversation.id == conversation_id).first()
    
    if not experts or not conversation:
        raise ValueError("Experts or conversation not found")
    
    # Add user message to conversation
    user_msg = ExpertMessage(
        conversation_id=conversation_id,
        role="user",
        content=user_message
    )
    db.add(user_msg)
    db.flush()
    
    # Load context
    context_parts = []
    if conversation.context_entity_type and conversation.context_entity_id:
        entity_context = inject_entity_context(
            db, conversation.context_entity_type, conversation.context_entity_id
        )
        if entity_context:
            context_parts.append(f"Context: {entity_context}")
    
    # Create CrewAI crew
    crew = create_expert_crew(experts, process="hierarchical")
    
    # Create collaborative task
    task_description = user_message
    if context_parts:
        task_description = "\n\n".join(context_parts) + f"\n\nPanel Question: {user_message}"
    
    task = Task(
        description=task_description,
        expected_output="A comprehensive panel response with expert perspectives and consensus recommendations"
    )
    
    # Execute crew
    crew.tasks = [task]
    result = crew.kickoff()
    
    # Store individual expert responses (CrewAI provides aggregated result)
    expert_responses = []
    for expert in experts:
        # Store individual response (simplified for MVP)
        expert_message = ExpertMessage(
            conversation_id=conversation_id,
            expert_id=expert.id,
            role="assistant",
            content=f"[Panel Response] {str(result)}",
            prompt_version_used=expert.prompt_version,
            token_count=len(str(result)) // 4
        )
        db.add(expert_message)
        
        expert_responses.append({
            "expert_id": str(expert.id),
            "expert_name": expert.name,
            "expert_role": expert.role,
            "content": str(result),
            "message_id": str(expert_message.id)
        })
    
    # Update conversation
    conversation.last_message_at = datetime.utcnow()
    
    db.commit()
    
    return {
        "consensus": "collaborative",
        "confidence_score": 0.8,
        "primary_recommendation": str(result),
        "expert_responses": expert_responses,
        "disagreements": []
    }


# ============================================================================
# CONTEXT INJECTION
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
    feedback_type: str,
    feedback_value: Any,
    comments: Optional[str] = None
) -> ExpertFeedback:
    """Store feedback on expert response."""
    # Convert feedback to 1-5 rating
    if feedback_type == "rating":
        rating = int(feedback_value)
        if not 1 <= rating <= 5:
            raise ValueError("Rating must be between 1 and 5")
    elif feedback_type == "thumbs":
        rating = 5 if feedback_value == "up" else 1
    else:
        rating = 3  # Default neutral
    
    feedback = ExpertFeedback(
        message_id=message_id,
        rating=rating,
        correction=comments,
        tags=[feedback_type] if feedback_type else []
    )
    
    db.add(feedback)
    db.commit()
    db.refresh(feedback)
    return feedback


def retrieve_expert_knowledge(
    db: Session,
    expert_id: UUID,
    query: str,
    k: int = 5
) -> List[ExpertKnowledgeDoc]:
    """Retrieve relevant knowledge docs for expert."""
    try:
        from amprenta_rag.utils.embeddings import get_embedding
        query_embedding = get_embedding(query)
        namespace = f"expert_{expert_id}"
        
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


def export_training_data(
    db: Session,
    expert_id: Optional[UUID] = None,
    format: str = 'jsonl'
) -> List[Dict]:
    """Export approved training examples for fine-tuning."""
    query = db.query(ExpertTrainingExample).filter(
        ExpertTrainingExample.is_approved.is_(True)
    )
    if expert_id:
        query = query.filter(ExpertTrainingExample.expert_id == expert_id)
    
    examples = query.all()
    
    if format == 'jsonl':
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
    else:
        # Raw format
        return [
            {
                "expert_name": example.expert.name,
                "question": example.question,
                "answer": example.ideal_answer,
                "prompt_version": example.prompt_version
            }
            for example in examples
        ]


# ============================================================================
# CONVERSATION MANAGEMENT
# ============================================================================

def create_conversation(
    db: Session,
    user_id: UUID,
    expert_ids: List[UUID],
    title: Optional[str] = None,
    context_entity_type: Optional[str] = None,
    context_entity_id: Optional[UUID] = None
) -> ExpertConversation:
    """Create new conversation with experts."""
    conversation = ExpertConversation(
        user_id=user_id,
        title=title or "Expert Consultation",
        context_entity_type=context_entity_type,
        context_entity_id=context_entity_id,
        is_panel=len(expert_ids) > 1,
        max_messages=50
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


# ============================================================================
# EXPERT UTILITIES
# ============================================================================

def get_experts(db: Session, active_only: bool = True) -> List[ExpertAgent]:
    """List all expert agents."""
    query = db.query(ExpertAgent)
    if active_only:
        query = query.filter(ExpertAgent.is_active.is_(True))
    return query.order_by(ExpertAgent.name).all()


def get_expert(db: Session, expert_id: UUID) -> Optional[ExpertAgent]:
    """Get expert by ID."""
    return db.query(ExpertAgent).filter(ExpertAgent.id == expert_id).first()


def get_expert_statistics(db: Session, expert_id: UUID) -> Dict[str, Any]:
    """Get performance statistics for an expert."""
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
    
    # Average rating from feedback
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
        "average_rating": avg_rating,
        "total_feedback": len(ratings),
        "prompt_version": expert.prompt_version,
        "specializations": expert.specializations,
    }
