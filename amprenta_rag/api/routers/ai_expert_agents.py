"""AI Expert Agents API endpoints."""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from amprenta_rag.api import schemas_ai_expert_agents as schemas
from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.database.models import User
from amprenta_rag.services import expert_agents as service

router = APIRouter()


# ============================================================================
# EXPERT BROWSING
# ============================================================================

@router.get("/", response_model=List[schemas.ExpertAgentResponse])
def list_experts(
    active_only: bool = True,
    db: Session = Depends(get_db),
):
    """List all available expert agents."""
    experts = service.get_experts(db, active_only=active_only)
    
    return [
        schemas.ExpertAgentResponse(
            id=expert.id,
            name=expert.name,
            role=expert.role,
            description=f"Specializes in: {', '.join(expert.specializations[:3])}",
            avatar_emoji="👨‍🔬" if "Dr." in expert.name else "🤖",
            specializations=expert.specializations,
            is_active=expert.is_active,
            prompt_version=expert.prompt_version,
        )
        for expert in experts
    ]


@router.get("/{expert_id}", response_model=schemas.ExpertAgentResponse)
def get_expert_details(
    expert_id: UUID,
    db: Session = Depends(get_db),
):
    """Get detailed information about a specific expert."""
    expert = service.get_expert(db, expert_id)
    if not expert:
        raise HTTPException(status_code=404, detail="Expert not found")
    
    return schemas.ExpertAgentResponse(
        id=expert.id,
        name=expert.name,
        role=expert.role,
        description=expert.system_prompt[:200] + "..." if len(expert.system_prompt) > 200 else expert.system_prompt,
        specializations=expert.specializations,
        is_active=expert.is_active,
        prompt_version=expert.prompt_version,
    )


@router.get("/{expert_id}/stats", response_model=schemas.ExpertStatistics)
def get_expert_statistics(
    expert_id: UUID,
    db: Session = Depends(get_db),
):
    """Get performance statistics for an expert."""
    stats = service.get_expert_statistics(db, expert_id)
    if not stats:
        raise HTTPException(status_code=404, detail="Expert not found")
    
    return schemas.ExpertStatistics(**stats)


# ============================================================================
# CONVERSATION MANAGEMENT
# ============================================================================

@router.post("/conversations", response_model=schemas.ConversationResponse, status_code=201)
def create_conversation(
    data: schemas.ConversationCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Start new conversation with expert(s)."""
    # Extract entity context if provided
    context_entity_type = None
    context_entity_id = None
    if data.entity_context:
        context_entity_type = data.entity_context.get("entity_type")
        context_entity_id = data.entity_context.get("entity_id")
        if context_entity_id:
            context_entity_id = UUID(context_entity_id)
    
    conversation = service.create_conversation(
        db=db,
        user_id=current_user.id,
        expert_ids=data.expert_ids,
        title=data.title,
        context_entity_type=context_entity_type,
        context_entity_id=context_entity_id,
    )
    
    # Get expert names
    expert_names = [expert.name for expert in conversation.participants]
    
    return schemas.ConversationResponse(
        id=conversation.id,
        title=conversation.title,
        expert_ids=[expert.id for expert in conversation.participants],
        expert_names=expert_names,
        is_panel=conversation.is_panel,
        created_at=conversation.created_at,
        last_message_at=conversation.last_message_at,
        message_count=len(conversation.messages),
    )


@router.get("/conversations", response_model=List[schemas.ConversationResponse])
def list_user_conversations(
    limit: int = 50,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """List user's conversations."""
    conversations = service.get_user_conversations(db, current_user.id, limit=limit)
    
    return [
        schemas.ConversationResponse(
            id=conv.id,
            title=conv.title,
            expert_ids=[expert.id for expert in conv.participants],
            expert_names=[expert.name for expert in conv.participants],
            is_panel=conv.is_panel,
            created_at=conv.created_at,
            last_message_at=conv.last_message_at,
            message_count=len(conv.messages),
        )
        for conv in conversations
    ]


@router.get("/conversations/{conversation_id}", response_model=schemas.ConversationDetail)
def get_conversation_detail(
    conversation_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get conversation with full message history."""
    conversation = service.get_conversation(db, conversation_id)
    if not conversation:
        raise HTTPException(status_code=404, detail="Conversation not found")
    
    # Check ownership
    if conversation.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    # Get messages
    messages = service.get_conversation_messages(db, conversation_id)
    
    message_responses = []
    for msg in messages:
        expert_name = None
        if msg.expert_id:
            expert = service.get_expert(db, msg.expert_id)
            expert_name = expert.name if expert else "Unknown Expert"
        
        message_responses.append(schemas.MessageResponse(
            id=msg.id,
            role=msg.role,
            content=msg.content,
            expert_id=msg.expert_id,
            expert_name=expert_name,
            reasoning=msg.reasoning,
            token_count=msg.token_count,
            created_at=msg.created_at,
        ))
    
    # Convert experts
    expert_responses = [
        schemas.ExpertAgentResponse(
            id=expert.id,
            name=expert.name,
            role=expert.role,
            specializations=expert.specializations,
            is_active=expert.is_active,
            prompt_version=expert.prompt_version,
        )
        for expert in conversation.participants
    ]
    
    return schemas.ConversationDetail(
        id=conversation.id,
        title=conversation.title,
        experts=expert_responses,
        messages=message_responses,
        is_panel=conversation.is_panel,
        context_entity_type=conversation.context_entity_type,
        context_entity_id=conversation.context_entity_id,
        created_at=conversation.created_at,
        last_message_at=conversation.last_message_at,
    )


@router.post("/conversations/{conversation_id}/messages", response_model=schemas.MessageResponse)
def send_message(
    conversation_id: UUID,
    data: schemas.MessageCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Send message to conversation and get expert response."""
    conversation = service.get_conversation(db, conversation_id)
    if not conversation:
        raise HTTPException(status_code=404, detail="Conversation not found")
    
    # Check ownership
    if conversation.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    # Add user message
    user_message = service.add_user_message(db, conversation_id, data.content)
    
    # Get response from first expert (single expert mode)
    if conversation.participants:
        expert = conversation.participants[0]
        try:
            expert_response = service.get_expert_response(
                db, conversation_id, expert.id, data.content
            )
            
            return schemas.MessageResponse(
                id=expert_response.id,
                role=expert_response.role,
                content=expert_response.content,
                expert_id=expert_response.expert_id,
                expert_name=expert.name,
                reasoning=expert_response.reasoning,
                token_count=expert_response.token_count,
                created_at=expert_response.created_at,
            )
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Expert response failed: {str(e)}")
    else:
        raise HTTPException(status_code=400, detail="No experts in conversation")


@router.delete("/conversations/{conversation_id}", status_code=204)
def archive_conversation(
    conversation_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Archive/delete conversation."""
    conversation = service.get_conversation(db, conversation_id)
    if not conversation:
        raise HTTPException(status_code=404, detail="Conversation not found")
    
    # Check ownership
    if conversation.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    # Soft delete
    conversation.is_active = False
    db.commit()


# ============================================================================
# PANEL DISCUSSIONS
# ============================================================================

@router.post("/panel", response_model=schemas.PanelResponse)
def panel_discussion(
    data: schemas.PanelRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """One-shot panel discussion without persistent conversation."""
    if len(data.expert_ids) < 2:
        raise HTTPException(status_code=400, detail="Panel requires at least 2 experts")
    
    # Create temporary conversation
    context_entity_type = None
    context_entity_id = None
    if data.entity_context:
        context_entity_type = data.entity_context.get("entity_type")
        context_entity_id = data.entity_context.get("entity_id")
        if context_entity_id:
            context_entity_id = UUID(context_entity_id)
    
    conversation = service.create_conversation(
        db=db,
        user_id=current_user.id,
        expert_ids=data.expert_ids,
        title="Panel Discussion",
        context_entity_type=context_entity_type,
        context_entity_id=context_entity_id,
        is_panel=True,
    )
    
    try:
        # Get panel response
        panel_result = service.get_panel_response(db, conversation.id, data.question)
        
        # Convert to response format
        message_responses = []
        for resp in panel_result["expert_responses"]:
            message_responses.append(schemas.MessageResponse(
                id=UUID(resp["message_id"]),
                role="assistant",
                content=resp["content"],
                expert_id=UUID(resp["expert_id"]),
                expert_name=resp["expert_name"],
                reasoning=None,
                token_count=None,
                created_at=datetime.utcnow(),
            ))
        
        return schemas.PanelResponse(
            responses=message_responses,
            consensus=panel_result["consensus"],
            confidence_score=panel_result["confidence_score"],
            disagreements=panel_result["disagreements"],
            primary_recommendation=panel_result["primary_recommendation"],
        )
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Panel discussion failed: {str(e)}")
    
    finally:
        # Clean up temporary conversation
        conversation.is_active = False
        db.commit()


# ============================================================================
# FEEDBACK
# ============================================================================

@router.post("/messages/{message_id}/feedback", response_model=schemas.FeedbackResponse, status_code=201)
def submit_feedback(
    message_id: UUID,
    data: schemas.FeedbackCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Submit feedback on expert response."""
    # Verify message exists
    from amprenta_rag.database.models import ExpertMessage
    message = db.query(ExpertMessage).filter(ExpertMessage.id == message_id).first()
    if not message:
        raise HTTPException(status_code=404, detail="Message not found")
    
    # Check if user owns the conversation
    if message.conversation.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    try:
        feedback = service.record_feedback(
            db=db,
            message_id=message_id,
            user_id=current_user.id,
            rating=data.rating,
            correction=data.correction,
            tags=data.tags,
        )
        
        return schemas.FeedbackResponse(
            id=feedback.id,
            message_id=feedback.message_id,
            rating=feedback.rating,
            correction=feedback.correction,
            tags=feedback.tags,
            created_at=feedback.created_at,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


# ============================================================================
# TRAINING & KNOWLEDGE MANAGEMENT
# ============================================================================

@router.post("/experts/{expert_id}/training-examples", response_model=schemas.TrainingExampleResponse, status_code=201)
def create_training_example(
    expert_id: UUID,
    data: schemas.TrainingExampleCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Add training example for expert (admin only)."""
    try:
        example = service.create_training_example(
            db=db,
            expert_id=expert_id,
            question=data.input_text,
            ideal_answer=data.ideal_output,
            created_by_id=current_user.id,
            source="manual"
        )
        
        return schemas.TrainingExampleResponse(
            id=example.id,
            expert_id=example.expert_id,
            question=example.question,
            ideal_answer=example.ideal_answer,
            is_approved=example.is_approved,
            prompt_version=example.prompt_version,
            source=example.source,
            created_at=example.created_at,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/experts/{expert_id}/training-examples", response_model=List[schemas.TrainingExampleResponse])
def list_training_examples(
    expert_id: UUID,
    approved_only: bool = False,
    limit: int = 100,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """List expert's training examples."""
    from amprenta_rag.database.models import ExpertTrainingExample
    
    query = db.query(ExpertTrainingExample).filter(
        ExpertTrainingExample.expert_id == expert_id
    )
    
    if approved_only:
        query = query.filter(ExpertTrainingExample.is_approved.is_(True))
    
    examples = query.order_by(ExpertTrainingExample.created_at.desc()).limit(limit).all()
    
    return [
        schemas.TrainingExampleResponse(
            id=example.id,
            expert_id=example.expert_id,
            question=example.question,
            ideal_answer=example.ideal_answer,
            is_approved=example.is_approved,
            prompt_version=example.prompt_version,
            source=example.source,
            created_at=example.created_at,
        )
        for example in examples
    ]


@router.post("/experts/{expert_id}/knowledge", response_model=schemas.KnowledgeDocResponse, status_code=201)
def upload_knowledge_document(
    expert_id: UUID,
    data: schemas.KnowledgeDocCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Upload knowledge document for expert (chunks and embeds automatically)."""
    try:
        doc = service.add_knowledge_doc(
            db=db,
            expert_id=expert_id,
            title=data.title,
            content=data.content,
            source_type=data.source_type,
            source_url=data.source_url,
        )
        
        return schemas.KnowledgeDocResponse(
            id=doc.id,
            expert_id=doc.expert_id,
            namespace=doc.namespace,
            title=doc.title,
            chunk_index=doc.chunk_index,
            source_type=doc.source_type,
            source_url=doc.source_url,
            created_at=doc.created_at,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Knowledge upload failed: {str(e)}")


@router.get("/experts/{expert_id}/feedback-summary", response_model=schemas.FeedbackSummary)
def get_feedback_summary(
    expert_id: UUID,
    db: Session = Depends(get_db),
):
    """Get aggregated feedback metrics for expert."""
    from amprenta_rag.database.models import ExpertFeedback, ExpertMessage
    
    expert = service.get_expert(db, expert_id)
    if not expert:
        raise HTTPException(status_code=404, detail="Expert not found")
    
    # Get all feedback for this expert's messages
    feedback_query = db.query(ExpertFeedback).join(ExpertMessage).filter(
        ExpertMessage.expert_id == expert_id
    )
    
    feedback_list = feedback_query.all()
    
    if not feedback_list:
        return schemas.FeedbackSummary(
            expert_id=expert_id,
            expert_name=expert.name,
            avg_rating=None,
            feedback_count=0,
            recent_corrections=[],
            rating_distribution={},
        )
    
    # Calculate metrics
    ratings = [f.rating for f in feedback_list]
    avg_rating = sum(ratings) / len(ratings)
    
    # Rating distribution
    rating_dist = {}
    for rating in ratings:
        rating_dist[rating] = rating_dist.get(rating, 0) + 1
    
    # Recent corrections
    recent_corrections = [
        f.correction for f in feedback_list[-5:] 
        if f.correction
    ]
    
    return schemas.FeedbackSummary(
        expert_id=expert_id,
        expert_name=expert.name,
        avg_rating=avg_rating,
        feedback_count=len(feedback_list),
        recent_corrections=recent_corrections,
        rating_distribution=rating_dist,
    )


@router.patch("/experts/{expert_id}", response_model=schemas.ExpertAgentResponse)
def update_expert(
    expert_id: UUID,
    data: schemas.ExpertUpdateRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Update expert configuration (admin only)."""
    try:
        expert = service.update_expert(
            db=db,
            expert_id=expert_id,
            system_prompt=data.system_prompt,
            specializations=data.specializations,
            bump_version=data.bump_version,
        )
        
        if not expert:
            raise HTTPException(status_code=404, detail="Expert not found")
        
        if data.is_active is not None:
            expert.is_active = data.is_active
            db.commit()
            db.refresh(expert)
        
        return schemas.ExpertAgentResponse(
            id=expert.id,
            name=expert.name,
            role=expert.role,
            specializations=expert.specializations,
            is_active=expert.is_active,
            prompt_version=expert.prompt_version,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Expert update failed: {str(e)}")


@router.post("/experts/{expert_id}/export-training", response_model=schemas.TrainingDataExport)
def export_training_data(
    expert_id: UUID,
    format: str = "jsonl",
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Export approved training examples for fine-tuning (admin only)."""
    try:
        training_data = service.export_training_data(db, expert_id=expert_id)
        
        return schemas.TrainingDataExport(
            format=format,
            expert_id=expert_id,
            example_count=len(training_data),
            data=training_data,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Training export failed: {str(e)}")
