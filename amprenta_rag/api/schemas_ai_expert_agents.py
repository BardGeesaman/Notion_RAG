"""Schemas for AI Expert Agents API."""

from __future__ import annotations

from datetime import datetime
from typing import List, Optional, Dict, Any
from uuid import UUID

from pydantic import BaseModel, Field


class ExpertAgentResponse(BaseModel):
    """Expert agent response schema."""
    id: UUID
    name: str
    role: str
    description: Optional[str] = None
    avatar_emoji: str = "👨‍🔬"
    specializations: List[str] = []
    is_active: bool
    prompt_version: str
    
    model_config = {"from_attributes": True}


class ConversationCreate(BaseModel):
    """Create conversation request."""
    expert_ids: List[UUID] = Field(..., min_length=1, description="List of expert IDs")
    title: Optional[str] = Field(None, max_length=500)
    entity_context: Optional[Dict[str, Any]] = Field(None, description="Optional entity context")


class ConversationResponse(BaseModel):
    """Conversation response schema."""
    id: UUID
    title: str
    expert_ids: List[UUID]
    expert_names: List[str]
    is_panel: bool
    created_at: datetime
    last_message_at: Optional[datetime]
    message_count: int
    
    model_config = {"from_attributes": True}


class MessageCreate(BaseModel):
    """Create message request."""
    content: str = Field(..., min_length=1, max_length=10000)
    entity_context: Optional[Dict[str, Any]] = None


class MessageResponse(BaseModel):
    """Message response schema."""
    id: UUID
    role: str  # user, assistant
    content: str
    expert_id: Optional[UUID]
    expert_name: Optional[str]
    reasoning: Optional[str]
    token_count: Optional[int]
    created_at: datetime
    feedback_summary: Optional[Dict[str, Any]] = None
    
    model_config = {"from_attributes": True}


class PanelRequest(BaseModel):
    """Panel discussion request."""
    expert_ids: List[UUID] = Field(..., min_length=2, description="Multiple expert IDs for panel")
    question: str = Field(..., min_length=1, max_length=10000)
    entity_context: Optional[Dict[str, Any]] = None


class PanelResponse(BaseModel):
    """Panel discussion response."""
    responses: List[MessageResponse]
    consensus: Optional[str]
    confidence_score: float = Field(..., ge=0.0, le=1.0)
    disagreements: List[Dict[str, Any]] = []
    primary_recommendation: str


class FeedbackCreate(BaseModel):
    """Create feedback request."""
    rating: int = Field(..., ge=1, le=5, description="Rating from 1-5")
    correction: Optional[str] = Field(None, max_length=2000)
    tags: Optional[List[str]] = Field(None, description="Feedback tags")


class FeedbackResponse(BaseModel):
    """Feedback response schema."""
    id: UUID
    message_id: UUID
    rating: int
    correction: Optional[str]
    tags: Optional[List[str]]
    created_at: datetime
    
    model_config = {"from_attributes": True}


class ExpertStatistics(BaseModel):
    """Expert performance statistics."""
    expert_name: str
    expert_role: str
    message_count: int
    knowledge_doc_count: int
    training_example_count: int
    average_rating: Optional[float]
    total_feedback: int
    prompt_version: str
    specializations: List[str]


class ConversationDetail(BaseModel):
    """Detailed conversation with messages."""
    id: UUID
    title: str
    experts: List[ExpertAgentResponse]
    messages: List[MessageResponse]
    is_panel: bool
    context_entity_type: Optional[str]
    context_entity_id: Optional[UUID]
    created_at: datetime
    last_message_at: Optional[datetime]
