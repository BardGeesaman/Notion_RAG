from datetime import datetime
from typing import List, Literal
from uuid import UUID

from pydantic import BaseModel


class ChatMessage(BaseModel):
    role: Literal["user", "assistant", "system"]
    content: str
    timestamp: datetime


class ChatTurn(BaseModel):
    messages: List[ChatMessage]  # usually [user_msg, assistant_msg]


class ChatSessionState(BaseModel):
    id: UUID
    created_at: datetime
    turns: List[ChatTurn] = []
    # optional: current_focus (e.g. dataset_id, etc)
