# Quick Start Guide for New Chat Sessions

**Purpose**: Use this when starting a new Cursor chat to quickly restore context.

---

## 🚀 Fast Setup (2 Steps)

### Step 1: Upload Context File

In your new Cursor chat, upload this file:
```
context/MASTER_CONTEXT_FOR_NEW_CHAT.md
```

### Step 2: Paste This Message

```
I'm continuing work on the Amprenta Multi-Omics Platform. I've uploaded 
the master context document (MASTER_CONTEXT_FOR_NEW_CHAT.md).

Please:
1. Read and confirm you understand the system architecture
2. Confirm you understand the current implementation status  
3. Tell me what you see as the next logical steps based on the roadmap

The system uses Notion as canonical source of truth, all ingestion is 
idempotent, and errors should be non-blocking. All work should align 
with HIGH_LEVEL_STRATEGY_CONTEXT.md principles.
```

---

## 📋 Alternative: Quick Context Paste

If you prefer not to upload a file, paste this in the new chat:

```
CONTEXT: I'm working on the Amprenta Multi-Omics Platform. Here's the 
essential context:

KEY PRINCIPLES:
- Notion is canonical source of truth
- All operations are idempotent
- Error handling is non-blocking
- Consistent logging prefixes

CURRENT STATUS:
- Multi-omics ingestion pipelines: ✅ Complete (all 4 types)
- Feature linking: ✅ Complete (all 4 types)
- Multi-omics signatures: ✅ Complete
- Signature scoring: ✅ Complete
- Cross-omics RAG reasoning: ✅ Complete

KEY FILES:
- HIGH_LEVEL_STRATEGY_CONTEXT.md - Strategic principles
- UNIFIED_STRATEGIC_ROADMAP.md - Complete roadmap (5 tiers)
- amprenta_rag/config.py - Configuration

NEXT STEPS: See UNIFIED_STRATEGIC_ROADMAP.md Tier 1 items:
1. Feature Caching (performance)
2. Batch Ingestion (operations)
3. Enhanced Cross-Omics Reasoning

Please acknowledge understanding and confirm what you see as next steps.
```

---

## ✅ What You Should See

The AI should:
- ✅ Acknowledge understanding of the system
- ✅ Reference the strategic context
- ✅ Understand current implementation state
- ✅ Propose next steps aligned with roadmap
- ✅ Maintain architectural patterns

---

## 📄 Full Context Document

For comprehensive context, always reference:
- `context/MASTER_CONTEXT_FOR_NEW_CHAT.md` - Complete master context
- `context/HIGH_LEVEL_STRATEGY_CONTEXT.md` - Strategic principles
- `context/UNIFIED_STRATEGIC_ROADMAP.md` - Complete roadmap

---

## 🔄 When Chat Crashes or Needs Restart

1. Open new chat
2. Upload `MASTER_CONTEXT_FOR_NEW_CHAT.md`
3. Use the message template above
4. Continue seamlessly!

---

**That's it! You're ready to continue work in any new chat session.**

