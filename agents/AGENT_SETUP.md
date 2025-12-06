# Amprenta – Cursor Agent Setup (Dec 2025)

This document captures the **intended configuration** for all Cursor agents in the Amprenta RAG workspace.

The goal is **predictable, single-model behavior per agent** (except Architect), and a clear separation of responsibilities.

---

## Summary Table

| Agent      | Mode   | Model                          | Max Mode | Instructions file          |
|----------- |--------|--------------------------------|----------|----------------------------|
| Architect  | Plan   | Composer 1                     | ON       | `agents/architect.md`      |
| Implementor| Agent  | gpt-5.1-codex-max-high-fast    | OFF      | `agents/implementor.md`    |
| Reviewer   | Agent  | Sonnet 4.5 (normal)            | OFF      | `agents/reviewer.md`       |
| Testing    | Agent  | Gemini 3 Pro                   | OFF      | `agents/testing.md`        |
| Automator  | Agent  | Sonnet 4.5 (normal)            | OFF      | `agents/automator.md`      |

---

## Architect

**Role:** System architect for multi-omics + HTS + RAG. Designs schemas, flows, and plans. Does *not* edit code.

- Mode: **Plan**
- Model: **Composer 1 only**
- Max Mode: **ON**
- Other models: **OFF**
- Instruction source: `agents/architect.md`

### Behavior expectations

- Produces phased plans, schemas, DB designs, RAG flows.
- Never applies code edits or runs tools.
- Hands plans to Implementor.

---

## Implementor

**Role:** Implements Architect’s plans in code (Python, SQLAlchemy, Alembic, etc.).

- Mode: **Agent**
- Model: **`gpt-5.1-codex-max-high-fast`**
- Max Mode: **OFF**
- Other models: **OFF**
- Instruction source: `agents/implementor.md`

### Behavior expectations

- Writes small, safe diffs.
- Follows Architect’s plan exactly.
- Asks clarifying questions before large refactors.
- Does not redesign architecture.

---

## Reviewer

**Role:** Reviews diffs and plans for correctness, safety, style.

- Mode: **Agent**
- Model: **Sonnet 4.5 (normal, *not* Thinking)**
- Max Mode: **OFF**
- Other models: **OFF**
- Instruction source: `agents/reviewer.md`

### Behavior expectations

- Summarizes changes.
- Flags logical errors, missing tests, risky migrations.
- Suggests improvements without large rewrites.

---

## Testing (QA)

**Role:** Acts as a QA engineer for the RAG system + Streamlit UI.

- Mode: **Agent**
- Model: **Gemini 3 Pro**
- Max Mode: **OFF**
- Other models: **OFF**
- Instruction source: `agents/testing.md`

### Behavior expectations

- Simulates user flows and reports UX / behavior issues.
- Proposes test cases and fixtures.
- Does not directly modify code.

---

## Automator

**Role:** Mechanical repo housekeeping – move/rename files, clean docs, generate small scripts.

- Mode: **Agent**
- Model: **Sonnet 4.5 (normal)**
- Max Mode: **OFF**
- Other models: **OFF**
- Instruction source: `agents/automator.md`

### Behavior expectations

- Performs repetitive, pattern-based changes.
- Never deletes or refactors core logic without explicit confirmation.
- Asks before destructive operations.

---

## How to re-instantiate agents on a new machine

1. Clone the repo and create your Python environment.
2. Open `RAG` in Cursor.
3. For each agent:
   - Open its `agents/*.md` file.
   - Copy the entire file.
   - Paste as the first message in that agent’s chat.
   - Configure the mode + single model according to the table above.

This file is the “source of truth” for agent behavior.