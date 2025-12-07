# Reset Context Template

This file defines how the Architect should reset or refresh working context while preserving the Agent Team Charter and all persistent memory files.

The reset process must **never delete**, overwrite, or alter:

* `agents/session-memory.md`
* `agents/*` agent instruction files
* `Agent Team Charter`
* Any project files not explicitly designated for reset

Resetting only affects **ephemeral working context**, not long-term memory.

---

## 1. Purpose of Reset

A context reset should be used when:

* The working memory has become noisy or too long
* The Architect’s understanding of the current state may be stale
* The user is switching machines, branches, or large work areas
* A major refactor or redesign requires a fresh, clean conceptual state
* The user explicitly requests a clean-slate perspective

---

## 2. Reset Command (User → Architect)

When the user wants a reset, they will say:

```
Architect:
reset context but preserve all agent definitions, the Agent Team Charter, and session-memory.md. Then rehydrate your context using the session memory and generate a fresh project summary.
```

Alternative phrasing the Architect should also accept:

```
Architect: context reset
```

or

```
Architect:
perform a clean context refresh
```

In all cases, the Architect must:

* Preserve all persistent information
* Drop only transient working memory
* Rebuild a clean internal model based on safe sources

---

## 3. Reset Procedure (Architect Responsibilities)

When receiving a reset request, the Architect must:

### Step 1 — Acknowledge

Confirm the reset is understood:

```
FROM: Architect
TO: User

Acknowledging context reset request. Preserving all persistent memory files and agent definitions.

END OF MESSAGE
FROM: Architect
TO: User
```

### Step 2 — Discard Ephemeral Context

Clear internal assumptions, cached plans, partial task sequences, stale states, or conversation-specific artifacts.

### Step 3 — Reload Persistent Memory

Re-read:

* `agents/session-memory.md`
* Any additional continuity files (e.g., `agents/history/*`)
* The Agent Team Charter (implicit)

### Step 4 — Rehydrate Working Model

Reconstruct an up-to-date internal state including:

* Current tasks
* Decisions
* Architectural state
* Open questions
* Roadmap
* Notes and continuity summary

### Step 5 — Generate a Fresh Summary for the User

Produce a clean, structured summary:

* Current known project state
* Tasks in progress
* Outstanding decisions
* Recommended next actions

This summary becomes the fresh starting point after reset.

---

## 4. Architect Reset Output Format

After a complete reset, Architect must send:

```
FROM: Architect
TO: User

Context Reset Completed

1. Summary of Preserved Information
   - Agent Team Charter
   - session-memory.md
   - Decision log
   - Active and completed tasks
   - Any relevant history files

2. Rehydrated Project State
   - Current architecture view
   - Active tasks
   - Dependencies or blockers
   - Pending decisions

3. Recommended Next Steps
   - <list>

4. Notes
   - <clarifications or options>

END OF MESSAGE
FROM: Architect
TO: User
```

---

## 5. Notes & Guarantees

* Resetting context **must not** erase progress.
* Resetting improves clarity, reduces noise, and restores a stable baseline.
* Architect must never guess about state not found in files—missing items must be surfaced to the user.
* Resetting is always *safe* and *reversible* by revisiting persistent memory.

---

## 6. Optional (but supported) Fast Reset Command

The Architect should support a minimal fast-reset command:

```
Architect:
hard reset context
```

This performs the same process but in a more compact form, suitable for rapid cleanup.

---

## 7. Future Extensions

This file may be expanded later to support:

* Partial resets (e.g., reset testing context only)
* Branch-specific resets
* Automated or scheduled resets
* Multi-project switching workflows

Architect should adapt accordingly when new templates are added.
