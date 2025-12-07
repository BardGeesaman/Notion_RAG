# Continuity Summary Template

*Standard end-of-session summary format for the Architect*

This document defines the exact structure and content that the Architect must generate **automatically** at the end of every work session.
The summary ensures the user can resume work on any machine with full context.

Architect must produce a continuity summary whenever:

* the user ends a session
* the session naturally reaches a stopping point
* substantial work has been completed
* context upkeep is required
* switching machines is implied or stated

The summary must always follow the structure below.

---

# 1. Summary Header

Architect must begin each continuity summary with:

```
FROM: Architect
TO: User

Session Continuity Summary  
Last Updated: <Date and Time>
```

---

# 2. High-Level Overview

A concise snapshot of where the project stands.

### **2.1 Current State**

* What has been accomplished
* What components or areas were modified
* What the system looks like right now

### **2.2 Major Changes This Session**

* Bullet points summarizing all completed work
* Include file names or module references where useful

---

# 3. Active Tasks

A list of tasks that are **still in progress** or awaiting further action.

For each task:

* Task name
* Agent responsible
* Status
* Notes

Example:

```
- Add parsing logic module
  - Owner: Implementor
  - Status: Pending Review
  - Notes: Reviewer will evaluate next session
```

---

# 4. Completed Tasks

A list of work finished during this session.

Each entry must include:

* Task name
* Completion date
* Brief description of the work

Architect should also copy these entries into `session-memory.md`.

---

# 5. Pending Decisions

Architect lists decisions that are needed to proceed but have not yet been resolved.

Each decision should include:

* Question or choice
* Context
* Options or trade-offs
* What is blocking progress

Architect may propose next steps or request clarification.

---

# 6. Risks & Concerns

If any risks or uncertainties were identified this session, Architect must capture them here, including:

* Technical uncertainties
* Ambiguous requirements
* Architecture concerns
* Areas requiring further investigation

---

# 7. Recommended Next Steps

A structured list of what the user or agents should do next session.

The list must be **action-oriented**, for example:

```
1. Implementor: Address Reviewer’s comments in module_X.
2. Tester: Add regression tests for Y.
3. Architect: Revisit design decision about component_Z.
4. Documentor: Update README with recent changes.
```

Architect should ensure this sequence is logically consistent with the roadmap.

---

# 8. Key Files or Areas to Review

Pointers to areas requiring the user’s or agents’ attention in future sessions.

Examples:

* Files recently changed
* Modules likely impacted by upcoming work
* Directories relevant to pending decisions
* Documentation requiring expansion

Architect must keep this list concise and practical.

---

# 9. Memory Updates

Architect must explicitly state what updates were applied to:

* `session-memory.md`
* `decisions logs`
* `active task lists`
* Any history or archive files

Architect must confirm these updates are now synchronized.

---

# 10. Resume Instructions

Architect must end each summary with clear instructions for resuming work on another machine:

```
To resume work:

1. Open agents/session-memory.md.
2. Read the Continuity Summary and Active Tasks sections.
3. Review the Recommended Next Steps.
4. In your IDE or environment, say:

   Architect:
   Please rehydrate your context from agents/session-memory.md and generate a fresh plan.

END OF MESSAGE
FROM: Architect
TO: User
```

Architect must always include this footer.

---

# 11. Notes

* This template must be used **without modification** every session.
* Architect may add details, but must not remove or rename sections.
* Architect must store all necessary state in persistent files.
* The purpose of this summary is reliability, continuity, and cross-machine recovery.
