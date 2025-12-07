# /agents Directory  
*Authoritative definitions for the 6-agent multi-agent system*

This directory contains all role definitions, behavioral specifications, templates, workflows, and continuity tools required to run the multi-agent development architecture.

The **Architect** loads these files to coordinate all other agents:
- Implementor  
- Reviewer  
- Tester  
- Automator  
- Documentor  

---

# 1. Purpose of This Folder

This folder provides:

- Stable definitions for all agents  
- Project-agnostic workflow rules  
- Session continuity management  
- Context rehydration  
- Mode management (Strict, Loose, Reset)  
- Drift correction rules  
- Replication instructions for new machines  

It is the **single source of truth** for agent behavior.

---

# 2. Contents Overview

| File | Purpose |
|------|---------|
| `architect.md` | Role + behavior of Architect |
| `implementor.md` | Role + behavior of Implementor |
| `reviewer.md` | Role + behavior of Reviewer |
| `tester.md` | Role + behavior of Tester |
| `automator.md` | Role + behavior of Automator |
| `documentor.md` | Role + behavior of Documentor |
| `ceo.md` | Instructions for the human operator |
| `agent-team-charter.md` | Full definition of roles + workflows |
| `session-memory.md` | Persistent memory across machines |
| `continuity-summary-template.md` | Template for end-of-session summaries |
| `reset-context.md` | How Architect performs controlled resets |
| `modes-of-operation.md` | Strict / Loose mode overview |
| `system-modes.md` | Internal mechanics of modes |
| `tech-stack.md` | Language, formatting, and code conventions |
| `glossary.md` | System terminology definitions |
| `agent-commands-reference.md` | Commands the user can send to Architect |
| `new-project-scaffolding.md` | How Architect scaffolds new projects |
| `agent-replication.md` | Instructions for recreating agents on new machines |
| `replication-checklist.md` | Fast startup checklist |

---

# 3. How These Files Are Used

- **Architect** directly reads and relies on these files.  
- **Other agents** follow instructions passed through Architect.  
- **You** (CEO) maintain these for clarity and control.  

Whenever behavior feels “off”, run:

```
Architect:
Perform a drift check and realign with all files in the /agents directory.
```

---

# 4. How to Update Agent Instructions

1. Update the appropriate `*.md` file in this directory  
2. Paste the updated content into that agent’s Cursor chat  
3. Tell Architect:  

```
Architect:
Re-anchor yourself and all agents to the updated instruction files.
```

This preserves consistency across machines and sessions.

---

# 5. Replicating Agents on a New Machine

See `agent-replication.md` or run this quick checklist:

```
1. Clone repo
2. Create 6 agent tabs in Cursor
3. Paste each agent's .md file into its tab
4. Tell Architect to load /agents
5. Rehydrate context
```

Your multi-agent environment is now fully portable and deterministic.

---

If you add new files here (like templates or workflows), notify Architect so it can use them appropriately.
