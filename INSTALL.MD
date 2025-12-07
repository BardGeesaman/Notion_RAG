# INSTALL.md  
*How to install and initialize the 6-agent system*

This guide explains how to install, configure, and activate the multi-agent system on any workstation using Cursor.

---

# 1. Requirements

- Cursor (latest version)  
- A project folder containing:  
  ```
  /agents
      architect.md
      implementor.md
      reviewer.md
      tester.md
      automator.md
      documentor.md
      ceo.md
      agent-team-charter.md
      session-memory.md
      reset-context.md
      modes-of-operation.md
      system-modes.md
      tech-stack.md
      glossary.md
      agent-commands-reference.md
      continuity-summary-template.md
      new-project-scaffolding.md
      agent-replication.md
  ```

These define roles, behaviors, workflows, and continuity rules.

---

# 2. Create the Six Agents in Cursor

Cursor v0.47+ uses **tabs as agents**.

For each agent:

1. Open a new chat  
2. Rename it: Architect / Implementor / Reviewer / Tester / Automator / Documentor  
3. Paste the content of the matching `/agents/<agent>.md` file into the chat as the first message  

This defines the agent’s role and protocol for the session.

---

# 3. Tell Architect to Load Definitions

After all agents are created, say:

```
Architect:
All agent instruction files under /agents have been loaded. 
Please treat them, along with the Agent Team Charter, as the authoritative definitions for all roles and workflows.
Confirm readiness.
```

---

# 4. Load Project State

If this is not a fresh project, run:

```
Architect:
Rehydrate your context from agents/session-memory.md and generate a fresh plan.
```

---

# 5. Validate System Health

Run the validation suite:

```
Architect:
Run a full multi-agent system validation.
```

This verifies:

- Roles  
- Delegation  
- Message protocol  
- Cross-agent workflow  
- Continuity summary discipline  

---

# 6. You Are Now Ready to Work

You interact **only with Architect**.

Commands begin with:

```
Architect:
```

You may now build, test, refactor, extend, automate, or document anything using the six-agent system.

```
END OF FILE
```
