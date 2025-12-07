# Agent Replication Guide  
*How to recreate the complete six-agent system on any workstation*

This guide ensures the Architect, Implementor, Reviewer, Tester, Automator, and Documentor agents can be recreated exactly — including their instructions, behaviors, workflows, and recommended model selections.

---

# 1. Requirements

Any workstation must have:

- Cursor (latest version)
- Access to your project repository with the `/agents` folder
- Availability of the recommended model families (e.g., GPT-5.1, Claude 4.5, Gemini 3 Pro, etc.)

Cursor does not bind models to prompts — you select the model when opening the agent chat.

---

# 2. Folder Structure (Must Exist on Every Machine)

Ensure this directory structure exists:

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
    agent-replication.md  ← (this file)
```

Names may vary slightly (e.g., `architect-instructions.md` vs `architect.md`), but each agent must have one authoritative file.

---

# 3. Recreating Each Agent in Cursor

Cursor’s newer versions do not have a UI for “Agent Settings”; tabs *are* the agents.

To recreate agents:

1. Open a new chat tab  
2. Rename the tab to the correct agent name:
   - Architect  
   - Implementor  
   - Reviewer  
   - Tester  
   - Automator  
   - Documentor  
3. Paste into the chat:
   - The **Universal Header**, and  
   - The contents of `agents/<agent>.md`

Repeat for all six.

This becomes the agent’s system prompt for the session.

---

# 4. Architect Initialization

After all agents are recreated, send this to Architect:

```
Architect:
All agent instruction files under /agents have been loaded. 
Please treat them, along with the Agent Team Charter, as the authoritative definitions 
for all roles, constraints, workflows, and communication formats going forward.
Confirm recognition of all six agents and system readiness.
```

Architect will:

- Load the Charter  
- Recognize all agent roles  
- Realign its internal model  
- Report readiness  

---

# 5. Recommended Model Assignments

You will manually select models in each agent tab.

Suggested defaults:

| Agent        | Recommended Model                    |
|--------------|--------------------------------------|
| Architect    | Claude 4.5 Sonnet **or** GPT-5.1     |
| Implementor  | GPT-5.1 or GPT-5.1 Codex Max         |
| Reviewer     | Gemini 3 Pro or Claude 4.5 Sonnet    |
| Tester       | GPT-5.1 Codex Max                    |
| Automator    | GPT-5.1 or GPT-5.1 Codex Max         |
| Documentor   | Claude 4.5 Sonnet                    |

This ensures stable, predictable multi-agent behavior.

---

# 6. Restoring Context & Project State

To load project state on a new workstation, run:

```
Architect:
Please rehydrate your context from agents/session-memory.md and generate a fresh plan for next steps.
```

Architect will reconstruct:

- Active tasks  
- Completed work  
- Architecture  
- Decisions  
- Roadmap  

Everything your previous workstation had.

---

# 7. System Verification After Replication

Run a full validation to ensure the agents work correctly:

```
Architect:
Run a full multi-agent system validation.
```

Architect will:

- Verify role alignment  
- Verify messaging protocol  
- Test delegation  
- Evaluate review → test → documentation workflow  
- Confirm system health  


---

# 8. Updating Agents in the Future

When you modify any `.md` files inside `/agents`:

1. Commit changes to your repo  
2. Re-paste updated content into each agent’s chat tab  
3. Tell Architect:

```
Architect:
Re-anchor yourself and all agents to the updated instruction files.
```

This prevents drift and preserves consistency.

---

# 9. Summary

To recreate this six-agent system anywhere:

1. Clone the repository  
2. Ensure `/agents` folder is present  
3. Open six Cursor tabs (one per agent)  
4. Paste each agent’s instructions into its own chat  
5. Tell Architect to load the agent definitions  
6. Rehydrate context  
7. Validate the system  

Your multi-agent environment will now be identical across all machines.

```
END OF FILE
```
