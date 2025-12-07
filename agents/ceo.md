# CEO Instructions  
*How to operate the multi-agent system effectively*

You interact **only with the Architect**.  
Architect coordinates all work across Implementor, Reviewer, Tester, Automator, and Documentor.

Your responsibility is to describe *what* you want — Architect figures out *how* to accomplish it.

---

# 1. How to Start Any Task

Begin every request with:

```
Architect:
<your request here>
```

Examples:

- “Architect: implement a new feature that does X.”
- “Architect: create a plan for Y.”
- “Architect: refactor module Z for clarity.”
- “Architect: coordinate a full multi-agent workflow for this task.”

Architect will always:
1. Interpret your intent  
2. Generate a plan  
3. Delegate tasks to the proper agents  
4. Summarize results back to you  

---

# 2. How to Ask for Status or Updates

To check what's going on:

```
Architect:
Give me a status update.
```

To see the roadmap:

```
Architect:
Show me the current roadmap and active tasks.
```

To summarize the session so far:

```
Architect:
Summarize all work completed this session.
```

To see what’s blocking progress:

```
Architect:
List all blockers, outstanding reviews, tests, or decisions.
```

---

# 3. Ending a Session

When you’re done for now:

```
Architect:
I am done. Please generate the continuity summary.
```

Architect will:

- Update `agents/session-memory.md`
- Follow `continuity-summary-template.md`
- Produce structured instructions for resuming on any machine

---

# 4. Resuming Work (On Any Machine)

After loading your repo, say:

```
Architect:
Please rehydrate your context from agents/session-memory.md and generate a fresh plan.
```

Architect will rebuild its internal state automatically.

---

# 5. Mode Controls (Strict / Loose / Reset)

You may shift the system’s operating behavior:

### Strict mode:
```
Architect:
enter strict mode
```

### Loose mode:
```
Architect:
enter loose mode
```

### Check current mode:
```
Architect:
report current mode
```

### Reset context safely:
```
Architect:
reset context but preserve all persistent memory and rebuild your internal model.
```

---

# 6. Drift Detection and Correction

If anything feels “off”, use:

```
Architect:
Perform a drift check and realign with the Agent Team Charter and all /agents/*.md files.
```

Architect will:
- Re-read the authoritative files  
- Correct its internal model  
- Report any deviations  

---

# 7. Full System Validation (Health Check)

Run this after major updates or if behavior seems inconsistent:

```
Architect:
Run a full multi-agent system validation.
```

Architect will:
- Verify roles  
- Verify message protocol  
- Verify delegation  
- Verify review/test/documentation workflows  
- Report system health  

---

# 8. Additional Useful Commands

### Show open questions:
```
Architect:
List all unresolved questions or pending decisions.
```

### Get options for a design choice:
```
Architect:
Provide alternative approaches for <topic>.
```

### Ask Architect to reflect on workflow or architecture:
```
Architect:
Evaluate whether our architecture or conventions need improvement.
```

---

# 9. What You Do NOT Need to Do

You do **not**:

- Talk directly to Implementor, Tester, Reviewer, Documentor, Automator  
- Orchestrate multi-agent workflows  
- Manage message protocol  
- Maintain context manually  
- Rewrite agent instructions  
- Remember session state  

Architect handles all internal logic, planning, delegation, coordination, and continuity.

---

# 10. Your Role as CEO

Your job is simple:

- Describe goals  
- Clarify requirements  
- Request summaries or updates  
- Approve or redirect plans  
- Initiate or close sessions  

The system handles everything else.

You are the **vision and direction**.  
Architect is the **execution and coordination**.
