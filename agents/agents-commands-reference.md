# Agent Commands Reference

*How the User Communicates with the Architect*

This reference defines the standard commands the user may issue to the Architect.
All commands must begin with:

```
Architect:
```

Architect will interpret the directive, plan accordingly, and coordinate the full agent team.

---

# 1. Core Operational Commands

### **Start a new task**

```
Architect:
Please begin a new task: <describe goal>.
```

### **Break a goal into a plan**

```
Architect:
Create a step-by-step plan for <describe goal>.
```

### **Request a full multi-agent workflow**

```
Architect:
Coordinate a full workflow to accomplish <goal>, including planning, implementation, review, testing, automation, and documentation.
```

### **Ask for clarification / rephrasing**

```
Architect:
Please restate my request more clearly and confirm understanding.
```

---

# 2. Status & Roadmap Commands

### **Show current status**

```
Architect:
Give me a status update.
```

### **Show the roadmap**

```
Architect:
Display the current roadmap and active tasks.
```

### **Summarize the session so far**

```
Architect:
Provide a summary of all work done in this session.
```

### **Summarize the entire project**

```
Architect:
Generate a high-level project summary.
```

### **What are we waiting on?**

```
Architect:
List all blockers, pending reviews, or unresolved decisions.
```

---

# 3. Implementation Commands

### **Implement something**

```
Architect:
Implement <feature or change>. Handle all subtasks.
```

### **Modify a specific file or module**

```
Architect:
Make the following changes to <file or module>:
<description>
```

### **Refactor**

```
Architect:
Refactor <area> for clarity and maintainability.
```

---

# 4. Review Commands

### **Request a review**

```
Architect:
Have the Reviewer evaluate the recent implementation for correctness and quality.
```

### **Ask for targeted review**

```
Architect:
Review <specific part> for <concern>.
```

---

# 5. Testing Commands

### **Request a test plan**

```
Architect:
Have the Tester create a test plan for <area>.
```

### **Request full test suite**

```
Architect:
Have the Tester design and/or update the test suite for this module.
```

### **Execute tests**

```
Architect:
Ask the Tester to run the tests and report results.
```

---

# 6. Documentation Commands

### **Request documentation**

```
Architect:
Have the Documentor create documentation for <topic or component>.
```

### **Explain something**

```
Architect:
Explain <topic> in a clear and structured way.
```

### **Generate a README**

```
Architect:
Have the Documentor create a README for <component or project>.
```

---

# 7. Automation Commands

### **Request a workflow**

```
Architect:
Have the Automator create a workflow that automates <process>.
```

### **Request operational tooling**

```
Architect:
Have Automator build a script or procedure for <task>.
```

---

# 8. Mode Commands

### **Strict Mode**

```
Architect:
enter strict mode
```

### **Loose Mode**

```
Architect:
enter loose mode
```

### **Check current mode**

```
Architect:
report current mode
```

### **Reset context**

```
Architect:
reset context but preserve all persistent memory and rebuild context.
```

---

# 9. Continuity Commands

### **End the session**

```
Architect:
I am done for now. Please generate the continuity summary.
```

### **Prepare for switching machines**

```
Architect:
Please update session-memory.md and provide a continuity summary.
```

### **Resume from memory**

```
Architect:
Please rehydrate your context from agents/session-memory.md and generate a fresh plan.
```

---

# 10. Safety & Clarification Commands

### **Stop and reassess**

```
Architect:
Pause all ongoing tasks and reassess the plan.
```

### **Request alternative options**

```
Architect:
Provide alternative approaches for <goal>.
```

### **Identify risks**

```
Architect:
List risks or concerns with the current plan.
```

### **Re-plan**

```
Architect:
Revise the plan based on the following changes: <changes>.
```

---

# 11. Meta-Commands (System-Level)

### **Optimize workflow**

```
Architect:
Based on our history, suggest workflow improvements.
```

### **Reflect on current structure**

```
Architect:
Evaluate whether our architecture or conventions need updating.
```

### **Evaluate agent roles**

```
Architect:
Assess whether any agent roles or instructions need revision.
```

---

# 12. Notes

* These commands are **guidelines**, not strict requirements.
* Architect must interpret user intent even if commands are phrased in natural language.
* All commands trigger Architect → multi-agent workflows, never direct agent calls.
* This file may grow over time with additional shorthand or domain-specific commands.
