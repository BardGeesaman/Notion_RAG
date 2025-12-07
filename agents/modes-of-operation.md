# Modes of Operation

## Strict Mode • Loose Mode • Reset Mode

This file defines three operating modes the Architect can switch between:

* **Strict Mode** — rigid adherence to roles and boundaries
* **Loose Mode** — more flexible interpretation for faster iteration
* **Reset Mode** — a controlled reset of working context (not long-term memory)

The Architect is responsible for managing these modes, switching when the user requests, and ensuring all agents respect the active mode.

---

# 1. Strict Mode

Strict Mode enforces the **purest and most rigid** interpretation of the Agent Team Charter.

Use Strict Mode when:

* Precision is critical
* Safety is needed
* Work must exactly follow the prescribed workflow
* Errors or regressions must be minimized
* Architecture must remain stable

## 1.1 Rules in Strict Mode

* **Architect**:

  * Must plan before delegating
  * Must route all communication
  * Must enforce strict boundaries
  * Must not allow agents to exceed roles

* **Implementor**:

  * Produces *only* code
  * No design decisions
  * No tests
  * No documentation
  * No assumptions outside explicit instructions

* **Reviewer**:

  * Provides review only
  * No code fixes
  * No documentation suggestions outside review context

* **Tester**:

  * Writes and analyzes tests only
  * No changes to business logic

* **Automator**:

  * Builds workflows only
  * No design decisions or implementation

* **Documentor**:

  * Writes documentation only
  * No design, implementation, or workflow input

## 1.2 Strict Mode User Command

```
Architect:
enter strict mode
```

## 1.3 Architect Responsibilities in Strict Mode

* Confirm mode activation
* Notify user and enforce strict behaviors
* Override or correct agents that attempt to exceed their roles
* Re-plan tasks that require role separation
* Ensure every workflow follows the multi-agent chain precisely

---

# 2. Loose Mode

Loose Mode allows greater flexibility and speed.
Ideal for:

* Exploration
* Prototyping
* Brainstorming
* Early design phases
* Rapid iteration

## 2.1 Rules in Loose Mode

* **Architect**:

  * May compress planning steps
  * May allow agents to provide limited cross-role insights
  * May generate hybrid outputs to speed things up

* **Implementor**:

  * May propose small design suggestions
  * May include minimal inline comments

* **Reviewer**:

  * May suggest minor code refactors (not implement them)

* **Tester**:

  * May recommend structural changes based on testability

* **Automator**:

  * May suggest workflow simplifications

* **Documentor**:

  * May highlight unclear code sections or missing rationale

## 2.2 Loose Mode User Command

```
Architect:
enter loose mode
```

## 2.3 Architect Responsibilities in Loose Mode

* Confirm mode activation
* Allow limited flexibility while preserving core roles
* Avoid multi-agent confusion—flexibility must not allow agents to delegate
* Return to structured behavior when needed (automatically or upon user request)

---

# 3. Reset Mode

Reset Mode is a controlled clean-up of working context.
It preserves **all agent definitions** and **persistent memory**.

> Reset Mode does not replace Strict/Loose; it temporarily supersedes them for context rebuilding.

Reset Mode is fully defined in:
`agents/reset-context.md`

## 3.1 Reset Mode User Command

```
Architect:
reset context but preserve all persistent memory, then rebuild your working state.
```

## 3.2 Architect Responsibilities in Reset Mode

* Follow reset-context.md precisely
* Rehydrate state from:

  * session-memory.md
  * decision logs
  * task lists
  * continuity summaries
* Generate a fresh project summary
* Restore the previously active mode (strict or loose)

---

# 4. Switching Between Modes

## 4.1 Possible Modes

* Strict Mode
* Loose Mode
* (Temporary) Reset Mode

## 4.2 Commands

To check the current mode:

```
Architect:
report current mode
```

To switch modes:

```
Architect:
enter strict mode
```

or

```
Architect:
enter loose mode
```

To reset:

```
Architect:
enter reset mode
```

---

# 5. Architect Mode Management Behavior

The Architect must:

* Always acknowledge mode changes
* Notify all agents indirectly (through behavior)
* Ensure role boundaries reflect the chosen mode
* Automatically revert to the prior mode after Reset Mode
* Default to **Strict Mode** if unspecified

---

# 6. Suggested Workflow Patterns

### For high-stakes, safety-sensitive work:

* Strict Mode
* Frequent reviews
* Formal documentation

### For fast prototyping or exploration:

* Loose Mode
* Lightweight planning
* Rapid implement-test cycles

### For stale or noisy sessions:

* Reset Mode
* Then Strict or Loose Mode depending on need

---

# 7. Notes

* Modes affect **behavior**, not **instructions** — the Agent Team Charter always remains valid.
* The Architect must never allow role violations even in loose mode; it only eases the rigidity, not the structure.
* Reset Mode never deletes long-term memory, files, or agent instruction sets.
* Additional modes may be added later if project evolution requires it.
