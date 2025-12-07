# Technology Stack & Conventions

*Reference for Architect and all agents*

This document defines the preferred technology stack, coding conventions, structural expectations, and general development preferences that the Architect and all agents should follow.

It is intentionally **project-agnostic** and may be expanded as new preferences emerge.

---

# 1. Preferred Languages & Frameworks

### **Primary Language**

* Python (default for general-purpose development unless otherwise specified)

### **Secondary Languages** *(if needed)*

* JavaScript / TypeScript
* Bash / Shell scripting
* Markdown for documentation

### **When to Select a Language**

* Architect should choose Python unless the user explicitly requests another language or the task clearly requires a different ecosystem.

---

# 2. Directory & File Structure Expectations

### **General Structure**

* `src/` for application code
* `tests/` for test files
* `scripts/` for automation or CLI tools
* `docs/` for documentation (if applicable)
* `agents/` for multi-agent coordination files (this folder)

### **Architect Responsibilities**

* Maintain structural consistency
* Ensure Implementor places new files in appropriate directories
* Keep the project organized and predictable

---

# 3. Coding Conventions

### **Python Conventions**

* PEP 8 style unless otherwise stated
* Use type hints consistently
* Use clear, descriptive function and variable names
* Avoid large monolithic functions; prefer modularity
* Prefer dataclasses or typed objects as appropriate

### **Testing Conventions**

* Use pytest for all test suites
* Place unit tests under `tests/unit/`
* Place integration tests under `tests/integration/`
* Use descriptive test names and clear assertions

### **Documentation Conventions**

* Use Markdown for all documentation artifacts
* Include short usage examples where appropriate
* Use headings, lists, and code blocks for clarity

---

# 4. Workflow Conventions

### **Reviews**

* All non-trivial code must pass through Reviewer before considered complete
* Reviewer identifies issues but does not implement fixes

### **Testing**

* Tester must validate any new or modified functionality
* Tests must be deterministic and isolated

### **Automation**

* Automator creates reusable workflows using bash, Python scripts, or Makefiles
* Workflows should be documented by Documentor if user-facing

### **Documentation**

* Documentor provides:

  * High-level overviews
  * Internal architecture explanations
  * Usage instructions
  * Notes on important decisions

---

# 5. Version Control Preferences

Even if the project does not explicitly reference Git, assume:

* Meaningful commit messages (if Architect is asked to generate them)
* Logical commit grouping (per feature or fix)
* Avoid committing generated files unless necessary

Architect may propose commit groupings when appropriate.

---

# 6. Error Handling & Logging Preferences

### **Error Handling**

* Fail loudly when it helps debugging
* Provide clear error messages indicating:

  * cause
  * location
  * expected behavior

### **Logging**

* Minimal logging in small scripts
* Structured logging for major workflows

---

# 7. Environment & Dependency Conventions

### **Python Environments**

* Use a virtual environment or environment manager
* Keep dependency lists minimal and necessary
* Use `requirements.txt` or equivalent if requested

### **CLI Tools**

* Scripts should accept flags rather than hard-coded paths
* Use environment variables where appropriate

---

# 8. Decision-Making Preferences

When the Architect faces multiple design options, it should evaluate based on:

* Simplicity
* Maintainability
* Explicitness
* Readability
* Modularity
* Long-term clarity
* Ease of testing
* Familiarity to the user

Architect should record major decisions in `session-memory.md`.

---

# 9. Extensions (Optional)

You may later add:

* Preferred libraries
* Architectural patterns
* Database conventions
* API conventions
* Infrastructure preferences

Architect should treat this document as **dynamic** and update references as the project evolves.

---

# 10. Summary

This file is a foundational reference defining *how* the system should behave technically.
Architect should:

* Use it when designing plans
* Ensure Implementor respects conventions
* Ensure Tester and Reviewer follow appropriate patterns
* Ensure Documentor writes consistent, structured documentation

It is not a list of restrictions but a **guide for coherent, predictable development** across sessions and machines.
