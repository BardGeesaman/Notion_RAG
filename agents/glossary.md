# Glossary of Terms

*Shared terminology reference for the Architect and all agents.*

This document defines terms that the Architect and all agents may encounter.
It ensures consistent interpretation across sessions, machines, and contexts.

Terms here may be generic or user-defined.
You may extend or modify this list at any time.

---

## 1. General Terms

### **Task**

A unit of work assigned by the Architect to an agent.

### **Plan**

A structured, numbered sequence of steps created by the Architect before delegation.

### **Roadmap**

A living, high-level list of current tasks, milestones, and future work maintained by the Architect.

### **Session**

A continuous working period; ends when the Architect generates a continuity summary.

### **Continuity Summary**

A structured description of:

* current state
* completed work
* pending items
* next steps
  Generated automatically at session end.

### **Context Reset**

A controlled refresh of the Architect’s working memory.
Defined fully in `reset-context.md`.

---

## 2. Agent-Specific Terms

### **Architect**

Master planner, coordinator, and source of truth.
Owns delegation, decisions, roadmap, and summaries.

### **Implementor**

Writes and modifies code exactly as specified.
Never plans, reviews, tests, or documents.

### **Reviewer**

Evaluates code for correctness, clarity, and risk.
Does not modify code.

### **Tester**

Designs and writes tests.
Analyzes failures and reports findings.

### **Automator**

Creates repeatable workflows, scripts, or procedural sequences.

### **Documentor**

Writes clear, structured documentation and explanations.

---

## 3. Work Product Terms

### **Code Diff**

A representation of code changes.
Implementor uses it to communicate modifications.

### **Review Report**

A structured analysis from the Reviewer with strengths, issues, and recommendations.

### **Test Plan**

A description of tests the Tester proposes or will create.

### **Test Report**

A structured output summarizing test execution and results.

### **Workflow**

A sequence of steps created by Automator representing a repeatable process.

### **Documentation Artifact**

Any written explanation, guide, README, or conceptual overview produced by Documentor.

---

## 4. Status Terms

### **In Progress**

Work currently being performed by an agent.

### **Blocked**

Work that cannot continue due to missing information or unresolved dependencies.

### **Pending Review**

Awaiting Reviewer evaluation.

### **Pending Tests**

Awaiting Tester validation.

### **Ready for Documentation**

Implementation and review completed; Documentor should finalize user-facing explanation.

### **Completed**

Work fully finished and recorded in session memory.

---

## 5. Decision Terms

### **Design Decision**

A choice made by the Architect affecting architecture, structure, or approach.

### **Constraint**

A restriction on implementations (performance, style, dependencies, safety, etc.).

### **Trade-Off**

A documented reasoning between multiple possible directions.

### **Revisit Flag**

Indication that a decision should be reconsidered later.

---

## 6. User Communication Terms

### **Directive**

A command from the user beginning with “Architect: ...”.

### **Clarification Request**

A message from an agent to the Architect asking for missing context or specificity.

### **Mode Request**

A command changing Architect's operational mode (Strict, Loose, Reset).

---

## 7. Optional Project-Specific Extensions

You may add new entries (terms, definitions, categories) at any time, such as:

* Domain terminology
* Common file or module names
* Architectural patterns
* Preferred conventions
* Reusable schemas or workflows
* Repeated problem types

The Architect should use this glossary as a **source of truth** for term interpretation.
