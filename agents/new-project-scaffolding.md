# Project Scaffolding Generator

*Architect instructions for creating new project structures on demand*

This document defines how the Architect should generate new project scaffolding when the user requests it.
The goal is to produce clean, conventional, maintainable, and well-documented starter structures for any kind of project.

Scaffolding should always be:

* minimal
* conventional
* easy to navigate
* well-organized
* ready for expansion

---

# 1. Scaffolding Command (User → Architect)

The user may initiate scaffolding with commands like:

```
Architect:
Create a new project scaffold for <language> that accomplishes <goal>.
```

or:

```
Architect:
Set up a clean, conventional starter structure for a <language> project.
```

or:

```
Architect:
Generate a standard project layout with tests, scripts, and documentation placeholders.
```

Architect must interpret these commands and produce a scaffolding plan and structure.

---

# 2. Scaffolding Process Overview

Whenever the user requests new scaffolding, Architect must:

### **Step 1 — Clarify Inputs (If Needed)**

Architect may ask the user for additional details if the request lacks clarity:

* Target language
* Rough purpose of the project
* Whether tests should be included
* Whether automation should be included

Architect should never guess large structural decisions without confirming.

### **Step 2 — Create a File/Directory Plan**

Architect must create a **structured, numbered plan** including:

* Directory layout
* File placeholders
* Initial scripts or modules
* Where tests will go
* Where docs will go
* Any metadata files needed (README, setup, config, etc.)

### **Step 3 — Delegate to Implementor**

Architect sends tasks to Implementor to generate:

* File structure
* Initial code templates
* Basic imports or boilerplate
* Test skeletons
* Script templates

### **Step 4 — Optional Enhancements**

If requested by the user or generally beneficial:

* Reviewer evaluates structure
* Tester creates test stubs
* Documentor generates initial README or docs
* Automator creates setup/run workflows

### **Step 5 — Architect Finalizes**

Architect collects all outputs and summarizes:

* What was created
* How it's organized
* Where to begin

All summarized into session-memory.md if a new project began.

---

# 3. Standard Scaffolding Layouts

The Architect should default to the following generic structures unless the user specifies otherwise.

## **3.1 Python Project**

```
project/
    src/
        __init__.py
        main.py
    tests/
        __init__.py
        test_main.py
    scripts/
        run.sh
    docs/
        README.md
    requirements.txt
```

## **3.2 JavaScript / TypeScript Project**

```
project/
    src/
        index.js (or index.ts)
    tests/
        index.test.js
    scripts/
        dev.sh
    docs/
        README.md
    package.json
```

## **3.3 Shell Utility Project**

```
project/
    bin/
        tool.sh
    docs/
        README.md
    tests/
        test_tool.sh
```

## **3.4 Mixed / Language-Agnostic**

```
project/
    src/
    tests/
    scripts/
    docs/
    README.md
```

Architect should adapt structures as needed.

---

# 4. Naming & Structure Conventions

* Directory names must be lowercase with hyphens or underscores.
* File names should reflect content purpose clearly.
* Test files should mirror source file structure.
* Scripts should be simple and standalone.
* Documentation should begin with a basic usage guide.

Architect should avoid exotic or unusual structures unless explicitly requested.

---

# 5. README Generation

When scaffolding a new project, the Architect may instruct Documentor to generate a README with:

* Project overview
* How to run
* How to test
* How to extend
* Structure explanation
* Notes about conventions

README formatting must be clean Markdown.

---

# 6. Roadmap Initialization

When creating a new project, Architect should update `session-memory.md` with:

* Initial project purpose
* First tasks
* Known constraints
* Future ideas or extensions

---

# 7. Example Workflow (Architect Internal Steps)

### User:

```
Architect:
Create a Python project scaffold for a tool that processes text files.
```

### Architect:

1. Ask for clarifications (if needed)
2. Generate a directory and file plan
3. Delegate code generation to Implementor
4. Assign documentation creation to Documentor
5. Request optional test skeletons from Tester
6. Request workflow automation from Automator
7. Summarize results and update memory

Architect must follow this pattern consistently unless the user overrides it.

---

# 8. Notes

* Architect must never overwrite existing files unless the user explicitly requests it.
* Scaffolding should always be *minimal but complete*.
* Architect may propose optional enhancements (venv, setup.py, CI, Makefile) based on context, but must not assume user approval.
* This file may be extended later with domain-specific scaffolding patterns.

---

# 9. Purpose of this Template

This ensures the Architect can:

* Initialize new projects cleanly
* Provide predictable structure
* Maintain conventions
* Coordinate agents effectively
* Produce maintainable boilerplate

Scaffolding is a core capability for quickly starting or reorganizing projects.
