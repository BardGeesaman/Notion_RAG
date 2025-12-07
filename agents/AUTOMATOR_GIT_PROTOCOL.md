# Automator Agent - Git Commit & Push Protocol

**From:** Reviewer Agent  
**To:** Automator Agent  
**Date:** December 7, 2025  
**Subject:** Required Git Protocol for All Automation Tasks

---

## Background

On December 7, 2025, we experienced significant issues due to:
1. **50+ files modified without commits** - no checkpoints to roll back to
2. **Changes lost when things broke** - had to reset entire codebase
3. **Cross-workstation sync failures** - work done locally wasn't pushed

This document establishes the **required protocol** for all automation tasks.

---

## Required Git Workflow

### 1. Before Starting Any Task

```bash
# Always start with a clean, up-to-date state
cd /Users/bardgeesaman/Documents/RAG
git fetch origin
git status  # Verify clean working directory
git pull origin main  # Get latest changes
```

**If there are uncommitted changes:**
```bash
git stash  # Save them temporarily
# OR
git add -A && git commit -m "WIP: uncommitted changes before [task]"
```

### 2. During Task Execution

**Commit after EVERY successful batch of changes:**

```bash
# After each batch (2-5 files max)
git add -A
git commit -m "Batch X: [description of what changed]"
```

**Commit message format:**
```
[Component] Brief description

- Detail 1
- Detail 2
- Detail 3
```

**Examples:**
```bash
git commit -m "Fix MW search: add analytical_platform parameter"
git commit -m "Batch 1: Remove Notion imports from 5 ingestion files"
git commit -m "Add external_ids parameter to create_or_update_dataset_in_postgres"
```

### 3. Before Major Changes

**Create a checkpoint:**
```bash
git add -A
git commit -m "CHECKPOINT: Before [major change description]"
git push origin main
```

This ensures you can always roll back:
```bash
git reset --hard HEAD~1  # Undo last commit
# OR
git reset --hard <commit-hash>  # Go to specific commit
```

### 4. After Completing Task

**Always push to GitHub:**
```bash
git add -A
git commit -m "[Task] Complete: [summary]"
git push origin main
```

### 5. End of Session

**Final push with summary:**
```bash
git add -A
git commit -m "End of session: [date] - [summary of all changes]"
git push origin main
```

---

## Automation Script Template

For any automated task, use this template:

```python
import subprocess
import sys

def git_checkpoint(message: str) -> bool:
    """Create a git checkpoint with commit and push."""
    try:
        # Stage all changes
        subprocess.run(["git", "add", "-A"], check=True, cwd=PROJECT_ROOT)
        
        # Check if there are changes to commit
        result = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True, text=True, cwd=PROJECT_ROOT
        )
        
        if not result.stdout.strip():
            print("No changes to commit")
            return True
        
        # Commit
        subprocess.run(
            ["git", "commit", "-m", message],
            check=True, cwd=PROJECT_ROOT
        )
        
        # Push
        subprocess.run(
            ["git", "push", "origin", "main"],
            check=True, cwd=PROJECT_ROOT
        )
        
        print(f"✅ Checkpoint: {message}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Git error: {e}")
        return False


def run_batch_task():
    """Example batch task with proper checkpointing."""
    
    # Checkpoint before starting
    git_checkpoint("CHECKPOINT: Before batch task")
    
    for batch_num, batch in enumerate(batches):
        try:
            # Do the work
            process_batch(batch)
            
            # Test that it works
            test_imports()
            
            # Checkpoint after each batch
            git_checkpoint(f"Batch {batch_num + 1}: {batch.description}")
            
        except Exception as e:
            print(f"❌ Batch {batch_num + 1} failed: {e}")
            print("Rolling back to last checkpoint...")
            subprocess.run(["git", "checkout", "--", "."], cwd=PROJECT_ROOT)
            sys.exit(1)
    
    # Final checkpoint
    git_checkpoint("Task complete: all batches processed")
```

---

## Frequency Guidelines

| Task Type | Commit Frequency |
|-----------|------------------|
| File modifications | Every 2-5 files |
| Bug fixes | After each fix |
| Feature additions | After each working feature |
| Refactoring | After each logical unit |
| Batch operations | After each batch |
| Long-running tasks | Every 5-10 minutes of work |

---

## What NOT To Do

❌ **Don't modify 50+ files without committing**
❌ **Don't end a session without pushing**
❌ **Don't skip testing before committing**
❌ **Don't use vague commit messages like "updates" or "fixes"**
❌ **Don't assume the other workstation has your changes**

---

## Recovery Commands

If something goes wrong:

```bash
# See what changed
git diff

# Undo uncommitted changes
git checkout -- .

# Undo last commit (keep changes)
git reset --soft HEAD~1

# Undo last commit (discard changes)
git reset --hard HEAD~1

# Go back to specific commit
git log --oneline -10  # Find the commit hash
git reset --hard <commit-hash>

# Get back to remote state
git fetch origin
git reset --hard origin/main
```

---

## Cross-Workstation Sync

When switching workstations:

**On OLD workstation (before leaving):**
```bash
git add -A
git commit -m "WIP: switching workstations"
git push origin main
```

**On NEW workstation (before starting):**
```bash
cd /Users/bardgeesaman/Documents/RAG
git fetch origin
git pull origin main
```

---

## Summary

1. **Pull before starting**
2. **Commit after every batch**
3. **Push after every task**
4. **Push before leaving**

**The goal:** Anyone should be able to `git pull` and have a working codebase at any time.

---

*Protocol established: December 7, 2025*
*Based on lessons learned from codebase sync issues*

