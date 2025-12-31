# Agent Mailbox

This directory is used for inter-agent communication via the **mailbox protocol**.

## How It Works

- **File existence = pending message** (delete-after-read pattern)
- Architect writes delegations to `{agent}.md`
- Agents write responses to `architect.md`
- Files are deleted after being read

## Files

| File | Purpose |
|------|---------|
| `architect.md` | Agents write responses here |
| `implementor.md` | Delegations for Implementor |
| `reviewer.md` | Delegations for Reviewer |
| `tester.md` | Delegations for Tester |
| `debugger.md` | Delegations for Debugger |
| `documentor.md` | Delegations for Documentor |
| `automator.md` | Delegations for Automator |

## Protocol

### Sending a Delegation (Architect → Agent)
1. Architect writes to `agents/mailbox/{agent}.md`
2. Architect says: "Tell {Agent} to check mail"
3. Chairman switches chat, says: "Check mail"
4. Agent reads file, deletes it, executes task

### Sending a Response (Agent → Architect)
1. Agent writes to `agents/mailbox/architect.md`
2. Agent says: "Tell Architect to check mail"
3. Chairman switches chat, says: "Check mail"
4. Architect reads file, deletes it, proceeds

### Check Mail Behavior
- If file exists → Read content, delete file, process
- If no file → Respond "No pending mail"

## Note

The `.md` files in this directory (except this README) are gitignored as they are ephemeral operational data.

