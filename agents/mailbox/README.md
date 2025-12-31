# Agent Mailbox

This directory is used for inter-agent communication via the **mailbox protocol**.

## How It Works

- **File existence = pending message** (delete-after-read pattern)
- Architect writes delegations to `{agent}_MAIL.md`
- Agents write responses to `architect_MAIL.md`
- Files are deleted after being read

## Files

| File | Purpose |
|------|---------|
| `architect_MAIL.md` | Agents write responses here |
| `implementor_MAIL.md` | Delegations for Implementor |
| `reviewer_MAIL.md` | Delegations for Reviewer |
| `tester_MAIL.md` | Delegations for Tester |
| `debugger_MAIL.md` | Delegations for Debugger |
| `documentor_MAIL.md` | Delegations for Documentor |
| `automator_MAIL.md` | Delegations for Automator |

## Protocol

### Sending a Delegation (Architect → Agent)
1. Architect writes to `agents/mailbox/{agent}_MAIL.md`
2. Architect says: "Tell {Agent}: !"
3. Chairman switches chat, types: "!"
4. Agent reads file, deletes it, executes task

### Sending a Response (Agent → Architect)
1. Agent writes to `agents/mailbox/architect_MAIL.md`
2. Agent says: "Tell Architect: !"
3. Chairman switches chat, types: "!"
4. Architect reads file, deletes it, proceeds

### Check Mail Behavior
- If file exists → Read content, delete file, process
- If no file → Respond "No pending mail"

## Chairman Visibility

When agents read their mailbox, they MUST display the full delegation content at the END of their response. This allows Chairman to:
- See exactly what instructions were given
- Interject modifications before agent proceeds
- Maintain oversight of inter-agent communication

## Safety Guards

### 1. One Recipient at a Time
- Each agent ONLY reads their own mailbox file: `{agent}_MAIL.md`
- If an agent is asked to check mail but finds a file for a DIFFERENT agent, respond: "ROUTING ERROR: Found {other}_MAIL.md but I am {self}. Chairman should redirect."
- Do NOT read or process messages intended for other agents

### 2. No Overwriting
- Architect MUST check if `{agent}_MAIL.md` exists before writing
- If file already exists: "BLOCKED: {agent}_MAIL.md already exists. Previous delegation not yet received. Tell {Agent}: ! first."
- Do NOT overwrite existing delegation files
- This prevents lost messages

### 3. Error Responses
When agents encounter these errors, respond clearly:
- Routing error: "ROUTING ERROR: [details]"
- Overwrite blocked: "BLOCKED: [details]"

## Note

The `.md` files in this directory (except this README) are gitignored as they are ephemeral operational data.

