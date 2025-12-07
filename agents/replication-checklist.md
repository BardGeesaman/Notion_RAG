# Quick Replication Checklist  
*Fast steps to recreate the agent system on any workstation*

1. Clone the repository and ensure the `/agents` directory is present.  
2. Open six new Cursor tabs and rename them: Architect, Implementor, Reviewer, Tester, Automator, Documentor.  
3. Paste the contents of each corresponding `/agents/*.md` file into each tab as the first message.  
4. Tell Architect:  
   ```
   Architect:
   Load all agent instruction files from /agents and confirm full system readiness.
   ```
5. Rehydrate context:  
   ```
   Architect:
   Rehydrate your context from agents/session-memory.md and generate a fresh plan.
   ```
