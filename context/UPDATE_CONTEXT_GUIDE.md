# How to Update Context Files

This guide explains how to update the context documentation, especially the master context file used for new chat sessions.

---

## 🚀 Quick Update

Simply run:

```bash
python scripts/update_context.py
```

Or say to Cursor: **"update context"**

This will:
1. ✅ Archive any extraneous .md files from root to `md_archive/`
2. ✅ Update the timestamp in `MASTER_CONTEXT_FOR_NEW_CHAT.md`
3. ✅ Verify all context files are present

---

## 📋 What Gets Updated

### Automatic Updates

The script automatically:
- Archives new .md files (keeps only `README.md` in root)
- Updates timestamp in master context file
- Verifies context file presence

### Manual Updates Needed

You should manually update `MASTER_CONTEXT_FOR_NEW_CHAT.md` when:

1. **Major features completed**
   - Add to "Recent Work Completed" section
   - Update "Current System State"

2. **Configuration changes**
   - Update database IDs
   - Add new environment variables
   - Update configuration structure

3. **Architecture decisions**
   - New modules or patterns
   - Significant refactoring
   - New integration points

4. **Roadmap changes**
   - Priorities shift
   - New features added
   - Timelines updated

5. **New databases or schemas**
   - Notion database additions
   - Schema changes
   - New relations

---

## 🔄 Update Workflow

### Regular Updates (Weekly/Monthly)

```bash
# 1. Run the update script
python scripts/update_context.py

# 2. Review master context file
# Edit: context/MASTER_CONTEXT_FOR_NEW_CHAT.md

# 3. Update roadmap if needed
# Edit: context/UNIFIED_STRATEGIC_ROADMAP.md
```

### After Major Changes

1. **Complete the feature/work**
2. **Run update script**: `python scripts/update_context.py`
3. **Manually update** relevant sections in:
   - `context/MASTER_CONTEXT_FOR_NEW_CHAT.md`
   - `context/UNIFIED_STRATEGIC_ROADMAP.md` (if roadmap changed)
4. **Verify** context files are current

---

## 📝 Sections to Update in Master Context

When manually updating `MASTER_CONTEXT_FOR_NEW_CHAT.md`:

### Section 6: Recent Work Completed
Add new completed features:
```markdown
## ✅ Fully Implemented and Working

1. **New Feature Name**
   - Description
   - Key components
   - Status: Complete
```

### Section 7: Current System State
Update implementation status:
```markdown
### ✅ Complete (Production Ready)
- New feature (just completed)

### 🚧 In Progress
- Feature currently being worked on
```

### Section 4: Key Configuration
Add new environment variables or database IDs if needed.

### Section 8: Key File Locations
Add new modules or scripts if significant.

---

## 🗄️ Archive Management

The update script automatically archives:
- Implementation reports (after completion)
- Status updates (after resolved)
- Historical documentation

Files are moved to `md_archive/` with:
- Original filename preserved
- Timestamp appended if duplicate exists

### Manual Archive Review

Periodically review `md_archive/` to:
- Verify important files are preserved
- Remove truly obsolete files (if needed)
- Organize by category (optional)

---

## ✅ Verification Checklist

After updating context, verify:

- [ ] Master context timestamp is current
- [ ] All recent work is documented
- [ ] Current system state is accurate
- [ ] Configuration information is up-to-date
- [ ] File locations are correct
- [ ] Roadmap status is current
- [ ] Extraneous files are archived

---

## 💡 Tips

1. **Update regularly** - Don't let context get stale
2. **Be concise** - Master context should be comprehensive but not overwhelming
3. **Keep it current** - Outdated context is worse than no context
4. **Use the script** - Automate what you can, update manually what matters
5. **Review before major work** - Check context is current before starting new features

---

## 🔧 Customizing Updates

To customize the update script, edit:
- `scripts/update_context.py` - Main update logic
- Add new checks or updates as needed

---

**Remember**: The master context file is critical for new chat sessions. Keep it current!

