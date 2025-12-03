# Context Update System

**Purpose**: Automated system to keep context files current and archive extraneous documentation.

---

## 🚀 Quick Start

### Update Context Command

Simply say to Cursor:
```
update context
```

Or run manually:
```bash
python scripts/update_context.py
```

---

## ✅ What It Does

The update script automatically:

1. **Archives Extraneous Files**
   - Moves any new .md files from root to `md_archive/`
   - Keeps only `README.md` in root
   - Preserves files if duplicates exist (appends timestamp)

2. **Updates Timestamps**
   - Updates "Last Updated" in `MASTER_CONTEXT_FOR_NEW_CHAT.md`
   - Sets current date

3. **Verifies Context Files**
   - Checks all expected context files exist
   - Reports any missing files

---

## 📋 When to Run

### Automatic Triggers

Run "update context" when:
- ✅ After completing major features
- ✅ After significant changes
- ✅ Periodically (weekly/monthly)
- ✅ Before starting new work
- ✅ When root directory gets cluttered

### After Major Changes

1. Complete the feature/work
2. Say "update context" or run the script
3. Manually update `MASTER_CONTEXT_FOR_NEW_CHAT.md` sections:
   - Recent Work Completed
   - Current System State
   - Any configuration changes

---

## 📝 Manual Updates Still Needed

The script handles:
- ✅ File archiving
- ✅ Timestamp updates
- ✅ File verification

You still need to manually update:
- 📝 Recent work descriptions
- 📝 Current system state
- 📝 Configuration details
- 📝 Roadmap status

See `context/UPDATE_CONTEXT_GUIDE.md` for detailed manual update instructions.

---

## 🔧 Script Location

- **Script**: `scripts/update_context.py`
- **Executable**: Yes (can run directly)
- **Dependencies**: Python 3 (standard library only)

---

## 📂 Files Managed

### Context Files (in `context/`)
- `MASTER_CONTEXT_FOR_NEW_CHAT.md` - Main context (timestamp updated)
- `NEW_CHAT_QUICK_START.md`
- `HIGH_LEVEL_STRATEGY_CONTEXT.md`
- `UNIFIED_STRATEGIC_ROADMAP.md`
- All other context files (verified)

### Archive Files (in `md_archive/`)
- Implementation reports
- Status updates
- Historical documentation
- All .md files except `README.md`

### Root Directory
- Only `README.md` remains
- All other .md files archived automatically

---

## 💡 Usage Tips

1. **Run regularly** - Don't let context get stale
2. **Review after running** - Check what was archived
3. **Manual updates** - Still edit master context for content
4. **Say "update context"** - Easy command for Cursor

---

## 🔄 Workflow Example

```bash
# 1. Complete a major feature
# ... do work ...

# 2. Update context (automated)
python scripts/update_context.py
# or say: "update context"

# 3. Manual review and update
# Edit: context/MASTER_CONTEXT_FOR_NEW_CHAT.md
# - Add to "Recent Work Completed"
# - Update "Current System State"

# 4. Verify
# Check context files are current
# Review archived files
```

---

**The context update system keeps your documentation organized and current!**

