# Dashboard Import Pattern - Preventing Streamlit Caching Issues

**Date:** 2025-12-05  
**Issue:** Streamlit module caching can cause `AttributeError` when importing models directly  
**Solution:** Use direct module import pattern

---

## Problem

When importing SQLAlchemy models in Streamlit dashboard pages, you may encounter:
```
AttributeError: module 'amprenta_rag.database.models' has no attribute 'Compound'
```

This happens because:
- Streamlit caches modules aggressively
- Direct imports (`from module import Class`) can fail if the module isn't fully loaded
- SQLAlchemy models need the module to be fully initialized

---

## Solution: Direct Module Import Pattern

**✅ CORRECT Pattern (use this):**

```python
# Import the models module directly
import amprenta_rag.database.models as db_models

# Access models via module attribute
Compound = db_models.Compound
HTSCampaign = db_models.HTSCampaign
BiochemicalResult = db_models.BiochemicalResult
```

**❌ AVOID This Pattern:**

```python
# This can fail with Streamlit caching
from amprenta_rag.database.models import (
    Compound,
    HTSCampaign,
    BiochemicalResult,
)
```

---

## When to Use This Pattern

Use the direct module import pattern for:

1. **SQLAlchemy models** in dashboard pages
2. **Any module that might not be fully loaded** when Streamlit caches
3. **Modules with complex initialization** (like database models)

**Examples:**
- `amprenta_rag.database.models` - All SQLAlchemy models
- Any module with `__init__.py` that does complex setup

---

## Standard Pattern for Dashboard Pages

Here's the standard import pattern for new dashboard pages:

```python
"""Your dashboard page."""

from __future__ import annotations

import pandas as pd
import streamlit as st

# For SQLAlchemy models - use direct module import
import amprenta_rag.database.models as db_models
YourModel = db_models.YourModel

# For other imports - direct import is fine
from scripts.dashboard.db_session import db_session
from amprenta_rag.some_other_module import some_function
```

---

## Current Implementation

**File:** `scripts/dashboard/pages/chemistry.py`

**Pattern Used:**
```python
import amprenta_rag.database.models as db_models

Compound = db_models.Compound
HTSCampaign = db_models.HTSCampaign
BiochemicalResult = db_models.BiochemicalResult
```

---

## Testing the Pattern

To verify imports work:

```bash
python -c "import sys; sys.path.insert(0, '.'); from scripts.dashboard.pages.chemistry import Compound; print('✅ Import works')"
```

---

## Future Dashboard Pages

When creating new dashboard pages that use database models:

1. **Always use** the direct module import pattern for models
2. **Test imports** before running Streamlit
3. **Document** if you find other modules that need this pattern

---

## Related Files

- `scripts/dashboard/pages/chemistry.py` - Example implementation
- `scripts/dashboard/pages/overview.py` - Uses direct imports (works fine)
- `scripts/dashboard/pages/datasets.py` - Uses direct imports (works fine)

---

**Last Updated:** 2025-12-05
