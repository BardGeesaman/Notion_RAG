#!/usr/bin/env python
"""Validate all critical imports work."""
import sys

try:
    from amprenta_rag.api.main import app
    from amprenta_rag.database.models import Program, Experiment, Dataset, Feature, Signature
    from amprenta_rag.models.chemistry import HTSCampaign, HTSResult, Compound
    print("✅ All imports valid")
    sys.exit(0)
except ImportError as e:
    print(f"❌ Import error: {e}")
    sys.exit(1)

