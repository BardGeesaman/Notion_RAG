#!/usr/bin/env python3
"""
Test automatic signature detection from content.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from amprenta_rag.ingestion.signature_detection import (
    detect_signature_keywords,
    extract_embedded_signature_table,
    extract_signature_from_text_table,
)

# Test content with signature
test_content = """Lipid Signature Panel Results

We identified the following ceramide signature components:

Cer(d18:1/16:0) ↑ 1.0
Cer(d18:1/22:0) ↑ 0.8
SM(d18:1/24:1) ↓ 0.5

These lipid species form a core biomarker panel for disease detection.
"""

print("Testing Automatic Signature Detection")
print("=" * 80)

# Test 1: Keyword detection
print("\n1. Keyword Detection:")
has_keywords = detect_signature_keywords(test_content)
print(f"   Result: {has_keywords}")
print(f"   Expected: True")

# Test 2: Embedded table extraction
print("\n2. Embedded Table Extraction:")
embedded_table = extract_embedded_signature_table(test_content)
if embedded_table:
    print(f"   Result: Found {len(embedded_table)} rows")
    for row in embedded_table:
        print(f"     - {row}")
else:
    print(f"   Result: None found")

# Test 3: Text table extraction
print("\n3. Text Table Extraction:")
text_table = extract_signature_from_text_table(test_content)
if text_table and text_table.get("components"):
    print(f"   Result: Found {len(text_table['components'])} components")
    for comp in text_table['components']:
        print(f"     - {comp}")
else:
    print(f"   Result: None found")

print("\n" + "=" * 80)

