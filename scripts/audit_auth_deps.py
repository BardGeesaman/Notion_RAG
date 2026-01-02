"""Audit API endpoints for missing authentication dependencies."""

import os
import re
import csv
from pathlib import Path

ROUTERS_DIR = Path("amprenta_rag/api/routers")
OUTPUT_FILE = "docs/endpoint_auth_matrix.csv"

# Patterns
ENDPOINT_PATTERN = re.compile(r'@router\.(get|post|put|delete|patch)\s*\(\s*["\']([^"\']+)["\']')
AUTH_PATTERN = re.compile(r'get_current_user|get_current_user_with_company')

def audit_file(filepath: Path) -> list[dict]:
    """Audit single router file for auth dependencies."""
    results = []
    content = filepath.read_text()
    lines = content.split('\n')
    
    for i, line in enumerate(lines):
        match = ENDPOINT_PATTERN.search(line)
        if match:
            method = match.group(1).upper()
            path = match.group(2)
            
            # Check next 20 lines for auth dependency
            func_block = '\n'.join(lines[i:i+20])
            has_auth = bool(AUTH_PATTERN.search(func_block))
            
            results.append({
                'file': filepath.name,
                'method': method,
                'path': path,
                'line': i + 1,
                'has_auth': has_auth,
                'needs_review': not has_auth
            })
    
    return results

def main():
    all_results = []
    
    for router_file in ROUTERS_DIR.glob("*.py"):
        if router_file.name.startswith("__"):
            continue
        results = audit_file(router_file)
        all_results.extend(results)
    
    # Write CSV
    with open(OUTPUT_FILE, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['file', 'method', 'path', 'line', 'has_auth', 'needs_review'])
        writer.writeheader()
        writer.writerows(all_results)
    
    # Summary
    total = len(all_results)
    missing = sum(1 for r in all_results if r['needs_review'])
    print(f"Audited {total} endpoints across {len(list(ROUTERS_DIR.glob('*.py')))} files")
    print(f"Missing auth: {missing} ({missing/total*100:.1f}%)")
    print(f"Output: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
