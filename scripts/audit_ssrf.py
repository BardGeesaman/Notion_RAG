"""Audit HTTP calls for SSRF vulnerabilities."""

import re
from pathlib import Path

HTTP_PATTERN = re.compile(r'requests\.(get|post|put|delete)\s*\(([^)]+)\)')

def audit_file(filepath: Path) -> list[dict]:
    """Check if HTTP calls use hardcoded or variable URLs."""
    results = []
    
    try:
        content = filepath.read_text()
    except Exception:
        return results
    
    for match in HTTP_PATTERN.finditer(content):
        method = match.group(1)
        args = match.group(2)
        
        # Check if URL is hardcoded string or variable
        is_hardcoded = args.strip().startswith('"') or args.strip().startswith("'")
        is_fstring = 'f"' in args or "f'" in args
        
        results.append({
            'file': str(filepath),
            'method': method,
            'url_type': 'hardcoded' if is_hardcoded else ('f-string' if is_fstring else 'variable'),
            'snippet': args[:80],
            'needs_review': not is_hardcoded
        })
    
    return results

def main():
    """Run SSRF audit on all Python files."""
    all_results = []
    
    # Check all Python files in amprenta_rag
    for py_file in Path("amprenta_rag").rglob("*.py"):
        results = audit_file(py_file)
        all_results.extend(results)
    
    # Summary
    total = len(all_results)
    needs_review = sum(1 for r in all_results if r['needs_review'])
    
    print(f"SSRF Audit Results:")
    print(f"- Total HTTP calls: {total}")
    print(f"- Variable URLs (needs review): {needs_review}")
    print(f"- Hardcoded URLs (safe): {total - needs_review}")
    
    if needs_review > 0:
        print(f"\nFiles needing review:")
        for result in all_results:
            if result['needs_review']:
                print(f"  {result['file']}:{result['method']} - {result['url_type']}")

if __name__ == "__main__":
    main()
