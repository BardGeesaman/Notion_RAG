#!/usr/bin/env python3
"""
Check Remaining Notion Imports.

Scans the codebase for remaining imports from deleted Notion modules.
Used to track Phase 2 progress and identify files that still need cleanup.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Notion modules that were deleted in Phase 1.2
NOTION_MODULES = {
    "notion_client": [
        r"from amprenta_rag\.clients\.notion_client",
        r"import.*notion_client",
    ],
    "notion_pages": [
        r"from amprenta_rag\.ingestion\.notion_pages",
        r"import.*notion_pages",
    ],
    "dataset_notion_utils": [
        r"from amprenta_rag\.ingestion\.dataset_notion_utils",
        r"import.*dataset_notion_utils",
    ],
    "signature_notion_crud": [
        r"from amprenta_rag\.ingestion\.signature_notion_crud",
        r"import.*signature_notion_crud",
    ],
    "pathway_notion_integration": [
        r"from amprenta_rag\.analysis\.pathway_notion_integration",
        r"import.*pathway_notion_integration",
    ],
    "chemistry_notion_integration": [
        r"from amprenta_rag\.chemistry\.notion_integration",
        r"from amprenta_rag\.chemistry import.*notion",
    ],
    "dual_write": [
        r"from amprenta_rag\.migration\.dual_write",
        r"import.*dual_write",
    ],
}

# Directories to search
SEARCH_DIRS = [
    Path("amprenta_rag"),
    Path("scripts"),
]

# Directories/files to exclude
EXCLUDE_PATTERNS = [
    "__pycache__",
    ".pyc",
    "docs/",
    "tests/",
    "test_",
]


def should_exclude_file(file_path: Path) -> bool:
    """Check if file should be excluded from search."""
    path_str = str(file_path)
    return any(pattern in path_str for pattern in EXCLUDE_PATTERNS)


def find_notion_imports() -> Dict[str, List[Dict[str, str]]]:
    """Find all remaining Notion imports in the codebase."""
    results = defaultdict(list)
    
    for search_dir in SEARCH_DIRS:
        if not search_dir.exists():
            continue
            
        for py_file in search_dir.rglob("*.py"):
            if should_exclude_file(py_file):
                continue
                
            try:
                content = py_file.read_text(encoding="utf-8")
                lines = content.split("\n")
                
                for module_name, patterns in NOTION_MODULES.items():
                    for pattern in patterns:
                        for line_num, line in enumerate(lines, 1):
                            if re.search(pattern, line):
                                results[module_name].append({
                                    "file": str(py_file),
                                    "line": line_num,
                                    "content": line.strip(),
                                })
                                break  # Only report once per file per module
                                
            except Exception as e:
                print(f"Warning: Could not read {py_file}: {e}", file=sys.stderr)
    
    return dict(results)


def generate_report(results: Dict[str, List[Dict[str, str]]], verbose: bool = False) -> None:
    """Generate a formatted report of remaining imports."""
    total_files = set()
    total_imports = 0
    
    print("=" * 80)
    print("NOTION IMPORT REMOVAL PROGRESS REPORT")
    print("=" * 80)
    print()
    
    if not results:
        print("✅ SUCCESS: No Notion imports found!")
        print()
        return
    
    print(f"Found imports from {len(results)} deleted Notion modules:\n")
    
    for module_name, imports in sorted(results.items()):
        files = {imp["file"] for imp in imports}
        total_files.update(files)
        total_imports += len(imports)
        
        print(f"📦 {module_name}")
        print(f"   Files: {len(files)}")
        print(f"   Import statements: {len(imports)}")
        
        if verbose:
            for imp in imports:
                print(f"      - {imp['file']}:{imp['line']}")
                print(f"        {imp['content']}")
        else:
            # Show first 5 files
            for imp in imports[:5]:
                print(f"      - {imp['file']}")
            if len(imports) > 5:
                print(f"      ... and {len(imports) - 5} more")
        print()
    
    print("=" * 80)
    print(f"SUMMARY")
    print("=" * 80)
    print(f"Total files with Notion imports: {len(total_files)}")
    print(f"Total import statements: {total_imports}")
    print()
    
    if total_files:
        print("⚠️  Action required: Remove imports from these files in Phase 2")
        print()
        print("Files by category:")
        
        # Categorize files
        core_files = [f for f in total_files if "amprenta_rag/" in f]
        script_files = [f for f in total_files if "scripts/" in f]
        
        print(f"  Core modules: {len(core_files)}")
        print(f"  Scripts: {len(script_files)}")
        print()
        
        if verbose:
            print("All files:")
            for f in sorted(total_files):
                print(f"  - {f}")


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Check for remaining Notion imports in the codebase"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Show detailed file listings"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON"
    )
    
    args = parser.parse_args()
    
    results = find_notion_imports()
    
    if args.json:
        import json
        print(json.dumps(results, indent=2))
    else:
        generate_report(results, verbose=args.verbose)
        
        # Exit with error code if imports found
        if results:
            sys.exit(1)
        else:
            sys.exit(0)


if __name__ == "__main__":
    main()

