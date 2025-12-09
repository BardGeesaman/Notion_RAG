#!/usr/bin/env python3
"""
Notion Removal Progress Tracker.

Generates a detailed progress report comparing current state to Phase 1.1 baseline.
Shows which Phase 2 tasks are complete and which remain.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Baseline from Phase 1.1 dependency report
PHASE_1_1_BASELINE = {
    "notion_client": 48,
    "notion_pages": 10,
    "dataset_notion_utils": 11,
    "signature_notion_crud": 2,
    "pathway_notion_integration": 1,
    "chemistry_notion_integration": 4,
    "dual_write": 1,
}

# Phase 2 task mappings (from NOTION_REMOVAL_AGENT_DELEGATION_PLAN.md)
PHASE_2_TASKS = {
    "amprenta_rag/query/rag/query.py": "Phase 2.1 ✅",
    "amprenta_rag/ingestion/dataset_ingestion.py": "Phase 2.2",
    "amprenta_rag/ingestion/metabolomics/embedding.py": "Phase 2.3",
    "amprenta_rag/ingestion/proteomics/embedding.py": "Phase 2.4",
    "amprenta_rag/ingestion/lipidomics/embedding.py": "Phase 2.5",
    "amprenta_rag/ingestion/transcriptomics/embedding.py": "Phase 2.6",
    "amprenta_rag/query/cross_omics/helpers.py": "Phase 2.7",
    "amprenta_rag/query/cross_omics/program_summary.py": "Phase 2.8",
    "amprenta_rag/query/cross_omics/signature_summary.py": "Phase 2.9",
    "amprenta_rag/query/cross_omics/dataset_summary.py": "Phase 2.10",
    "amprenta_rag/rag/postgres_resolver.py": "Phase 2.11",
    "amprenta_rag/analysis/enrichment.py": "Phase 2.12",
}

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

SEARCH_DIRS = [Path("amprenta_rag"), Path("scripts")]
EXCLUDE_PATTERNS = ["__pycache__", ".pyc", "docs/", "tests/", "test_"]


def should_exclude_file(file_path: Path) -> bool:
    """Check if file should be excluded."""
    path_str = str(file_path)
    return any(pattern in path_str for pattern in EXCLUDE_PATTERNS)


def find_notion_imports() -> Dict[str, Set[str]]:
    """Find all files with Notion imports."""
    results = defaultdict(set)
    
    for search_dir in SEARCH_DIRS:
        if not search_dir.exists():
            continue
            
        for py_file in search_dir.rglob("*.py"):
            if should_exclude_file(py_file):
                continue
                
            try:
                content = py_file.read_text(encoding="utf-8")
                
                for module_name, patterns in NOTION_MODULES.items():
                    for pattern in patterns:
                        if re.search(pattern, content):
                            results[module_name].add(str(py_file))
                            break
                                
            except Exception:
                pass
    
    return dict(results)


def generate_progress_report() -> None:
    """Generate comprehensive progress report."""
    current_state = find_notion_imports()
    
    # Calculate metrics
    current_counts = {mod: len(files) for mod, files in current_state.items()}
    total_current = len(set().union(*current_state.values())) if current_state else 0
    total_baseline = sum(PHASE_1_1_BASELINE.values())
    
    progress_pct = ((total_baseline - total_current) / total_baseline * 100) if total_baseline > 0 else 100
    
    print("=" * 80)
    print("NOTION REMOVAL PROGRESS REPORT")
    print("=" * 80)
    print()
    print(f"Overall Progress: {progress_pct:.1f}% ({total_baseline - total_current}/{total_baseline} files cleaned)")
    print()
    
    # Module-by-module breakdown
    print("Module-by-Module Progress:")
    print("-" * 80)
    
    for module in sorted(PHASE_1_1_BASELINE.keys()):
        baseline = PHASE_1_1_BASELINE[module]
        current = current_counts.get(module, 0)
        reduction = baseline - current
        pct = (reduction / baseline * 100) if baseline > 0 else 100
        
        status = "✅" if current == 0 else "⚠️"
        print(f"{status} {module:30s} {baseline:3d} → {current:3d} files ({pct:5.1f}% reduction)")
    
    print()
    
    # Phase 2 task status
    print("Phase 2 Task Status:")
    print("-" * 80)
    
    completed_tasks = []
    remaining_tasks = []
    
    for file_path, task in PHASE_2_TASKS.items():
        file_obj = Path(file_path)
        if file_obj.exists():
            has_imports = False
            for module_files in current_state.values():
                if file_path in module_files:
                    has_imports = True
                    break
            
            if not has_imports:
                completed_tasks.append((file_path, task))
            else:
                remaining_tasks.append((file_path, task))
        else:
            remaining_tasks.append((file_path, task))
    
    print(f"✅ Completed: {len(completed_tasks)}/{len(PHASE_2_TASKS)} tasks")
    for file_path, task in completed_tasks:
        print(f"   {task:15s} {file_path}")
    
    print()
    print(f"⏳ Remaining: {len(remaining_tasks)}/{len(PHASE_2_TASKS)} tasks")
    for file_path, task in remaining_tasks[:10]:  # Show first 10
        print(f"   {task:15s} {file_path}")
    if len(remaining_tasks) > 10:
        print(f"   ... and {len(remaining_tasks) - 10} more")
    
    print()
    
    # Files by category
    all_files = set().union(*current_state.values()) if current_state else set()
    core_files = [f for f in all_files if "amprenta_rag/" in f]
    script_files = [f for f in all_files if "scripts/" in f]
    
    print("Files Still Needing Cleanup:")
    print("-" * 80)
    print(f"Core modules: {len(core_files)}")
    print(f"Scripts: {len(script_files)}")
    print()
    
    if core_files:
        print("Core modules (priority):")
        for f in sorted(core_files)[:15]:
            task = PHASE_2_TASKS.get(f, "Phase 5+")
            print(f"  - {f} ({task})")
        if len(core_files) > 15:
            print(f"  ... and {len(core_files) - 15} more")
    
    print()
    print("=" * 80)


def main():
    """Main entry point."""
    generate_progress_report()


if __name__ == "__main__":
    main()

