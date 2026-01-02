#!/usr/bin/env python3
"""
Custom pre-commit checks for Amprenta RAG.

Catches common issues before they're committed:
1. Module-level celery imports in router files
2. Deprecated datetime.utcnow() usage
3. Schemas inheriting BaseModel instead of BaseSchema

Usage:
    python scripts/pre_commit_checks.py [files...]
    
Exit codes:
    0 = All checks passed
    1 = Errors found (blocks commit)
"""

import re
import sys
from pathlib import Path
from typing import List, Tuple


def check_celery_imports(filepath: Path, content: str) -> List[str]:
    """
    Block module-level celery imports in router files.
    
    Celery imports should be inside functions to avoid import chain issues
    that break test environments.
    """
    errors = []
    
    # Only check router files
    if '/routers/' not in str(filepath) and '\\routers\\' not in str(filepath):
        return errors
    
    # Pattern for celery task imports
    celery_import_pattern = re.compile(
        r'^(?:from\s+amprenta_rag\.jobs\.tasks\S*\s+import|'
        r'import\s+amprenta_rag\.jobs\.tasks)',
        re.MULTILINE
    )
    
    # Check if import exists
    if celery_import_pattern.search(content):
        # Verify it's at module level (not inside a function)
        # Simple heuristic: if the import line has no leading whitespace, it's module-level
        for i, line in enumerate(content.split('\n'), 1):
            if celery_import_pattern.match(line):
                errors.append(
                    f"{filepath}:{i}: Module-level celery import detected. "
                    f"Move import inside the function that uses it to avoid import chain issues."
                )
    
    return errors


def check_utcnow(filepath: Path, content: str) -> List[str]:
    """
    Block deprecated datetime.utcnow() usage.
    
    datetime.utcnow() is deprecated in Python 3.12+.
    Use datetime.now(timezone.utc) instead.
    """
    errors = []
    
    pattern = re.compile(r'datetime\.utcnow\s*\(\)')
    
    for i, line in enumerate(content.split('\n'), 1):
        if pattern.search(line):
            errors.append(
                f"{filepath}:{i}: datetime.utcnow() is deprecated. "
                f"Use datetime.now(timezone.utc) instead."
            )
    
    return errors


def check_schema_inheritance(filepath: Path, content: str) -> List[Tuple[str, bool]]:
    """
    Warn on schemas inheriting BaseModel instead of BaseSchema.
    
    BaseSchema includes from_attributes=True which is needed for
    serializing dataclasses and ORM objects.
    
    Returns list of (message, is_error) tuples. Schema inheritance
    issues are warnings, not errors.
    """
    warnings = []
    
    # Only check schema files
    if 'schemas.py' not in str(filepath):
        return warnings
    
    # Skip BaseSchema and StrictBaseSchema definitions themselves
    pattern = re.compile(r'class\s+(\w+Schema)\(BaseModel\)')
    
    for i, line in enumerate(content.split('\n'), 1):
        match = pattern.search(line)
        if match:
            class_name = match.group(1)
            # Skip the base classes themselves
            if class_name in ('BaseSchema', 'StrictBaseSchema'):
                continue
            warnings.append((
                f"{filepath}:{i}: {class_name} inherits from BaseModel. "
                f"Consider using BaseSchema for automatic from_attributes=True.",
                False  # Warning, not error
            ))
    
    return warnings


def main(files: List[str]) -> int:
    """
    Run all checks on provided files.
    
    Returns:
        0 if no errors, 1 if errors found
    """
    errors = []
    warnings = []
    
    for filepath_str in files:
        filepath = Path(filepath_str)
        
        # Skip non-Python files
        if filepath.suffix != '.py':
            continue
        
        # Skip if file doesn't exist (might be deleted)
        if not filepath.exists():
            continue
        
        try:
            content = filepath.read_text(encoding='utf-8')
        except Exception as e:
            print(f"Warning: Could not read {filepath}: {e}", file=sys.stderr)
            continue
        
        # Run checks
        errors.extend(check_celery_imports(filepath, content))
        errors.extend(check_utcnow(filepath, content))
        
        for msg, is_error in check_schema_inheritance(filepath, content):
            if is_error:
                errors.append(msg)
            else:
                warnings.append(msg)
    
    # Print warnings (don't block commit)
    for warning in warnings:
        print(f"⚠️  WARNING: {warning}", file=sys.stderr)
    
    # Print errors (block commit)
    for error in errors:
        print(f"❌ ERROR: {error}", file=sys.stderr)
    
    if errors:
        print(f"\n{len(errors)} error(s) found. Commit blocked.", file=sys.stderr)
        return 1
    
    if warnings:
        print(f"\n{len(warnings)} warning(s) found. Commit allowed.", file=sys.stderr)
    
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

