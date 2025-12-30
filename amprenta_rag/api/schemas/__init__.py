"""API schemas package."""

# Import all schemas from the main schemas.py file
# Import everything from the parent schemas.py file
import importlib.util
import sys
from pathlib import Path

# Get path to schemas.py in parent directory
schemas_path = Path(__file__).parent / "schemas.py"

if schemas_path.exists():
    spec = importlib.util.spec_from_file_location("schemas", schemas_path)
    schemas_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(schemas_module)
    
    # Import all public attributes from schemas module
    for attr_name in dir(schemas_module):
        if not attr_name.startswith('_'):
            globals()[attr_name] = getattr(schemas_module, attr_name)
