#!/bin/bash
# Verify environment before running tests
set -e

echo "Checking environment..."

# Check if conda environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "myenv" ]]; then
    echo "ERROR: Conda environment 'myenv' is not activated"
    echo "Run: conda activate myenv"
    exit 1
fi

# Check if critical packages are installed
python -c "import fastapi, sqlalchemy, pytest, playwright" 2>/dev/null || {
    echo "ERROR: Missing required packages"
    echo "Run: pip install -r requirements.txt"
    exit 1
}

echo "✓ Environment OK"

