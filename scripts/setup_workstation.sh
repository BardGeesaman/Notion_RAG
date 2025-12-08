#!/bin/bash

# Automated workstation setup for Amprenta RAG Platform

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "🔧 Setting up Amprenta workstation..."

# Detect username
USERNAME=$(whoami)
echo "Detected username: $USERNAME"

# Create .env from template if not exists
if [ ! -f .env ]; then
    cp .env.template .env
    sed -i '' "s/YOUR_USERNAME_HERE/$USERNAME/g" .env
    echo "✅ Created .env with username: $USERNAME"
else
    echo "⚠️  .env already exists, skipping"
fi

# Install Python dependencies
echo "Installing Python packages..."
pip install -r requirements.txt

# Install Playwright browsers
echo "Installing Playwright browsers..."
playwright install chromium

# Clear Python caches
echo "Clearing Python caches..."
find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null
find . -name "*.pyc" -delete

# Create logs directory
mkdir -p logs

echo ""
echo "✅ Setup complete!"
echo ""
echo "Next steps:"
echo "  1. Edit .env if needed (database credentials, API keys)"
echo "  2. Run tests: pytest amprenta_rag/tests/dashboard/test_pages_import.py"
echo "  3. Start dashboard: streamlit run scripts/run_dashboard.py"
echo "  4. Start API: uvicorn amprenta_rag.api.main:app --reload"


