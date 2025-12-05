#!/bin/bash
# Quick import script for omics data from repositories
# Usage: ./scripts/quick_import_omics.sh

set -e

echo "=========================================="
echo "Omics Repository Import - Quick Start"
echo "=========================================="
echo ""

# Check if Python script exists
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMPORT_SCRIPT="$SCRIPT_DIR/import_all_omics_repositories.py"

if [ ! -f "$IMPORT_SCRIPT" ]; then
    echo "Error: $IMPORT_SCRIPT not found"
    exit 1
fi

# Default options
DISEASE=""
OMICS_TYPE=""
KEYWORDS=""
MAX_RESULTS=20
INGEST=true
CREATE_NOTION=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --disease)
            DISEASE="$2"
            shift 2
            ;;
        --omics-type)
            OMICS_TYPE="$2"
            shift 2
            ;;
        --keywords)
            KEYWORDS="$2"
            shift 2
            ;;
        --max-results)
            MAX_RESULTS="$2"
            shift 2
            ;;
        --no-ingest)
            INGEST=false
            shift
            ;;
        --create-notion)
            CREATE_NOTION=true
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --disease DISEASE          Filter by disease (e.g., 'ALS', 'Alzheimer')"
            echo "  --omics-type TYPE          Filter by omics type (transcriptomics, proteomics, metabolomics, lipidomics)"
            echo "  --keywords KEYWORDS        Search keywords (comma-separated)"
            echo "  --max-results N            Maximum results per repository (default: 20)"
            echo "  --no-ingest                Don't ingest to Pinecone"
            echo "  --create-notion            Also create Notion pages"
            echo "  --help                     Show this help"
            echo ""
            echo "Examples:"
            echo "  $0 --disease ALS --omics-type metabolomics"
            echo "  $0 --keywords 'Alzheimer,amyloid' --omics-type transcriptomics"
            echo "  $0 --disease 'Parkinson' --max-results 10"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Build command
CMD="python $IMPORT_SCRIPT"

if [ -n "$DISEASE" ]; then
    CMD="$CMD --disease \"$DISEASE\""
fi

if [ -n "$OMICS_TYPE" ]; then
    CMD="$CMD --omics-type $OMICS_TYPE"
fi

if [ -n "$KEYWORDS" ]; then
    # Split comma-separated keywords
    IFS=',' read -ra KW_ARRAY <<< "$KEYWORDS"
    CMD="$CMD --keywords"
    for kw in "${KW_ARRAY[@]}"; do
        CMD="$CMD \"${kw// /}\""
    done
fi

CMD="$CMD --max-results $MAX_RESULTS"

CMD="$CMD --import"

if [ "$INGEST" = true ]; then
    CMD="$CMD --ingest"
fi

if [ "$CREATE_NOTION" = true ]; then
    CMD="$CMD --create-notion"
fi

echo "Running: $CMD"
echo ""
eval $CMD

