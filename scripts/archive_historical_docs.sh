#!/bin/bash
# Archive historical documentation files
# Created: December 9, 2025
# Purpose: Move completed/historical docs from docs/ to md_archive/

set -e

DOCS_DIR="/Users/bard/Documents/RAG/docs"
ARCHIVE_DIR="/Users/bard/Documents/RAG/md_archive"

echo "📦 Archiving historical documentation..."
echo ""

# Create archive directory if it doesn't exist
mkdir -p "$ARCHIVE_DIR"

# Counter
archived=0

# Migration/Status Files
echo "📁 Archiving migration files..."
files=(
    "COMPLETE_FEATURE_MIGRATION_COMPARISON.md"
    "COMPLETE_MIGRATION_SUMMARY.md"
    "COMPLETE_NOTION_MIGRATION_COMPLETE.md"
    "COMPLETE_NOTION_MIGRATION_PLAN.md"
    "COMPLETE_REMAINING_ROADMAP_ITEMS.md"
    "NOTION_MIGRATION_100_PERCENT_COMPLETE.md"
    "NOTION_MIGRATION_FINAL_STATUS.md"
    "NOTION_MIGRATION_PLAN.md"
    "NOTION_MIGRATION_STATUS.md"
    "NOTION_REMOVAL_AGENT_DELEGATION_PLAN.md"
    "NOTION_REMOVAL_IMPLEMENTATION_INSTRUCTIONS.md"
    "FINAL_NOTION_MIGRATION_AUDIT.md"
    "FINAL_NOTION_REMOVAL_STATUS.md"
    "FINAL_NOTION_REMOVAL_VERIFICATION_REPORT.md"
    "POSTGRES_MIGRATION_GUIDE.md"
    "POSTGRES_MIGRATION_STRATEGY.md"
    "POSTGRES_SOT_TRANSITION.md"
    "MIGRATION_STRATEGY.md"
    "MIGRATION_TIMEOUT_DIAGNOSIS.md"
    "MIGRATION_TIMEOUT_FIX.md"
    "MIGRATION_TIMEOUT_RESOLUTION.md"
    "MIGRATION_TIMEOUT_RESOLVED.md"
)

for file in "${files[@]}"; do
    if [ -f "$DOCS_DIR/$file" ]; then
        mv "$DOCS_DIR/$file" "$ARCHIVE_DIR/"
        echo "  ✓ $file"
        ((archived++))
    fi
done

# Completion/Status Reports
echo ""
echo "📁 Archiving completion reports..."
files=(
    "GAP_FILLING_COMPLETE.md"
    "GAP_FILLING_PROGRESS.md"
    "FEATURE_LINKING_COMPLETE.md"
    "GENOMICS_PIPELINE_COMPLETE.md"
    "GENOMICS_PIPELINE_IMPLEMENTATION.md"
    "MASTER_PROTOCOL_IMPLEMENTATION_COMPLETE.md"
    "METABOLOMICS_WORKBENCH_EXTRACTION_COMPLETE.md"
    "POSTGRES_SIGNATURE_CREATION_COMPLETE.md"
    "PRIDE_OPTIMIZATION_COMPLETE.md"
    "PRIORITY_3_COMPLETE.md"
    "PRIORITY_4_COMPLETE.md"
    "PRIORITY_4_PROGRESS.md"
    "SALMON_IMPLEMENTATION_COMPLETE.md"
    "SIGNATURE_INTEGRATION_COMPLETE.md"
    "TIER3_COMPLETE.md"
    "OPTIONAL_ENHANCEMENTS_COMPLETE.md"
    "OPTIONAL_ENHANCEMENTS_IMPLEMENTATION_PLAN.md"
    "OPTIONAL_ENHANCEMENTS_STATUS.md"
    "FINAL_IMPLEMENTATION_SUMMARY.md"
    "IMPLEMENTATION_VERIFICATION_SUMMARY.md"
    "TIER3_INGESTION_INTEGRATION.md"
    "TIER3_ARCHITECTURE_EVOLUTION.md"
)

for file in "${files[@]}"; do
    if [ -f "$DOCS_DIR/$file" ]; then
        mv "$DOCS_DIR/$file" "$ARCHIVE_DIR/"
        echo "  ✓ $file"
        ((archived++))
    fi
done

# Test Results
echo ""
echo "📁 Archiving test results..."
files=(
    "ENA_INGESTION_TEST_RESULTS.md"
    "GEO_API_KEY_TEST_RESULTS.md"
    "GEO_REPOSITORY_TEST_RESULTS.md"
    "TEST_IMPORT_RESULTS.md"
    "TIER3_TESTING_CONTINUED.md"
    "TIER3_TESTING_FINAL_RESULTS.md"
    "TIER3_TESTING_GUIDE.md"
    "TIER3_TESTING_RESULTS.md"
    "TIER3_TESTING_STATUS.md"
    "TIER3_TESTING_SUMMARY.md"
    "REPOSITORY_INGESTION_TEST_RESULTS.md"
    "REPOSITORY_INGESTION_TEST_PLAN.md"
    "POST_INCIDENT_VERIFICATION_REPORT.md"
)

for file in "${files[@]}"; do
    if [ -f "$DOCS_DIR/$file" ]; then
        mv "$DOCS_DIR/$file" "$ARCHIVE_DIR/"
        echo "  ✓ $file"
        ((archived++))
    fi
done

# Historical Implementation Docs
echo ""
echo "📁 Archiving historical implementation docs..."
files=(
    "ALS_GEO_STUDIES_FOUND.md"
    "ALS_STUDY_OMICS_DATA_INSPECTION.md"
    "GEO_CONFIGURATION_COMPLETE.md"
    "GEO_FEATURE_EXTRACTION_SUCCESS.md"
    "GEO_GEOparse_MIGRATION_COMPLETE.md"
    "GEO_GEOparse_RECOMMENDATION.md"
    "GEO_IMPORT_SUCCESS.md"
    "GEO_PRIDE_FIXES.md"
    "GEO_REPOSITORY_REWRITE_COMPLETE.md"
    "METABOLIGHTS_EXTRACTION_STATUS.md"
    "METABOLIGHTS_MAF_EXTRACTION.md"
    "METABOLIGHTS_MAF_SEARCH_STATUS.md"
    "METABOLIGHTS_FILES_ENDPOINT_LIMITATION.md"
    "FEATURE_EXTRACTION_OPTIMIZATION_LEARNINGS.md"
    "FEATURE_EXTRACTION_STATUS.md"
    "PATHWAY_ANALYSIS_ID_MAPPING_STATUS.md"
    "INSTRUCTIONS_AUTOMATOR_PHASE1_1.md"
)

for file in "${files[@]}"; do
    if [ -f "$DOCS_DIR/$file" ]; then
        mv "$DOCS_DIR/$file" "$ARCHIVE_DIR/"
        echo "  ✓ $file"
        ((archived++))
    fi
done

# Redundant Status/Summary Files
echo ""
echo "📁 Archiving redundant status files..."
files=(
    "ACTUAL_STATUS.md"
    "CURRENT_STATUS.md"
    "COMPREHENSIVE_ROADMAP_STATUS.md"
    "IMPROVEMENTS_SUMMARY.md"
    "LESSONS_LEARNED_DEC_2025.md"
    "API_IMPLEMENTATION_STATUS.md"
    "API_KEY_CONFIGURATION_STATUS.md"
    "AUTHENTICATION_SUMMARY.md"
    "FEATURE_EXTRACTION_STATUS.md"
    "REPOSITORY_STATUS.md"
    "TESTING_STATUS.md"
    "REPOSITORY_INGESTION_ENHANCEMENTS.md"
    "REPOSITORY_OPTIMIZATION_WORKFLOW.md"
    "TIER3_DEPENDENCIES_INSTALLED.md"
    "UI_IMPROVEMENT_ROADMAP.md"
    "AI_AGENT_TOOLS_UI_PLAN.md"
)

for file in "${files[@]}"; do
    if [ -f "$DOCS_DIR/$file" ]; then
        mv "$DOCS_DIR/$file" "$ARCHIVE_DIR/"
        echo "  ✓ $file"
        ((archived++))
    fi
done

# Optional: Archive some Notion-specific docs (now that Notion removal is complete)
echo ""
echo "📁 Archiving Notion-specific implementation docs..."
files=(
    "NOTION_DATABASE_SETUP.md"
    "NOTION_DISABLED_BY_DEFAULT.md"
)

for file in "${files[@]}"; do
    if [ -f "$DOCS_DIR/$file" ]; then
        mv "$DOCS_DIR/$file" "$ARCHIVE_DIR/"
        echo "  ✓ $file"
        ((archived++))
    fi
done

echo ""
echo "✅ Archive complete!"
echo "📊 Archived $archived files to md_archive/"
echo ""
echo "📁 Remaining active docs in docs/:"
ls -1 "$DOCS_DIR" | wc -l

