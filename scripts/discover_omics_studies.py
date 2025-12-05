#!/usr/bin/env python3
"""
Unified discovery script for omics studies across all repositories.

Searches for studies in GEO, PRIDE, MetaboLights, and Metabolomics Workbench
with a unified interface.
"""

import argparse
import json
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repositories.discovery import (
    discover_studies,
    list_available_repositories,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Discover omics studies across public repositories (GEO, PRIDE, MetaboLights, MW)"
    )
    parser.add_argument(
        "--keywords",
        nargs="+",
        required=False,
        help="Search keywords (e.g., 'ALS' 'Alzheimer')",
    )
    parser.add_argument(
        "--omics-type",
        choices=["transcriptomics", "proteomics", "metabolomics", "lipidomics"],
        help="Filter by omics type",
    )
    parser.add_argument(
        "--repository",
        help="Search specific repository (GEO, PRIDE, MetaboLights, MW, MW_LIPIDOMICS, MW_METABOLOMICS)",
    )
    parser.add_argument(
        "--disease",
        help="Filter by disease (e.g., 'ALS', 'Alzheimer')",
    )
    parser.add_argument(
        "--organism",
        help="Filter by organism (e.g., 'Homo sapiens', 'Mus musculus')",
    )
    parser.add_argument(
        "--sample-type",
        help="Filter by sample type (e.g., 'CSF', 'plasma', 'tissue')",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=50,
        help="Maximum results per repository (default: 50)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file for study results",
    )
    parser.add_argument(
        "--list-repositories",
        action="store_true",
        help="List all available repositories and exit",
    )

    args = parser.parse_args()

    # List repositories if requested
    if args.list_repositories:
        repos = list_available_repositories()
        print("\nAvailable repositories:")
        for repo in repos:
            print(f"  - {repo}")
        print()
        return

    # Require keywords if not listing repositories
    if not args.keywords:
        parser.error("--keywords is required (or use --list-repositories)")

    # Build filters
    filters: dict[str, any] = {}
    if args.disease:
        filters["disease"] = args.disease
    if args.organism:
        filters["organism"] = args.organism
    if args.sample_type:
        filters["sample_type"] = args.sample_type

    # Discover studies
    try:
        results = discover_studies(
            keywords=args.keywords,
            omics_type=args.omics_type,
            repository=args.repository,
            filters=filters if filters else None,
            max_results=args.max_results,
        )

        # Print results
        print("\n" + "=" * 80)
        print("OMICS STUDY DISCOVERY RESULTS")
        print("=" * 80)
        print(f"\nKeywords: {', '.join(args.keywords)}")
        if args.omics_type:
            print(f"Omics Type: {args.omics_type}")
        if args.repository:
            print(f"Repository: {args.repository}")
        if filters:
            print(f"Filters: {filters}")

        total_studies = sum(len(ids) for ids in results.values())
        print(f"\nTotal studies found: {total_studies}\n")

        for repo_name, study_ids in results.items():
            if study_ids:
                print(f"{repo_name}: {len(study_ids)} studies")
                for i, study_id in enumerate(study_ids[:10], 1):
                    print(f"  {i:2d}. {study_id}")
                if len(study_ids) > 10:
                    print(f"  ... and {len(study_ids) - 10} more")
                print()

        # Save to file if requested
        if args.output:
            output_data = {
                "keywords": args.keywords,
                "omics_type": args.omics_type,
                "repository": args.repository,
                "filters": filters,
                "results": results,
                "total_studies": total_studies,
            }
            args.output.write_text(json.dumps(output_data, indent=2))
            print(f"✅ Results saved to {args.output}")

        print("=" * 80)

    except Exception as e:
        logger.error("[DISCOVER] Error during discovery: %r", e)
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

