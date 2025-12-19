#!/usr/bin/env python3
"""
Batch import script for repository studies.

Imports multiple studies from public repositories (GEO, PRIDE, MetaboLights, MW)
in a single run. Supports:
- JSON input file (from discovery script)
- Command-line study IDs
- Progress tracking
- Error handling and reporting
"""

import argparse
import json
import sys
from pathlib import Path
from typing import List, Dict, Optional

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.harvest_repository_study import harvest_study
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def load_studies_from_json(json_file: Path) -> List[Dict[str, str]]:
    """
    Load study list from JSON file (from discovery script output).

    Args:
        json_file: Path to JSON file

    Returns:
        List of dicts with 'study_id' and 'repository' keys
    """
    try:
        data = json.loads(json_file.read_text())

        studies = []

        # Handle discovery script output format
        if "results" in data:
            for repo_name, study_ids in data["results"].items():
                for study_id in study_ids:
                    studies.append({
                        "study_id": study_id,
                        "repository": repo_name,
                    })

        # Handle simple list format
        elif isinstance(data, list):
            studies = data

        logger.info(
            "[BATCH][IMPORT] Loaded %d studies from %s",
            len(studies),
            json_file,
        )

        return studies

    except Exception as e:
        logger.error(
            "[BATCH][IMPORT] Error loading studies from %s: %r",
            json_file,
            e,
        )
        raise


def parse_study_list(study_list: List[str]) -> List[Dict[str, str]]:
    """
    Parse study list from command line.

    Format: REPOSITORY:STUDY_ID or just STUDY_ID (if repository specified)

    Args:
        study_list: List of study identifiers

    Returns:
        List of dicts with 'study_id' and 'repository' keys
    """
    studies = []

    for study_spec in study_list:
        if ":" in study_spec:
            # Format: REPOSITORY:STUDY_ID
            parts = study_spec.split(":", 1)
            repository = parts[0].strip()
            study_id = parts[1].strip()
            studies.append({
                "study_id": study_id,
                "repository": repository,
            })
        else:
            # Just study ID - requires --repository argument
            studies.append({
                "study_id": study_spec.strip(),
                "repository": None,  # Will be filled from --repository
            })

    return studies


def batch_import_studies(
    studies: List[Dict[str, str]],
    default_repository: Optional[str] = None,
    create_notion: bool = False,  # DEPRECATED: Notion integration removed - these options are no-ops
    create_postgres: bool = True,
    ingest: bool = False,
    stop_on_error: bool = False,
) -> Dict[str, any]:
    """
    Batch import multiple studies.

    Args:
        studies: List of dicts with 'study_id' and 'repository' keys
        default_repository: Default repository if not specified in study dict
        create_notion: Create Notion pages
        create_postgres: Create Postgres datasets
        ingest: Trigger ingestion
        stop_on_error: Stop on first error (default: continue)

    Returns:
        Dict with results: {'succeeded': [...], 'failed': [...], 'skipped': [...]}
    """
    results = {
        "succeeded": [],
        "failed": [],
        "skipped": [],
    }

    total = len(studies)

    logger.info(
        "[BATCH][IMPORT] Starting batch import of %d studies",
        total,
    )

    for i, study_info in enumerate(studies, 1):
        study_id = study_info["study_id"]
        repository = study_info.get("repository") or default_repository

        if not repository:
            logger.error(
                "[BATCH][IMPORT] No repository specified for study %s",
                study_id,
            )
            results["failed"].append({
                "study_id": study_id,
                "repository": None,
                "error": "No repository specified",
            })
            if stop_on_error:
                break
            continue

        logger.info(
            "[BATCH][IMPORT] [%d/%d] Importing %s study %s",
            i,
            total,
            repository,
            study_id,
        )

        try:
            dataset_id, page_id = harvest_study(
                study_id=study_id,
                repository=repository,
                create_notion=create_notion,
                create_postgres=create_postgres,
                ingest=ingest,
                dry_run=False,
            )

            if dataset_id or page_id:
                results["succeeded"].append({
                    "study_id": study_id,
                    "repository": repository,
                    "dataset_id": dataset_id,
                    "page_id": page_id,
                })
                logger.info(
                    "[BATCH][IMPORT] [%d/%d] ✅ Success: %s %s",
                    i,
                    total,
                    repository,
                    study_id,
                )
            else:
                results["skipped"].append({
                    "study_id": study_id,
                    "repository": repository,
                    "reason": "No dataset or page created",
                })
                logger.warning(
                    "[BATCH][IMPORT] [%d/%d] ⚠️  Skipped: %s %s",
                    i,
                    total,
                    repository,
                    study_id,
                )

        except Exception as e:
            error_msg = str(e)
            logger.error(
                "[BATCH][IMPORT] [%d/%d] ❌ Failed: %s %s - %s",
                i,
                total,
                repository,
                study_id,
                error_msg,
            )
            results["failed"].append({
                "study_id": study_id,
                "repository": repository,
                "error": error_msg,
            })

            if stop_on_error:
                logger.error(
                    "[BATCH][IMPORT] Stopping due to error (--stop-on-error)",
                )
                break

    # Print summary
    print("\n" + "=" * 80)
    print("BATCH IMPORT SUMMARY")
    print("=" * 80)
    print(f"Total studies: {total}")
    print(f"✅ Succeeded: {len(results['succeeded'])}")
    print(f"⚠️  Skipped: {len(results['skipped'])}")
    print(f"❌ Failed: {len(results['failed'])}")
    print("=" * 80)

    if results["failed"]:
        print("\nFailed studies:")
        for failed in results["failed"]:
            print(f"  - {failed['repository']} {failed['study_id']}: {failed['error']}")

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch import studies from public repositories",
    )

    # Input methods (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--json-file",
        type=Path,
        help="JSON file with studies (from discovery script output)",
    )
    input_group.add_argument(
        "--studies",
        nargs="+",
        help="Study IDs (format: REPOSITORY:STUDY_ID or just STUDY_ID with --repository)",
    )

    parser.add_argument(
        "--repository",
        help="Default repository name (required if using --studies without REPOSITORY: prefix)",
    )

    # Import options
    parser.add_argument(
        "--create-notion",
        action="store_true",
        help="Create Notion pages",
    )
    parser.add_argument(
        "--no-postgres",
        action="store_true",
        help="Skip Postgres dataset creation",
    )
    parser.add_argument(
        "--ingest",
        action="store_true",
        help="Trigger ingestion after creating datasets",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop on first error (default: continue)",
    )

    # Output
    parser.add_argument(
        "--output",
        type=Path,
        help="Save results to JSON file",
    )

    args = parser.parse_args()

    # Load studies
    if args.json_file:
        studies = load_studies_from_json(args.json_file)
    else:
        studies = parse_study_list(args.studies)
        if not args.repository:
            # Check if all studies have repository specified
            if any(not s.get("repository") for s in studies):
                parser.error(
                    "--repository is required when using --studies without REPOSITORY: prefix",
                )

    if not studies:
        print("❌ No studies to import")
        sys.exit(1)

    # Default repository
    default_repo = args.repository
    for study in studies:
        if not study.get("repository") and default_repo:
            study["repository"] = default_repo

    # Run batch import
    try:
        results = batch_import_studies(
            studies=studies,
            default_repository=default_repo,
            create_notion=args.create_notion,
            create_postgres=not args.no_postgres,
            ingest=args.ingest,
            stop_on_error=args.stop_on_error,
        )

        # Save results if requested
        if args.output:
            output_data = {
                "summary": {
                    "total": len(studies),
                    "succeeded": len(results["succeeded"]),
                    "skipped": len(results["skipped"]),
                    "failed": len(results["failed"]),
                },
                "results": results,
            }
            args.output.write_text(json.dumps(output_data, indent=2))
            print(f"\n✅ Results saved to {args.output}")

        # Exit code based on results
        if results["failed"]:
            sys.exit(1)
        else:
            sys.exit(0)

    except Exception as e:
        logger.error("[BATCH][IMPORT] Batch import failed: %r", e)
        print(f"\n❌ Batch import failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

