#!/usr/bin/env python3
"""
Comprehensive script to discover and import all types of omics data from online repositories.

This script provides a unified interface to:
1. Discover studies from all repositories (GEO, PRIDE, MetaboLights, MW)
2. Filter by omics type, disease, organism, etc.
3. Import selected studies into Postgres
4. Optionally trigger ingestion to Pinecone

Usage:
    # Discover all omics types for a disease
    python scripts/import_all_omics_repositories.py \
        --disease "ALS" \
        --discover-only \
        --output discovered_studies.json
    
    # Import all discovered studies
    python scripts/import_all_omics_repositories.py \
        --input discovered_studies.json \
        --import \
        --ingest
    
    # Discover and import in one step
    python scripts/import_all_omics_repositories.py \
        --keywords "Alzheimer" "amyloid" \
        --omics-type metabolomics \
        --import \
        --ingest
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional
from tqdm import tqdm

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repositories.discovery import (
    discover_studies,
    list_available_repositories,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def discover_all_repositories(
    keywords: Optional[List[str]] = None,
    omics_type: Optional[str] = None,
    disease: Optional[str] = None,
    organism: Optional[str] = None,
    sample_type: Optional[str] = None,
    repository: Optional[str] = None,
    max_results: int = 50,
) -> Dict[str, List[str]]:
    """
    Discover studies from all or selected repositories.
    
    Args:
        keywords: Search keywords
        omics_type: Filter by omics type
        disease: Filter by disease
        organism: Filter by organism
        sample_type: Filter by sample type
        repository: Specific repository to search (or None for all)
        max_results: Maximum results per repository
        
    Returns:
        Dictionary mapping repository names to lists of study IDs
    """
    # Map omics types to repositories
    repository_map = {
        "transcriptomics": ["GEO"],
        "proteomics": ["PRIDE"],
        "metabolomics": ["MetaboLights", "MW_METABOLOMICS"],
        "lipidomics": ["MW_LIPIDOMICS"],
    }
    
    # Determine which repositories to search
    repositories_to_search = []
    if repository:
        repositories_to_search = [repository]
    elif omics_type:
        repositories_to_search = repository_map.get(omics_type, list_available_repositories())
    else:
        # Search all repositories
        repositories_to_search = None  # None = all repositories
    
    # Build filters
    filters: Dict[str, str] = {}
    if disease:
        filters["disease"] = disease
    if organism:
        filters["organism"] = organism
    if sample_type:
        filters["sample_type"] = sample_type
    
    logger.info(
        "[IMPORT][ALL] Discovering studies with keywords=%s, omics_type=%s, filters=%s",
        keywords,
        omics_type,
        filters,
    )
    
    # Discover studies
    # If multiple repositories, search all and combine results
    all_results: Dict[str, List[str]] = {}
    
    if repositories_to_search:
        # Search specific repositories
        for repo_name in repositories_to_search:
            repo_results = discover_studies(
                keywords=keywords or [],
                repository=repo_name,
                omics_type=omics_type,
                filters=filters if filters else None,
                max_results=max_results,
            )
            all_results.update(repo_results)
    else:
        # Search all repositories
        all_results = discover_studies(
            keywords=keywords or [],
            repository=None,
            omics_type=omics_type,
            filters=filters if filters else None,
            max_results=max_results,
        )
    
    results = all_results
    
    return results


def save_discovery_results(results: Dict[str, List[str]], output_file: Path) -> None:
    """Save discovery results to JSON file."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Format for batch import script
    output_data = {
        "results": results,
        "summary": {
            "total_repositories": len(results),
            "total_studies": sum(len(studies) for studies in results.values()),
            "by_repository": {
                repo: len(studies) for repo, studies in results.items()
            },
        },
    }
    
    output_file.write_text(json.dumps(output_data, indent=2))
    logger.info(
        "[IMPORT][ALL] Saved discovery results to %s: %d repositories, %d total studies",
        output_file,
        output_data["summary"]["total_repositories"],
        output_data["summary"]["total_studies"],
    )


def print_discovery_summary(results: Dict[str, List[str]]) -> None:
    """Print a summary of discovered studies."""
    total_studies = sum(len(studies) for studies in results.values())
    
    print("\n" + "=" * 60)
    print("DISCOVERY SUMMARY")
    print("=" * 60)
    print(f"Total repositories: {len(results)}")
    print(f"Total studies found: {total_studies}")
    print("\nBreakdown by repository:")
    
    for repo, studies in sorted(results.items()):
        print(f"  {repo}: {len(studies)} studies")
        if studies:
            print(f"    Examples: {', '.join(studies[:5])}")
            if len(studies) > 5:
                print(f"    ... and {len(studies) - 5} more")
    
    print("=" * 60 + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Discover and import omics studies from all repositories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Discover studies for a disease across all repositories
  %(prog)s --disease "ALS" --discover-only --output results.json
  
  # Discover transcriptomics studies
  %(prog)s --omics-type transcriptomics --keywords "Alzheimer" --discover-only
  
  # Import from discovery file
  %(prog)s --input results.json --import --ingest
  
  # Discover and import in one step
  %(prog)s --disease "Parkinson" --import --ingest
  
  # Discover from specific repository
  %(prog)s --repository MW_LIPIDOMICS --keywords "ceramide" --discover-only
        """,
    )
    
    # Discovery arguments
    discovery_group = parser.add_argument_group("Discovery Options")
    discovery_group.add_argument(
        "--keywords",
        nargs="+",
        help="Search keywords (e.g., 'ALS' 'Alzheimer')",
    )
    discovery_group.add_argument(
        "--omics-type",
        choices=["transcriptomics", "proteomics", "metabolomics", "lipidomics"],
        help="Filter by omics type",
    )
    discovery_group.add_argument(
        "--repository",
        help="Search specific repository (GEO, PRIDE, MetaboLights, MW, MW_LIPIDOMICS, MW_METABOLOMICS). "
             "If not specified, searches all repositories.",
    )
    discovery_group.add_argument(
        "--disease",
        help="Filter by disease (e.g., 'ALS', 'Alzheimer', 'Parkinson')",
    )
    discovery_group.add_argument(
        "--organism",
        help="Filter by organism (e.g., 'Homo sapiens', 'Mus musculus')",
    )
    discovery_group.add_argument(
        "--sample-type",
        help="Filter by sample type (e.g., 'CSF', 'plasma', 'tissue')",
    )
    discovery_group.add_argument(
        "--max-results",
        type=int,
        default=50,
        help="Maximum results per repository (default: 50)",
    )
    
    # Import arguments
    import_group = parser.add_argument_group("Import Options")
    import_group.add_argument(
        "--input",
        type=Path,
        help="JSON file with study list (from discovery or manual)",
    )
    import_group.add_argument(
        "--discover-only",
        action="store_true",
        help="Only discover studies, don't import",
    )
    import_group.add_argument(
        "--import",
        dest="do_import",
        action="store_true",
        help="Import discovered/loaded studies",
    )
    import_group.add_argument(
        "--output",
        type=Path,
        help="Output JSON file for discovery results",
    )
    import_group.add_argument(
        "--ingest",
        action="store_true",
        help="Trigger ingestion to Pinecone after import",
    )
    import_group.add_argument(
        "--create-notion",
        action="store_true",
        help="Also create/update Notion pages (optional)",
    )
    import_group.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop batch import on first error",
    )
    
    # Other options
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
    
    # Discovery phase
    discovery_results: Optional[Dict[str, List[str]]] = None
    
    if args.input:
        # Load from file
        logger.info("[IMPORT][ALL] Loading studies from %s", args.input)
        try:
            data = json.loads(args.input.read_text())
            if "results" in data:
                discovery_results = data["results"]
            elif isinstance(data, dict):
                discovery_results = data
            else:
                logger.error("[IMPORT][ALL] Invalid input file format")
                return
        except Exception as e:
            logger.error("[IMPORT][ALL] Error loading input file: %r", e)
            return
    elif args.keywords or args.disease or args.omics_type or args.repository:
        # Discover studies
        logger.info("[IMPORT][ALL] Starting discovery phase...")
        discovery_results = discover_all_repositories(
            keywords=args.keywords,
            omics_type=args.omics_type,
            disease=args.disease,
            organism=args.organism,
            sample_type=args.sample_type,
            repository=args.repository,
            max_results=args.max_results,
        )
        
        # Print summary
        print_discovery_summary(discovery_results)
        
        # Save to file if requested
        if args.output:
            save_discovery_results(discovery_results, args.output)
    
    # Check if we should proceed with import
    if args.discover_only:
        if not args.output:
            logger.warning(
                "[IMPORT][ALL] Discovery-only mode but no --output specified. "
                "Results will not be saved."
            )
        logger.info("[IMPORT][ALL] Discovery complete. Use --input to import later.")
        return
    
    if not args.do_import:
        logger.info("[IMPORT][ALL] No import requested. Use --import to import studies.")
        return
    
    if not discovery_results:
        logger.error("[IMPORT][ALL] No studies to import. Run discovery first or provide --input file.")
        return
    
    # Import phase
    logger.info("[IMPORT][ALL] Starting import phase...")
    
    # Prepare batch import arguments
    # Convert discovery results to batch import format
    studies_to_import = []
    for repo, study_ids in discovery_results.items():
        for study_id in study_ids:
            studies_to_import.append({
                "study_id": study_id,
                "repository": repo,
            })
    
    logger.info(
        "[IMPORT][ALL] Importing %d studies from %d repositories",
        len(studies_to_import),
        len(discovery_results),
    )
    
    # Create temporary JSON file for batch import
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tmp_file:
        tmp_path = Path(tmp_file.name)
        json.dump({"results": discovery_results}, tmp_file, indent=2)
    
    try:
        # Call batch import script
        batch_args = [
            str(tmp_path),
            "--create-postgres",
        ]
        
        if args.ingest:
            batch_args.append("--ingest")
        
        if args.create_notion:
            batch_args.append("--create-notion")
        
        if args.stop_on_error:
            batch_args.append("--stop-on-error")
        
        logger.info("[IMPORT][ALL] Running batch import with args: %s", batch_args)
        
        # Import batch_import_main function and call it with our args
        from scripts.batch_import_repository_studies import load_studies_from_json
        from scripts.harvest_repository_study import harvest_study
        
        studies = load_studies_from_json(tmp_path)
        
        # Process each study
        success_count = 0
        error_count = 0
        
        for study in tqdm(studies, desc="Importing studies"):
            try:
                result = harvest_study(
                    study_id=study["study_id"],
                    repository=study["repository"],
                    create_postgres=True,
                    create_notion=args.create_notion,
                    ingest=args.ingest,
                )
                
                if result:
                    success_count += 1
                    logger.info(
                        "[IMPORT][ALL] ✓ Successfully imported %s:%s",
                        study["repository"],
                        study["study_id"],
                    )
                else:
                    error_count += 1
                    logger.warning(
                        "[IMPORT][ALL] ✗ Failed to import %s:%s",
                        study["repository"],
                        study["study_id"],
                    )
                
                if args.stop_on_error and not result:
                    logger.error("[IMPORT][ALL] Stopping due to error (--stop-on-error)")
                    break
                    
            except Exception as e:
                error_count += 1
                logger.error(
                    "[IMPORT][ALL] ✗ Error importing %s:%s: %r",
                    study["repository"],
                    study["study_id"],
                    e,
                )
                
                if args.stop_on_error:
                    logger.error("[IMPORT][ALL] Stopping due to error (--stop-on-error)")
                    break
        
        logger.info(
            "[IMPORT][ALL] Import complete: %d succeeded, %d failed",
            success_count,
            error_count,
        )
        
    finally:
        # Clean up temporary file
        tmp_path.unlink()


if __name__ == "__main__":
    main()

