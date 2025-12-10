#!/usr/bin/env python3
"""
Discover signature candidates from datasets.

Analyzes multiple datasets to find recurring patterns and generate
candidate signatures based on feature co-occurrence and direction consistency.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.discovery import (
    discover_signature_candidates,
    export_candidates_to_tsv,
)

logger = get_logger(__name__)


def list_datasets():
    """List all datasets from Postgres."""
    db = next(get_db())
    try:
        datasets_query = db.query(DatasetModel).all()
        return [{"id": str(dataset.id), "name": dataset.name} for dataset in datasets_query]
    except Exception as e:
        logger.error("[SIG-DISC] Error listing datasets: %r", e)
        return []
    finally:
        db.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Discover signature candidates from datasets"
    )
    parser.add_argument(
        "--dataset-ids",
        nargs="+",
        help="Specific dataset page IDs to analyze (with dashes)",
    )
    parser.add_argument(
        "--all-datasets",
        action="store_true",
        help="Analyze all datasets in the Experimental Data Assets database",
    )
    parser.add_argument(
        "--min-support",
        type=int,
        default=3,
        help="Minimum number of datasets that must contain the pattern (default: 3)",
    )
    parser.add_argument(
        "--min-features",
        type=int,
        default=3,
        help="Minimum number of features in a signature (default: 3)",
    )
    parser.add_argument(
        "--max-features",
        type=int,
        default=20,
        help="Maximum number of features in a signature (default: 20)",
    )
    parser.add_argument(
        "--min-co-occurrence",
        type=int,
        default=2,
        help="Minimum co-occurrence count for feature pairs (default: 2)",
    )
    parser.add_argument(
        "--min-confidence",
        type=float,
        default=0.5,
        help="Minimum confidence score to include candidate (default: 0.5)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("discovered_signatures.tsv"),
        help="Output file path (default: discovered_signatures.tsv). Extension determines format (.tsv or .json)",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Use parallel processing for feature extraction",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=5,
        help="Maximum number of parallel workers (default: 5)",
    )
    parser.add_argument(
        "--ingest",
        action="store_true",
        help="Ingest discovered signatures into Notion after discovery",  # DEPRECATED: Notion integration removed - these options are no-ops
    )

    args = parser.parse_args()

    # Get dataset IDs
    if args.dataset_ids:
        dataset_ids = args.dataset_ids
    elif args.all_datasets:
        logger.info("[SIG-DISC] Fetching all datasets...")
        try:
            all_datasets = list_datasets()
            dataset_ids = [d["id"] for d in all_datasets]
            logger.info("[SIG-DISC] Found %d datasets", len(dataset_ids))
        except Exception as e:
            logger.error("[SIG-DISC] Error listing datasets: %r", e)
            print(f"\n❌ Error: Could not list datasets: {e}\n", file=sys.stderr)
            sys.exit(1)
    else:
        parser.error("Either --dataset-ids or --all-datasets must be specified")

    if not dataset_ids:
        print("\n❌ No datasets to analyze\n", file=sys.stderr)
        sys.exit(1)

    print("\n" + "=" * 80)
    print("SIGNATURE DISCOVERY")
    print("=" * 80)
    print(f"\nAnalyzing {len(dataset_ids)} dataset(s)...")
    print(f"Parameters:")
    print(f"  Min support: {args.min_support}")
    print(f"  Min features: {args.min_features}")
    print(f"  Max features: {args.max_features}")
    print(f"  Min co-occurrence: {args.min_co_occurrence}")
    print(f"  Min confidence: {args.min_confidence}")

    try:
        # Discover candidates
        candidates = discover_signature_candidates(
            dataset_ids=dataset_ids,
            min_support=args.min_support,
            min_features=args.min_features,
            max_features=args.max_features,
            min_co_occurrence=args.min_co_occurrence,
            min_confidence=args.min_confidence,
        )

        if not candidates:
            print("\n⚠️  No signature candidates found with the given parameters.")
            print("   Try lowering --min-support or --min-confidence")
            return

        # Display results
        print(f"\n{'=' * 80}")
        print(f"DISCOVERED {len(candidates)} SIGNATURE CANDIDATE(S)")
        print(f"{'=' * 80}\n")

        for i, candidate in enumerate(candidates, 1):
            print(f"{i}. {candidate.name}")
            print(f"   Features: {len(candidate.features)}")
            print(f"   Types: {', '.join(sorted(candidate.feature_types))}")
            print(f"   Support: {candidate.support_count} dataset(s)")
            print(f"   Confidence: {candidate.confidence:.3f}")
            print(f"   Co-occurrence: {candidate.co_occurrence_score:.3f}")
            print(f"   Direction consistency: {candidate.direction_consistency:.3f}")
            print(f"   Top features:")
            for feat in candidate.features[:5]:
                direction_str = f" {feat.direction}" if feat.direction else ""
                print(f"     - {feat.feature_name} ({feat.feature_type}){direction_str}")
            if len(candidate.features) > 5:
                print(f"     ... and {len(candidate.features) - 5} more")
            print()

        # Export to file (format based on extension)
        output_path_str = str(args.output)
        if output_path_str.lower().endswith(".json"):
            from amprenta_rag.signatures.discovery import export_candidates_to_json
            export_candidates_to_json(candidates, output_path_str)
        else:
            export_candidates_to_tsv(candidates, output_path_str)
        print(f"✅ Exported {len(candidates)} candidates to {args.output}")

        # Optionally ingest into Notion
        if args.ingest:
            print("\n📥 Ingesting candidates into Notion...")
            from amprenta_rag.ingestion.signature_ingestion import (
                ingest_signature_from_file,
            )

            ingested_count = 0
            for candidate in candidates:
                # Create a temporary TSV file for this candidate
                import tempfile

                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".tsv", delete=False
                ) as tmp_file:
                    # Write candidate to TSV
                    tmp_file.write(
                        "feature_type\tfeature_name\tdirection\tweight\n"
                    )
                    for feat in candidate.features:
                        weight = candidate.confidence
                        tmp_file.write(
                            f"{feat.feature_type}\t"
                            f"{feat.feature_name}\t"
                            f"{feat.direction or ''}\t"
                            f"{weight:.3f}\n"
                        )
                    tmp_file_path = tmp_file.name

                try:
                    # Ingest the signature
                    inferred_metadata = {
                        "name": candidate.name,
                        "signature_type": "Discovered",
                        "description": (
                            f"Automatically discovered signature with "
                            f"{len(candidate.features)} features, "
                            f"supported by {candidate.support_count} dataset(s), "
                            f"confidence: {candidate.confidence:.3f}"
                        ),
                    }
                    ingest_signature_from_file(
                        tmp_file_path,
                        inferred_metadata=inferred_metadata,
                    )
                    ingested_count += 1
                    print(f"  ✅ Ingested: {candidate.name}")
                except Exception as e:
                    logger.warning(
                        "[SIG-DISC] Error ingesting candidate %s: %r",
                        candidate.name,
                        e,
                    )
                    print(f"  ⚠️  Failed to ingest: {candidate.name}")
                finally:
                    # Clean up temp file
                    Path(tmp_file_path).unlink()

            print(f"\n✅ Ingested {ingested_count}/{len(candidates)} candidates")

        print(f"\n{'=' * 80}\n")

    except Exception as e:
        logger.error("[SIG-DISC] Discovery failed: %r", e)
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

