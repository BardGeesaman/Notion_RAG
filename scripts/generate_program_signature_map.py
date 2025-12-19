#!/usr/bin/env python3
"""
Generate program-signature maps for dashboard intelligence.

Computes Program × Signature scoring matrices, coverage analysis,
and convergence indicators.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.analysis.program_signature_maps import (
    generate_program_map_report,
    generate_program_signature_map,
    update_notion_with_program_map,  # DEPRECATED: Notion integration removed - these options are no-ops
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate program-signature maps"
    )

    parser.add_argument(
        "--program-id",
        required=True,
        help="Notion page ID of program (with dashes)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of top signatures to include (default: 10)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output file path for map report (markdown format)",
    )
    parser.add_argument(
        "--update-notion",
        action="store_true",
        help="Update Notion program page with map summary",  # DEPRECATED: Notion integration removed - these options are no-ops
    )
    parser.add_argument(
        "--include-all",
        action="store_true",
        help="Include all signatures in report (not just top N)",
    )

    args = parser.parse_args()

    print("\n" + "=" * 80)
    print("PROGRAM-SIGNATURE MAP GENERATION")
    print("=" * 80)
    print(f"\nProgram ID: {args.program_id}")
    print(f"Top N signatures: {args.top_n}")
    print()

    try:
        # Generate program-signature map
        program_map = generate_program_signature_map(
            program_page_id=args.program_id,
            top_n=args.top_n,
            use_cache=True,
        )

        print(f"\n{'=' * 80}")
        print("PROGRAM-SIGNATURE MAP RESULTS")
        print(f"{'=' * 80}\n")
        print(f"Program: {program_map.program_name}")

        if program_map.omics_coverage:
            print(f"\nOmics Coverage:")
            print(f"  Total Datasets: {program_map.omics_coverage.total_datasets}")
            for omics_type, count in sorted(program_map.omics_coverage.datasets_by_omics.items()):
                features_count = program_map.omics_coverage.features_by_omics.get(omics_type, 0)
                print(f"  {omics_type.capitalize()}: {count} dataset(s), {features_count} features")

        print(f"\nSignature Scores: {len(program_map.signature_scores)} matching signature(s)")

        if program_map.top_signatures:
            print(f"\nTop {len(program_map.top_signatures)} Signatures:")
            for i, score in enumerate(program_map.top_signatures, 1):
                print(f"  {i}. {score.signature_name}")
                print(f"     Score: {score.overall_score:.3f}")
                print(f"     Coverage: {score.coverage_fraction:.1%} ({len(score.matching_datasets)} datasets)")
                if score.score_by_omics:
                    omics_str = ", ".join(f"{k}: {v:.3f}" for k, v in sorted(score.score_by_omics.items()))
                    print(f"     By Omics: {omics_str}")

        if program_map.convergence_indicators:
            print(f"\nCross-Omics Convergence:")
            print(f"  Multi-Omics Signatures: {program_map.convergence_indicators.get('multi_omics_signature_count', 0)}")
            print(f"  Convergence Fraction: {program_map.convergence_indicators.get('convergence_fraction', 0.0):.3f}")
            print(f"  Avg Omics per Signature: {program_map.convergence_indicators.get('avg_omics_per_signature', 0.0):.2f}")

        # Generate report
        report = generate_program_map_report(
            program_map,
            include_all_signatures=args.include_all,
        )

        if args.output:
            args.output.write_text(report, encoding="utf-8")
            print(f"\n✅ Map report saved to {args.output}")
        else:
            print(f"\n{'=' * 80}")
            print("FULL PROGRAM-SIGNATURE MAP REPORT")
            print(f"{'=' * 80}\n")
            print(report)

        # Update Notion if requested
        if args.update_notion:
            print("\n📥 Updating Notion program page...")
            update_notion_with_program_map(program_map)
            print("✅ Notion page updated")

        print(f"\n{'=' * 80}\n")

    except Exception as e:
        logger.error("[ANALYSIS][PROGRAM-MAPS] Map generation failed: %r", e)
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

