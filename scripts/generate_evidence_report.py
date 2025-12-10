#!/usr/bin/env python3
"""
CLI script for generating evidence reports.

Generates comprehensive cross-omics evidence reports for Programs, Experiments,
Datasets, Signatures, and Features.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.reporting.evidence_report import (
    format_evidence_report,
    generate_dataset_evidence_report,
    generate_experiment_evidence_report,
    generate_feature_evidence_report,
    generate_program_evidence_report,
    generate_signature_evidence_report,
    write_evidence_report_to_notion,  # DEPRECATED: Notion integration removed - these options are no-ops
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate evidence reports for Notion entities"
    )
    
    # Entity type selection
    entity_group = parser.add_mutually_exclusive_group(required=True)
    entity_group.add_argument(
        "--program-id",
        help="Notion page ID of program (with dashes)",
    )
    entity_group.add_argument(
        "--experiment-id",
        help="Notion page ID of experiment (with dashes)",
    )
    entity_group.add_argument(
        "--dataset-id",
        help="Notion page ID of dataset (with dashes)",
    )
    entity_group.add_argument(
        "--signature-id",
        help="Notion page ID of signature (with dashes)",
    )
    entity_group.add_argument(
        "--feature",
        help="Feature name",
    )
    entity_group.add_argument(
        "--feature-type",
        choices=["gene", "protein", "metabolite", "lipid"],
        help="Feature type (required with --feature)",
    )
    
    # Options
    parser.add_argument(
        "--top-k-per-omics",
        type=int,
        default=20,
        help="Number of top chunks per omics type (for programs, default: 20)",
    )
    parser.add_argument(
        "--top-k-chunks",
        type=int,
        default=100,
        help="Number of top chunks to retrieve (default: 100)",
    )
    parser.add_argument(
        "--top-k-datasets",
        type=int,
        default=20,
        help="Number of top datasets to include (for signatures/features, default: 20)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output file path for evidence report (markdown format)",
    )
    parser.add_argument(
        "--write-to-notion",
        action="store_true",
        help="Write report to Notion page (property: 'Evidence Report')",  # DEPRECATED: Notion integration removed - these options are no-ops
    )
    parser.add_argument(
        "--notion-property",
        default="Evidence Report",
        help="Notion property name to write to (default: 'Evidence Report')",
    )
    
    args = parser.parse_args()
    
    if args.feature and not args.feature_type:
        parser.error("--feature-type is required when using --feature")
    
    print("\n" + "=" * 80)
    print("EVIDENCE REPORT GENERATION")
    print("=" * 80)
    
    try:
        report = None
        
        if args.program_id:
            print(f"\nGenerating evidence report for program: {args.program_id}")
            report = generate_program_evidence_report(
                program_page_id=args.program_id,
                top_k_per_omics=args.top_k_per_omics,
            )
        
        elif args.experiment_id:
            print(f"\nGenerating evidence report for experiment: {args.experiment_id}")
            report = generate_experiment_evidence_report(
                experiment_page_id=args.experiment_id,
                top_k_chunks=args.top_k_chunks,
            )
        
        elif args.dataset_id:
            print(f"\nGenerating evidence report for dataset: {args.dataset_id}")
            report = generate_dataset_evidence_report(
                dataset_page_id=args.dataset_id,
                top_k_chunks=args.top_k_chunks,
            )
        
        elif args.signature_id:
            print(f"\nGenerating evidence report for signature: {args.signature_id}")
            report = generate_signature_evidence_report(
                signature_page_id=args.signature_id,
                top_k_datasets=args.top_k_datasets,
                top_k_chunks=args.top_k_chunks,
            )
        
        elif args.feature:
            print(f"\nGenerating evidence report for {args.feature_type} feature: {args.feature}")
            report = generate_feature_evidence_report(
                feature_name=args.feature,
                feature_type=args.feature_type,
                top_k_datasets=args.top_k_datasets,
                top_k_chunks=args.top_k_chunks,
            )
        
        if not report:
            print("\n❌ Error: Could not generate report", file=sys.stderr)
            sys.exit(1)
        
        # Format report
        formatted_report = format_evidence_report(report, include_metadata=True)
        
        print(f"\n{'=' * 80}")
        print("EVIDENCE REPORT")
        print(f"{'=' * 80}\n")
        print(formatted_report)
        
        # Write to file if requested
        if args.output:
            args.output.write_text(formatted_report, encoding="utf-8")
            print(f"\n✅ Evidence report saved to {args.output}")
        
        # Write to Notion if requested
        if args.write_to_notion:
            print(f"\n📥 Writing report to Notion...")
            success = write_evidence_report_to_notion(
                report=report,
                property_name=args.notion_property,
            )
            if success:
                print("✅ Evidence report written to Notion")
            else:
                print("⚠️  Warning: Could not write report to Notion (check logs)")
        
        print(f"\n{'=' * 80}\n")
        
    except Exception as e:
        logger.error("[REPORTING][EVIDENCE] Report generation failed: %r", e)
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
