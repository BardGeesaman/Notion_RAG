#!/usr/bin/env python3
"""
CLI script for pathway enrichment analysis.

Performs pathway enrichment for datasets or signatures and generates reports.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.analysis.pathway_summaries import (
    generate_pathway_aware_dataset_summary,
    generate_pathway_aware_signature_summary,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Perform pathway enrichment analysis"
    )
    
    parser.add_argument(
        "--dataset-id",
        help="Notion page ID of dataset (with dashes)",
    )
    parser.add_argument(
        "--signature-id",
        help="Notion page ID of signature (with dashes)",
    )
    parser.add_argument(
        "--p-value-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for significance (default: 0.05)",
    )
    parser.add_argument(
        "--top-pathways",
        type=int,
        default=10,
        help="Number of top pathways to display (default: 10)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output file path for enrichment report (markdown format)",
    )
    parser.add_argument(
        "--pathway-sources",
        nargs="+",
        choices=["KEGG", "Reactome"],
        help="Pathway sources to use (default: both)",
    )
    
    args = parser.parse_args()
    
    if not args.dataset_id and not args.signature_id:
        parser.error("Either --dataset-id or --signature-id is required")
    
    if args.dataset_id and args.signature_id:
        parser.error("Cannot specify both --dataset-id and --signature-id")
    
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT ANALYSIS")
    print("=" * 80)
    
    try:
        if args.dataset_id:
            print(f"\nDataset ID: {args.dataset_id}")
            print(f"P-value threshold: {args.p_value_threshold}")
            print(f"Top pathways: {args.top_pathways}")
            print()
            
            # Generate pathway-aware summary
            summary = generate_pathway_aware_dataset_summary(
                dataset_page_id=args.dataset_id,
                top_pathways=args.top_pathways,
                p_value_threshold=args.p_value_threshold,
            )
            
            print(f"\n{'=' * 80}")
            print("PATHWAY ENRICHMENT RESULTS")
            print(f"{'=' * 80}\n")
            print(summary)
            
            if args.output:
                args.output.write_text(summary, encoding="utf-8")
                print(f"\n✅ Enrichment report saved to {args.output}")
        
        elif args.signature_id:
            print(f"\nSignature ID: {args.signature_id}")
            print(f"P-value threshold: {args.p_value_threshold}")
            print(f"Top pathways: {args.top_pathways}")
            print()
            
            # Generate pathway-aware summary
            summary = generate_pathway_aware_signature_summary(
                signature_page_id=args.signature_id,
                top_pathways=args.top_pathways,
                p_value_threshold=args.p_value_threshold,
            )
            
            print(f"\n{'=' * 80}")
            print("PATHWAY ENRICHMENT RESULTS")
            print(f"{'=' * 80}\n")
            print(summary)
            
            if args.output:
                args.output.write_text(summary, encoding="utf-8")
                print(f"\n✅ Enrichment report saved to {args.output}")
        
        print(f"\n{'=' * 80}\n")
        
    except Exception as e:
        logger.error("[ANALYSIS][PATHWAY] Enrichment analysis failed: %r", e)
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

