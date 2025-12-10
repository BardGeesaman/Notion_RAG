#!/usr/bin/env python3
"""
CLI script for dataset comparison and clustering.

Compares datasets pairwise and optionally clusters them by similarity.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.analysis.dataset_comparison import (
    cluster_datasets_by_similarity,
    compare_datasets,
    compare_multiple_datasets,
    generate_clustering_report,
    generate_comparison_report,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare datasets and optionally cluster by similarity"
    )
    
    parser.add_argument(
        "--dataset1-id",
        help="Notion page ID of first dataset (with dashes)",  # DEPRECATED: Notion integration removed - these options are no-ops
    )
    parser.add_argument(
        "--dataset2-id",
        help="Notion page ID of second dataset (with dashes)",  # DEPRECATED: Notion integration removed - these options are no-ops
    )
    parser.add_argument(
        "--dataset-ids",
        nargs="+",
        help="List of dataset page IDs for multi-dataset comparison",
    )
    parser.add_argument(
        "--cluster",
        action="store_true",
        help="Cluster datasets by similarity",
    )
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.5,
        help="Similarity threshold for clustering (default: 0.5)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output file path for comparison report (markdown format)",
    )
    
    args = parser.parse_args()
    
    if not args.dataset1_id and not args.dataset_ids:
        parser.error("Either --dataset1-id/--dataset2-id or --dataset-ids is required")
    
    if args.dataset1_id and not args.dataset2_id:
        parser.error("--dataset2-id is required when using --dataset1-id")
    
    if args.dataset1_id and args.dataset_ids:
        parser.error("Cannot specify both --dataset1-id/--dataset2-id and --dataset-ids")
    
    print("\n" + "=" * 80)
    print("DATASET COMPARISON")
    print("=" * 80)
    
    try:
        if args.dataset1_id and args.dataset2_id:
            # Pairwise comparison
            print(f"\nComparing datasets:")
            print(f"  Dataset 1: {args.dataset1_id}")
            print(f"  Dataset 2: {args.dataset2_id}")
            print()
            
            comparison = compare_datasets(
                args.dataset1_id,
                args.dataset2_id,
                use_cache=True,
            )
            
            print(f"\n{'=' * 80}")
            print("COMPARISON RESULTS")
            print(f"{'=' * 80}\n")
            print(f"Overall Similarity: {comparison.overall_similarity:.3f}")
            print(f"Jaccard Similarity: {comparison.jaccard_similarity:.3f}")
            
            if comparison.similarity_by_omics:
                print(f"\nSimilarity by Omics Type:")
                for omics_type, similarity in sorted(comparison.similarity_by_omics.items()):
                    print(f"  {omics_type.capitalize()}: {similarity:.3f}")
            
            # Generate report
            report = generate_comparison_report(comparison)
            
            if args.output:
                args.output.write_text(report, encoding="utf-8")
                print(f"\n✅ Comparison report saved to {args.output}")
            else:
                print(f"\n{'=' * 80}")
                print("FULL COMPARISON REPORT")
                print(f"{'=' * 80}\n")
                print(report)
        
        elif args.dataset_ids:
            # Multi-dataset comparison
            print(f"\nComparing {len(args.dataset_ids)} datasets")
            print()
            
            comparisons = compare_multiple_datasets(
                args.dataset_ids,
                use_cache=True,
            )
            
            print(f"\n{'=' * 80}")
            print("COMPARISON RESULTS")
            print(f"{'=' * 80}\n")
            print(f"Completed {len(comparisons)} pairwise comparisons\n")
            
            # Show top similarities
            comparisons_sorted = sorted(comparisons, key=lambda c: c.overall_similarity, reverse=True)
            print("Top 10 Most Similar Pairs:")
            for i, comp in enumerate(comparisons_sorted[:10], 1):
                print(f"  {i}. {comp.dataset1_name} ↔ {comp.dataset2_name}: {comp.overall_similarity:.3f}")
            
            # Clustering if requested
            if args.cluster:
                print(f"\n{'=' * 80}")
                print("CLUSTERING")
                print(f"{'=' * 80}\n")
                
                clusters = cluster_datasets_by_similarity(
                    comparisons,
                    similarity_threshold=args.similarity_threshold,
                )
                
                print(f"Found {len(clusters)} cluster(s) (threshold={args.similarity_threshold})\n")
                
                for cluster in clusters:
                    print(f"Cluster {cluster.cluster_id + 1}:")
                    print(f"  Datasets: {len(cluster.dataset_ids)}")
                    print(f"  Average similarity: {cluster.average_similarity:.3f}")
                    print(f"  Representative: {cluster.representative_dataset}")
                    print()
                
                # Generate clustering report
                clustering_report = generate_clustering_report(clusters)
                
                if args.output:
                    # Combine comparison and clustering reports
                    comparison_summary = f"# Dataset Comparison Summary\n\n"
                    comparison_summary += f"Compared {len(args.dataset_ids)} datasets in {len(comparisons)} pairwise comparisons.\n\n"
                    full_report = comparison_summary + "\n\n" + clustering_report
                    args.output.write_text(full_report, encoding="utf-8")
                    print(f"✅ Clustering report saved to {args.output}")
                else:
                    print(f"\n{'=' * 80}")
                    print("FULL CLUSTERING REPORT")
                    print(f"{'=' * 80}\n")
                    print(clustering_report)
            else:
                # Just comparison report
                if args.output:
                    report = f"# Dataset Comparison Summary\n\n"
                    report += f"Compared {len(args.dataset_ids)} datasets in {len(comparisons)} pairwise comparisons.\n\n"
                    report += "## Pairwise Comparisons\n\n"
                    for comp in comparisons_sorted:
                        report += f"### {comp.dataset1_name} vs {comp.dataset2_name}\n\n"
                        report += f"- **Similarity**: {comp.overall_similarity:.3f}\n"
                        report += f"- **Jaccard**: {comp.jaccard_similarity:.3f}\n\n"
                    args.output.write_text(report, encoding="utf-8")
                    print(f"\n✅ Comparison report saved to {args.output}")
        
        print(f"\n{'=' * 80}\n")
        
    except Exception as e:
        logger.error("[ANALYSIS][DATASET-COMP] Comparison failed: %r", e)
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
