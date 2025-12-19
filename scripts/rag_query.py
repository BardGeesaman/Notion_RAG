#!/usr/bin/env python3
# scripts/rag_query.py

import argparse
import json
from dataclasses import asdict

from amprenta_rag.query.rag_query_engine import (MatchSummary, RAGQueryResult,
                                                 query_rag,
                                                 signature_similarity_query)
from amprenta_rag.query.cross_omics_reasoning import (
    cross_omics_dataset_summary,
    cross_omics_feature_summary,
    cross_omics_program_summary,
    cross_omics_signature_summary,
)


def format_match_for_cli(m: MatchSummary) -> str:
    tag_str = ", ".join(m.tags)
    return (
        f"[score={m.score:.3f}]\n"
        f"[{m.source}] {m.title}  [tags: {tag_str}]\n"
        f"  {m.snippet}"
    )


def print_result_cli(
    result: RAGQueryResult,
    *,
    show_context: bool,
    show_answer: bool,
) -> None:
    print(f"\n🔎 Query: {result.query}\n")

    if not result.matches:
        print(result.answer)
        return

    print("📌 Top matches:")
    for m in result.filtered_matches or result.matches:
        print()
        print(format_match_for_cli(m))

    if show_context and result.context_chunks:
        print("\n📚 Context chunks used:")
        for i, chunk in enumerate(result.context_chunks, start=1):
            print(f"\n--- Chunk {i} ---")
            print(chunk[:2000])

    if show_answer:
        print("\n🧠 Synthesized answer:\n")
        print(result.answer)
        print()


def main() -> None:
    parser = argparse.ArgumentParser(description="Query the Amprenta RAG engine.")
    parser.add_argument(
        "--query",
        required=False,
        help="User query text. Not required if --signature-score is provided.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=10,
        help="Number of results to retrieve from Pinecone.",
    )
    parser.add_argument(
        "--source-type",
        dest="source_types",
        nargs="*",
        choices=["Literature", "Email", "Experiment", "Dataset", "Signature"],
        default=["Literature"],
        help=(
            "Optional source type filter(s). "
            "You can pass multiple values, e.g. "
            "--source-type Literature Experiment Dataset Signature. "
            "Defaults to ['Literature']."
        ),
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="Optional tag substring filter.",
    )
    parser.add_argument(
        "--show-context",
        action="store_true",
        help="Print the context chunks used to answer the query.",
    )
    parser.add_argument(
        "--no-answer",
        action="store_true",
        help="Do not call OpenAI to synthesize an answer; show matches (and optionally context) only.",
    )
    parser.add_argument(
        "--raw-json",
        action="store_true",
        help="Output the full RAGQueryResult as JSON instead of pretty CLI output.",
    )

    parser.add_argument(
        "--disease",
        default=None,
        help="Optional disease filter (e.g. 'ALS', 'AD'). Matches 'diseases' metadata.",
    )
    parser.add_argument(
        "--target",
        default=None,
        help="Optional molecular target filter (e.g. 'SPTLC1', 'DEGS1'). Matches 'targets' metadata.",
    )
    parser.add_argument(
        "--lipid",
        default=None,
        help="Optional lipid filter (canonical ID or raw label). Matches 'lipids' or 'lipids_raw' metadata.",
    )
    parser.add_argument(
        "--signature",
        default=None,
        help="Optional lipid signature filter (Short ID). Matches 'lipid_signatures' metadata.",
    )
    parser.add_argument(
        "--signature-score",
        dest="signature_score_dataset_id",
        default=None,
        help=(
            "Dataset page ID for signature similarity query. "
            "When provided, computes and returns top matching signatures for the dataset. "
            "Ignores --query if provided."
        ),
    )
    parser.add_argument(
        "--cross-omics-program",
        dest="cross_omics_program_id",
        default=None,
        help="Program page ID for cross-omics summary. Ignores --query if provided.",
    )
    parser.add_argument(
        "--cross-omics-signature",
        dest="cross_omics_signature_id",
        default=None,
        help="Signature page ID for cross-omics summary. Ignores --query if provided.",
    )
    parser.add_argument(
        "--cross-omics-feature",
        dest="cross_omics_feature",
        default=None,
        help="Feature for cross-omics summary in format 'feature_type:feature_name' (e.g., 'gene:TP53', 'lipid:Cer(d18:1/16:0)'). Ignores --query if provided.",
    )
    parser.add_argument(
        "--cross-omics-dataset",
        dest="cross_omics_dataset_id",
        default=None,
        help="Dataset page ID for cross-omics summary. Ignores --query if provided.",
    )

    args = parser.parse_args()

    # Handle signature similarity query
    if args.signature_score_dataset_id:
        results = signature_similarity_query(
            dataset_page_id=args.signature_score_dataset_id,
            top_k=args.top_k,
        )

        if not results:
            print(
                f"\n🔎 Signature similarity for dataset {args.signature_score_dataset_id}\n"
            )
            print("No matching signatures found.")
            return

        print(
            f"\n🔎 Signature similarity for dataset {args.signature_score_dataset_id}\n"
        )

        for i, result in enumerate(results, start=1):
            signature_name = result.get("signature_name", "Unknown")
            score = result.get("score", 0.0)
            overlap = result.get("overlap_fraction", 0.0)
            matched = result.get("matched_components", [])
            missing = result.get("missing_components", [])
            conflicting = result.get("conflicting_components", [])

            # Calculate component counts
            total_components = len(matched) + len(missing)
            matched_count = len(matched)

            print(f"{i}. {signature_name}")
            print(f"   - Score: {score:.3f}")
            print(
                f"   - Overlap: {overlap:.2f} ({matched_count}/{total_components} components)"
            )

            if matched:
                matched_str = ", ".join(matched[:10])
                if len(matched) > 10:
                    matched_str += f", ... ({len(matched) - 10} more)"
                print(f"   - Matched: {matched_str}")
            else:
                print("   - Matched: (none)")

            if missing:
                missing_str = ", ".join(missing[:5])
                if len(missing) > 5:
                    missing_str += f", ... ({len(missing) - 5} more)"
                print(f"   - Missing: {missing_str}")
            else:
                print("   - Missing: (none)")

            if conflicting:
                conflicting_str = ", ".join(conflicting[:5])
                if len(conflicting) > 5:
                    conflicting_str += f", ... ({len(conflicting) - 5} more)"
                print(f"   - Conflicting: {conflicting_str}")
            else:
                print("   - Conflicting: (none)")

            # Show additional metadata if available
            if result.get("signature_short_id"):
                print(f"   - Short ID: {result['signature_short_id']}")
            if result.get("disease"):
                print(f"   - Disease: {', '.join(result['disease'])}")
            if result.get("matrix"):
                print(f"   - Matrix: {', '.join(result['matrix'])}")

            print()

        return

    # Handle cross-omics program summary
    if args.cross_omics_program_id:
        summary = cross_omics_program_summary(
            program_page_id=args.cross_omics_program_id,
            top_k_per_omics=args.top_k,
        )
        print(f"\n🔬 Cross-Omics Summary for Program: {args.cross_omics_program_id}\n")
        print(summary)
        print()
        return

    # Handle cross-omics signature summary
    if args.cross_omics_signature_id:
        summary = cross_omics_signature_summary(
            signature_page_id=args.cross_omics_signature_id,
            top_k_datasets=args.top_k,
            top_k_chunks=args.top_k * 5,  # More chunks for comprehensive summary
        )
        print(f"\n🔬 Cross-Omics Summary for Signature: {args.cross_omics_signature_id}\n")
        print(summary)
        print()
        return

    # Handle cross-omics feature summary
    if args.cross_omics_feature:
        if ":" not in args.cross_omics_feature:
            print(
                "Error: --cross-omics-feature must be in format 'feature_type:feature_name' "
                "(e.g., 'gene:TP53', 'lipid:Cer(d18:1/16:0)')"
            )
            return

        feature_type, feature_name = args.cross_omics_feature.split(":", 1)
        feature_type = feature_type.lower().strip()
        feature_name = feature_name.strip()

        if feature_type not in ["gene", "protein", "metabolite", "lipid"]:
            print(
                f"Error: Invalid feature_type '{feature_type}'. Must be one of: gene, protein, metabolite, lipid"
            )
            return

        summary = cross_omics_feature_summary(
            feature_name=feature_name,
            feature_type=feature_type,
            top_k_datasets=args.top_k,
            top_k_chunks=args.top_k * 5,
        )
        print(f"\n🔬 Cross-Omics Summary for {feature_type.title()}: {feature_name}\n")
        print(summary)
        print()
        return

    # Handle cross-omics dataset summary
    if args.cross_omics_dataset_id:
        summary = cross_omics_dataset_summary(
            dataset_page_id=args.cross_omics_dataset_id,
            top_k_chunks=args.top_k * 5,
        )
        print(f"\n🔬 Cross-Omics Summary for Dataset: {args.cross_omics_dataset_id}\n")
        print(summary)
        print()
        return

    # Require --query if no special mode is specified
    if not args.query:
        parser.error("--query is required unless one of --signature-score, --cross-omics-program, --cross-omics-signature, --cross-omics-feature, or --cross-omics-dataset is provided.")

    result = query_rag(
        user_query=args.query,
        top_k=args.top_k,
        source_types=args.source_types,
        tag=args.tag,
        generate_answer=not args.no_answer,
        disease=args.disease,
        target=args.target,
        lipid=args.lipid,
        signature=args.signature,
    )

    if args.raw_json:
        print(json.dumps(asdict(result), indent=2))
        return

    print_result_cli(
        result,
        show_context=args.show_context,
        show_answer=not args.no_answer,
    )


if __name__ == "__main__":
    main()
