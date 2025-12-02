#!/usr/bin/env python3
# scripts/rag_query.py

import argparse
import json
from dataclasses import asdict

from amprenta_rag.query.rag_query_engine import query_rag, MatchSummary, RAGQueryResult


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
    parser.add_argument("--query", required=True, help="User query text.")
    parser.add_argument("--top-k", type=int, default=10, help="Number of results to retrieve from Pinecone.")
    parser.add_argument(
        "--source-type",
        choices=["Email", "Note", "Literature", "Unknown"],
        default=None,
        help="Optional source_type filter.",
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
        help="Optional lipid signature filter (e.g. 'ALS-CSF-Core-6Ceramides'). Matches 'lipid_signatures' metadata.",
    )    

    args = parser.parse_args()

    result = query_rag(
        user_query=args.query,
        top_k=args.top_k,
        source_type=args.source_type,
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