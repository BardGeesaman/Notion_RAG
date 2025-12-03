#!/usr/bin/env python3
"""
Discover lipidomics studies in Metabolomics Workbench.

This script scans all MW studies and identifies those that contain lipidomics data,
then can optionally harvest them into Notion.
"""

import sys
from pathlib import Path
from typing import List, Dict, Any

sys.path.insert(0, str(Path(__file__).parent.parent))

import requests
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

MW_BASE_URL = "https://www.metabolomicsworkbench.org/rest"
MW_STUDY_SUMMARY_URL = f"{MW_BASE_URL}/study/study_id/ST/summary"


def discover_lipidomics_studies() -> List[Dict[str, Any]]:
    """
    Discover all lipidomics studies in Metabolomics Workbench.
    
    Returns:
        List of study dictionaries with study_id, title, disease, etc.
    """
    logger.info("[DISCOVER] Fetching all MW study summaries...")
    
    try:
        resp = requests.get(MW_STUDY_SUMMARY_URL, params={"format": "json"}, timeout=120)
        resp.raise_for_status()
        data = resp.json()
        
        # Normalize data structure
        if isinstance(data, dict):
            studies = []
            for sid, summary in data.items():
                if isinstance(summary, dict):
                    summary.setdefault("study_id", sid)
                    studies.append(summary)
            data = studies
        
        logger.info("[DISCOVER] Found %d total studies in MW", len(data))
        
        # Search for lipidomics studies
        lipidomics_keywords = [
            "lipid",
            "ceramide",
            "sphingomyelin",
            "sphingolipid",
            "triacylglycerol",
            "TAG",
            "fatty acid",
            "phospholipid",
            "glycerolipid",
            "glycerophospholipid",
        ]
        
        lipidomics_studies = []
        
        for study in data:
            title = str(study.get("study_title", "") or "").lower()
            summary = str(study.get("study_summary", "") or "").lower()
            study_id = study.get("study_id", "")
            
            if not study_id:
                continue
            
            # Check if study contains lipidomics keywords
            for keyword in lipidomics_keywords:
                if keyword in title or keyword in summary:
                    lipidomics_studies.append({
                        "study_id": study_id,
                        "title": study.get("study_title", ""),
                        "summary": study.get("study_summary", "") or "",
                        "disease": study.get("disease", "") or "",
                        "organism": study.get("organism", "") or "",
                    })
                    break
        
        logger.info("[DISCOVER] Found %d lipidomics studies", len(lipidomics_studies))
        return lipidomics_studies
    
    except Exception as e:
        logger.error("[DISCOVER] Error discovering lipidomics studies: %r", e)
        raise


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Discover lipidomics studies in Metabolomics Workbench"
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=50,
        help="Maximum number of results to display (default: 50)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output file path for study IDs (one per line)",
    )
    
    args = parser.parse_args()
    
    try:
        lipidomics_studies = discover_lipidomics_studies()
        
        print(f"\n{'=' * 80}")
        print(f"Lipidomics Studies in Metabolomics Workbench")
        print(f"{'=' * 80}")
        print(f"\nTotal lipidomics studies found: {len(lipidomics_studies)}\n")
        
        # Display results
        display_count = min(args.max_results, len(lipidomics_studies))
        for i, study in enumerate(lipidomics_studies[:display_count], 1):
            print(f"{i:3d}. {study['study_id']}: {study['title'][:70]}")
            if study['disease']:
                print(f"     Disease: {study['disease']}")
        
        if len(lipidomics_studies) > display_count:
            print(f"\n... and {len(lipidomics_studies) - display_count} more")
        
        # Write study IDs to file if requested
        if args.output:
            study_ids = [s['study_id'] for s in lipidomics_studies]
            args.output.write_text('\n'.join(study_ids) + '\n')
            print(f"\n✅ Wrote {len(study_ids)} study IDs to {args.output}")
        
        print(f"\n{'=' * 80}\n")
        print("To harvest a study, run:")
        print(f"  python scripts/harvest_mw_studies.py --study-id <STUDY_ID>")
        print("\nTo bulk harvest all lipidomics studies, use:")
        print(f"  python scripts/harvest_mw_studies.py --study-id-file {args.output or 'lipidomics_studies.txt'}")
        
        return 0
    
    except Exception as e:
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        logger.exception("Error in main")
        return 1


if __name__ == "__main__":
    sys.exit(main())

