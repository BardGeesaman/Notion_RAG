#!/usr/bin/env python3
"""
CLI script for ingesting screening data.

Ingests HTS campaign metadata, hit lists, and biochemical results into PostgreSQL,
and optionally promotes compounds to Notion.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.screening_ingestion import (
    ingest_biochemical_results,
    ingest_hts_campaign_metadata,
    ingest_hts_hit_list,
    promote_compounds_to_notion,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ingest screening data (HTS campaigns, hit lists, biochemical results)"
    )
    
    # Initialize database
    parser.add_argument(
        "--init-db",
        action="store_true",
        help="Initialize chemistry database (creates tables if they don't exist)",
    )
    
    # Campaign metadata
    parser.add_argument(
        "--campaign-metadata-file",
        type=Path,
        help="Path to HTS campaign metadata file (CSV/TSV)",
    )
    parser.add_argument(
        "--campaign-id",
        help="Campaign ID (auto-generated from metadata if not provided)",
    )
    
    # Hit list
    parser.add_argument(
        "--hit-list-file",
        type=Path,
        help="Path to HTS hit list file (CSV/TSV with SMILES and values)",
    )
    parser.add_argument(
        "--smiles-column",
        default="SMILES",
        help="Name of SMILES column (default: SMILES)",
    )
    parser.add_argument(
        "--value-column",
        default="value",
        help="Name of value/activity column (default: value)",
    )
    parser.add_argument(
        "--hit-threshold",
        type=float,
        help="Threshold for hit classification",
    )
    
    # Biochemical results
    parser.add_argument(
        "--biochemical-results-file",
        type=Path,
        help="Path to biochemical results file (CSV/TSV with SMILES and IC50/EC50/Ki/Kd)",
    )
    parser.add_argument(
        "--assay-name",
        help="Assay name (default: filename stem)",
    )
    
    # Promotion
    parser.add_argument(
        "--promote-compounds",
        action="store_true",
        help="Promote compounds from campaign to Notion",
    )
    parser.add_argument(
        "--program-id",
        help="Program ID to link promoted compounds to",
    )
    parser.add_argument(
        "--min-hit-value",
        type=float,
        help="Minimum hit value for promotion",
    )
    parser.add_argument(
        "--max-compounds",
        type=int,
        help="Maximum number of compounds to promote",
    )
    
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print("SCREENING DATA INGESTION")
    print("=" * 80)
    
    try:
        # Ingest campaign metadata
        if args.campaign_metadata_file:
            print(f"\n📋 Ingesting campaign metadata from {args.campaign_metadata_file}")
            campaign = ingest_hts_campaign_metadata(
                args.campaign_metadata_file,
                campaign_id=args.campaign_id,
            )
            print(f"✅ Campaign ingested: {campaign.campaign_id} - {campaign.campaign_name}")
            campaign_id = campaign.campaign_id
        else:
            campaign_id = args.campaign_id
        
        # Ingest hit list
        if args.hit_list_file:
            if not campaign_id:
                print("\n❌ Error: --campaign-id is required when ingesting hit list", file=sys.stderr)
                sys.exit(1)
            
            print(f"\n🧪 Ingesting hit list from {args.hit_list_file}")
            count = ingest_hts_hit_list(
                args.hit_list_file,
                campaign_id=campaign_id,
                smiles_column=args.smiles_column,
                value_column=args.value_column,
                hit_threshold=args.hit_threshold,
            )
            print(f"✅ Ingested {count} compounds")
        
        # Ingest biochemical results
        if args.biochemical_results_file:
            print(f"\n🔬 Ingesting biochemical results from {args.biochemical_results_file}")
            count = ingest_biochemical_results(
                args.biochemical_results_file,
                campaign_id=campaign_id,
                assay_name=args.assay_name,
            )
            print(f"✅ Ingested {count} biochemical results")
        
        # Promote compounds
        if args.promote_compounds:
            if not campaign_id:
                print("\n❌ Error: --campaign-id is required when promoting compounds", file=sys.stderr)
                sys.exit(1)
            
            print(f"\n📤 Promoting compounds from campaign {campaign_id} to Notion...")
            promoted = promote_compounds_to_notion(
                campaign_id=campaign_id,
                program_id=args.program_id,
                min_hit_value=args.min_hit_value,
                max_compounds=args.max_compounds,
            )
            print(f"✅ Promoted {len(promoted)} compounds (Notion integration pending)")
        
        if not any([
            args.init_db,
            args.campaign_metadata_file,
            args.hit_list_file,
            args.biochemical_results_file,
            args.promote_compounds,
        ]):
            parser.print_help()
            print("\n⚠️  No actions specified. Use --help for usage information.")
        
        print(f"\n{'=' * 80}\n")
        
    except Exception as e:
        logger.error("[INGEST][SCREENING] Ingestion failed: %r", e)
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

