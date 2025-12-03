#!/usr/bin/env python3
"""
List page IDs for cross-omics testing.

Lists available:
- Programs
- Signatures
- Datasets (Experimental Data Assets)
- Features (Genes, Proteins, Metabolites, Lipids)
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def list_programs():
    """
    List all Programs with page IDs.
    
    Note: If the Programs database has multiple data sources (connected/synced databases),
    direct queries won't work. This function will try direct query first, then fall back
    to finding programs via experiments.
    """
    cfg = get_config()
    
    programs_db_id = cfg.notion.programs_db_id if hasattr(cfg.notion, "programs_db_id") else None
    
    if not programs_db_id:
        print("⚠️  Programs database ID not configured (set NOTION_PROGRAMS_DB_ID in .env)")
        return []
    
    programs = []
    url = f"{cfg.notion.base_url}/databases/{programs_db_id}/query"
    
    try:
        has_more = True
        start_cursor = None
        
        while has_more:
            payload = {"page_size": 100}
            if start_cursor:
                payload["start_cursor"] = start_cursor
            
            import requests
            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            data = resp.json()
            for page in data.get("results", []):
                props = page.get("properties", {}) or {}
                # Try "Program" first (actual schema), fallback to "Name"
                name_prop = props.get("Program", {}).get("title", []) or []
                if not name_prop:
                    name_prop = props.get("Name", {}).get("title", []) or []
                name = name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                page_id = page.get("id", "")
                
                programs.append({"name": name, "page_id": page_id})
            
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")
        
        return programs
    
    except requests.exceptions.HTTPError as e:
        error_text = str(e)
        if "multiple_data_sources" in error_text.lower():
            logger.warning(
                "Programs database has multiple data sources - direct query not supported. "
                "Try finding programs via experiments or provide program page IDs directly."
            )
        else:
            logger.error("Error listing programs: %r", e)
        return []
    except Exception as e:
        logger.error("Error listing programs: %r", e)
        return []


def list_signatures():
    """List all Signatures with page IDs."""
    cfg = get_config()
    
    if not hasattr(cfg.notion, "signature_db_id") or not cfg.notion.signature_db_id:
        print("⚠️  Signatures database ID not configured")
        return []
    
    signatures = []
    url = f"{cfg.notion.base_url}/databases/{cfg.notion.signature_db_id}/query"
    
    try:
        has_more = True
        start_cursor = None
        
        while has_more:
            payload = {"page_size": 100}
            if start_cursor:
                payload["start_cursor"] = start_cursor
            
            import requests
            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            data = resp.json()
            for page in data.get("results", []):
                props = page.get("properties", {}) or {}
                name_prop = props.get("Name", {}).get("title", []) or []
                name = name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                
                short_id_prop = props.get("Short ID", {}).get("rich_text", []) or []
                short_id = short_id_prop[0].get("plain_text", "") if short_id_prop else None
                
                modalities_prop = props.get("Modalities", {}).get("multi_select", []) or []
                modalities = [m.get("name", "") for m in modalities_prop] if modalities_prop else []
                
                page_id = page.get("id", "")
                
                signatures.append({
                    "name": name,
                    "short_id": short_id,
                    "modalities": modalities,
                    "page_id": page_id,
                })
            
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")
        
        return signatures
    
    except Exception as e:
        logger.error("Error listing signatures: %r", e)
        return []


def list_experiments():
    """List all Experiments with page IDs."""
    cfg = get_config()
    
    experiments_db_id = cfg.notion.experiments_db_id if hasattr(cfg.notion, "experiments_db_id") else None
    
    if not experiments_db_id:
        print("⚠️  Experiments database ID not configured (set NOTION_EXPERIMENTS_DB_ID in .env)")
        return []
    
    experiments = []
    url = f"{cfg.notion.base_url}/databases/{experiments_db_id}/query"
    
    try:
        has_more = True
        start_cursor = None
        
        while has_more:
            payload = {"page_size": 100}
            if start_cursor:
                payload["start_cursor"] = start_cursor
            
            import requests
            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            data = resp.json()
            for page in data.get("results", []):
                props = page.get("properties", {}) or {}
                name_prop = props.get("Name", {}).get("title", []) or []
                name = name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                page_id = page.get("id", "")
                
                experiments.append({"name": name, "page_id": page_id})
            
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")
        
        return experiments
    
    except Exception as e:
        logger.error("Error listing experiments: %r", e)
        return []


def list_datasets():
    """List all Datasets (Experimental Data Assets) with page IDs."""
    cfg = get_config()
    
    # Use the exp_data_db_id property
    exp_data_db_id = None
    if hasattr(cfg.notion, "exp_data_db_id"):
        exp_data_db_id = cfg.notion.exp_data_db_id
    elif hasattr(cfg, "NOTION_EXP_DATA_DB_ID"):
        exp_data_db_id = getattr(cfg, "NOTION_EXP_DATA_DB_ID", None)
    
    if not exp_data_db_id:
        print("⚠️  Experimental Data Assets database ID not configured")
        return []
    
    datasets = []
    url = f"{cfg.notion.base_url}/databases/{exp_data_db_id}/query"
    
    try:
        has_more = True
        start_cursor = None
        
        while has_more:
            payload = {"page_size": 100}
            if start_cursor:
                payload["start_cursor"] = start_cursor
            
            import requests
            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            data = resp.json()
            for page in data.get("results", []):
                props = page.get("properties", {}) or {}
                name_prop = props.get("Name", {}).get("title", []) or []
                name = name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                
                omics_type_prop = props.get("Omics Type", {}).get("select")
                omics_type = omics_type_prop.get("name", "") if omics_type_prop else None
                
                data_origin_prop = props.get("Data Origin", {}).get("select")
                data_origin = data_origin_prop.get("name", "") if data_origin_prop else None
                
                page_id = page.get("id", "")
                
                datasets.append({
                    "name": name,
                    "omics_type": omics_type,
                    "data_origin": data_origin,
                    "page_id": page_id,
                })
            
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")
        
        return datasets
    
    except Exception as e:
        logger.error("Error listing datasets: %r", e)
        return []


def list_features(feature_type: str):
    """List features of a specific type."""
    cfg = get_config()
    
    db_map = {
        "gene": (cfg.notion.gene_features_db_id if hasattr(cfg.notion, "gene_features_db_id") else None, "Gene Features"),
        "protein": (cfg.notion.protein_features_db_id if hasattr(cfg.notion, "protein_features_db_id") else None, "Protein Features"),
        "metabolite": (cfg.notion.metabolite_features_db_id if hasattr(cfg.notion, "metabolite_features_db_id") else None, "Metabolite Features"),
        "lipid": (cfg.notion.lipid_species_db_id if hasattr(cfg.notion, "lipid_species_db_id") else None, "Lipid Species"),
    }
    
    db_id, db_name = db_map.get(feature_type.lower(), (None, None))
    
    if not db_id:
        print(f"⚠️  {db_name} database ID not configured")
        return []
    
    features = []
    url = f"{cfg.notion.base_url}/databases/{db_id}/query"
    
    try:
        has_more = True
        start_cursor = None
        
        while has_more:
            payload = {"page_size": 100}
            if start_cursor:
                payload["start_cursor"] = start_cursor
            
            import requests
            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            data = resp.json()
            for page in data.get("results", []):
                props = page.get("properties", {}) or {}
                name_prop = props.get("Name", {}).get("title", []) or []
                name = name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                page_id = page.get("id", "")
                
                features.append({"name": name, "page_id": page_id})
            
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")
        
        return features
    
    except Exception as e:
        logger.error("Error listing %s features: %r", feature_type, e)
        return []


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="List page IDs for cross-omics testing")
    parser.add_argument(
        "--type",
        choices=["programs", "experiments", "signatures", "datasets", "features", "all"],
        default="all",
        help="Type of pages to list",
    )
    parser.add_argument(
        "--feature-type",
        choices=["gene", "protein", "metabolite", "lipid"],
        help="Feature type (required if --type features)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=10,
        help="Maximum number of items to show per category (default: 10)",
    )
    
    args = parser.parse_args()
    
    if args.type == "features" and not args.feature_type:
        parser.error("--feature-type is required when --type is 'features'")
    
    print("=" * 80)
    print("🔍 CROSS-OMICS TEST ID FINDER")
    print("=" * 80)
    print()
    
    if args.type in ["programs", "all"]:
        print(f"📋 PROGRAMS (showing up to {args.limit}):")
        print("-" * 80)
        programs = list_programs()
        if programs:
            for i, prog in enumerate(programs[:args.limit], 1):
                print(f"{i}. {prog['name']}")
                print(f"   ID: {prog['page_id']}")
                print()
            if len(programs) > args.limit:
                print(f"   ... and {len(programs) - args.limit} more")
                print()
        else:
            print("   No programs found")
            print()
    
    if args.type in ["experiments", "all"]:
        print(f"🧪 EXPERIMENTS (showing up to {args.limit}):")
        print("-" * 80)
        experiments = list_experiments()
        if experiments:
            for i, exp in enumerate(experiments[:args.limit], 1):
                print(f"{i}. {exp['name']}")
                print(f"   ID: {exp['page_id']}")
                print()
            if len(experiments) > args.limit:
                print(f"   ... and {len(experiments) - args.limit} more")
                print()
        else:
            print("   No experiments found")
            print()
    
    if args.type in ["signatures", "all"]:
        print(f"🔖 SIGNATURES (showing up to {args.limit}):")
        print("-" * 80)
        signatures = list_signatures()
        if signatures:
            for i, sig in enumerate(signatures[:args.limit], 1):
                print(f"{i}. {sig['name']}")
                if sig['short_id']:
                    print(f"   Short ID: {sig['short_id']}")
                if sig['modalities']:
                    print(f"   Modalities: {', '.join(sig['modalities'])}")
                print(f"   ID: {sig['page_id']}")
                print()
            if len(signatures) > args.limit:
                print(f"   ... and {len(signatures) - args.limit} more")
                print()
        else:
            print("   No signatures found")
            print()
    
    if args.type in ["datasets", "all"]:
        print(f"📊 DATASETS (showing up to {args.limit}):")
        print("-" * 80)
        datasets = list_datasets()
        if datasets:
            for i, ds in enumerate(datasets[:args.limit], 1):
                print(f"{i}. {ds['name']}")
                if ds['omics_type']:
                    print(f"   Omics Type: {ds['omics_type']}")
                if ds['data_origin']:
                    print(f"   Data Origin: {ds['data_origin']}")
                print(f"   ID: {ds['page_id']}")
                print()
            if len(datasets) > args.limit:
                print(f"   ... and {len(datasets) - args.limit} more")
                print()
        else:
            print("   No datasets found")
            print()
    
    if args.type in ["features", "all"]:
        if args.feature_type or args.type == "all":
            feature_types = [args.feature_type] if args.feature_type else ["gene", "protein", "metabolite", "lipid"]
            
            for ft in feature_types:
                print(f"🧬 {ft.upper()} FEATURES (showing up to {args.limit}):")
                print("-" * 80)
                features = list_features(ft)
                if features:
                    for i, feat in enumerate(features[:args.limit], 1):
                        print(f"{i}. {feat['name']}")
                        print(f"   ID: {feat['page_id']}")
                        print()
                    if len(features) > args.limit:
                        print(f"   ... and {len(features) - args.limit} more")
                        print()
                else:
                    print(f"   No {ft} features found")
                    print()
    
    print("=" * 80)
    print()
    print("💡 Usage examples:")
    print()
    print("  # Test with a program:")
    print(f"  python test_cross_omics.py --program <program_id_from_above>")
    print()
    print("  # Test with a signature:")
    print(f"  python test_cross_omics.py --signature <signature_id_from_above>")
    print()
    print("  # Test with a dataset:")
    print(f"  python test_cross_omics.py --dataset <dataset_id_from_above>")
    print()
    print("  # Test with a feature:")
    print(f"  python test_cross_omics.py --feature 'metabolite:<feature_name_from_above>'")
    print()


if __name__ == "__main__":
    main()

