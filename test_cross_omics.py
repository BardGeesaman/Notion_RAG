#!/usr/bin/env python3
"""
Quick test script for cross-omics RAG reasoning functions.

Usage:
    python test_cross_omics.py --feature "lipid:Cer(d18:1/16:0)"
    python test_cross_omics.py --help
"""

import argparse
import sys

from amprenta_rag.query.cross_omics_reasoning import (
    cross_omics_dataset_summary,
    cross_omics_feature_summary,
    cross_omics_program_summary,
    cross_omics_signature_summary,
)


def test_feature_summary(feature_type: str, feature_name: str):
    """Test cross-omics feature summary."""
    print(f"\n🔬 Testing Cross-Omics Feature Summary")
    print(f"   Feature: {feature_type}:{feature_name}\n")
    
    try:
        summary = cross_omics_feature_summary(
            feature_name=feature_name,
            feature_type=feature_type,
            top_k_datasets=5,
            top_k_chunks=20,
        )
        
        print("=" * 80)
        print(summary)
        print("=" * 80)
        print("\n✅ Test completed successfully!")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def test_program_summary(program_page_id: str):
    """Test cross-omics program summary."""
    print(f"\n🔬 Testing Cross-Omics Program Summary")
    print(f"   Program ID: {program_page_id}\n")
    
    try:
        summary = cross_omics_program_summary(
            program_page_id=program_page_id,
            top_k_per_omics=5,
        )
        
        print("=" * 80)
        print(summary)
        print("=" * 80)
        print("\n✅ Test completed successfully!")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def test_signature_summary(signature_page_id: str):
    """Test cross-omics signature summary."""
    print(f"\n🔬 Testing Cross-Omics Signature Summary")
    print(f"   Signature ID: {signature_page_id}\n")
    
    try:
        summary = cross_omics_signature_summary(
            signature_page_id=signature_page_id,
            top_k_datasets=5,
            top_k_chunks=20,
        )
        
        print("=" * 80)
        print(summary)
        print("=" * 80)
        print("\n✅ Test completed successfully!")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def test_dataset_summary(dataset_page_id: str):
    """Test cross-omics dataset summary."""
    print(f"\n🔬 Testing Cross-Omics Dataset Summary")
    print(f"   Dataset ID: {dataset_page_id}\n")
    
    try:
        summary = cross_omics_dataset_summary(
            dataset_page_id=dataset_page_id,
            top_k_chunks=20,
        )
        
        print("=" * 80)
        print(summary)
        print("=" * 80)
        print("\n✅ Test completed successfully!")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Test cross-omics RAG reasoning functions")
    
    parser.add_argument(
        "--feature",
        help="Test feature summary in format 'type:name' (e.g., 'lipid:Cer(d18:1/16:0)', 'gene:TP53')",
    )
    parser.add_argument(
        "--program",
        help="Test program summary with Notion page ID",
    )
    parser.add_argument(
        "--signature",
        help="Test signature summary with Notion page ID",
    )
    parser.add_argument(
        "--dataset",
        help="Test dataset summary with Notion page ID",
    )
    
    args = parser.parse_args()
    
    if not any([args.feature, args.program, args.signature, args.dataset]):
        parser.print_help()
        print("\n💡 Examples:")
        print("   python test_cross_omics.py --feature 'lipid:Cer(d18:1/16:0)'")
        print("   python test_cross_omics.py --feature 'gene:TP53'")
        print("   python test_cross_omics.py --program <program_page_id>")
        print("   python test_cross_omics.py --signature <signature_page_id>")
        print("   python test_cross_omics.py --dataset <dataset_page_id>")
        sys.exit(1)
    
    if args.feature:
        if ":" not in args.feature:
            print("❌ Error: --feature must be in format 'type:name'")
            sys.exit(1)
        feature_type, feature_name = args.feature.split(":", 1)
        feature_type = feature_type.lower().strip()
        feature_name = feature_name.strip()
        
        if feature_type not in ["gene", "protein", "metabolite", "lipid"]:
            print(f"❌ Error: Invalid feature_type '{feature_type}'. Must be: gene, protein, metabolite, lipid")
            sys.exit(1)
        
        test_feature_summary(feature_type, feature_name)
    
    elif args.program:
        test_program_summary(args.program)
    
    elif args.signature:
        test_signature_summary(args.signature)
    
    elif args.dataset:
        test_dataset_summary(args.dataset)


if __name__ == "__main__":
    main()

