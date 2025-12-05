#!/usr/bin/env python3
"""
Test pathway enrichment with synthetic data.

Tests the pathway enrichment pipeline with a small set of known genes
to verify the ID mapping and pathway analysis works correctly.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.analysis.id_mapping import map_gene_to_kegg
from amprenta_rag.analysis.pathway_analysis import (
    map_features_to_kegg_pathways,
    perform_pathway_enrichment,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main():
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT TEST - SYNTHETIC DATA")
    print("=" * 80 + "\n")
    
    # Test with a small set of known cancer-related genes
    test_genes = {
        "TP53",  # Tumor protein p53
        "BRCA1",  # Breast cancer 1
        "BRCA2",  # Breast cancer 2
        "EGFR",  # Epidermal growth factor receptor
        "MYC",  # MYC proto-oncogene
    }
    
    print(f"Testing with {len(test_genes)} genes: {', '.join(sorted(test_genes))}\n")
    
    # Step 1: Test ID mapping
    print("Step 1: Testing Gene → KEGG ID Mapping")
    print("-" * 80)
    mapped_genes = {}
    for gene in test_genes:
        kegg_id = map_gene_to_kegg(gene, organism="hsa")
        if kegg_id:
            mapped_genes[gene] = kegg_id
            print(f"  ✅ {gene} → {kegg_id}")
        else:
            print(f"  ❌ {gene} → Not found")
    
    if not mapped_genes:
        print("\n❌ No genes could be mapped to KEGG IDs. Cannot proceed with pathway test.")
        return
    
    print(f"\n✅ Successfully mapped {len(mapped_genes)}/{len(test_genes)} genes\n")
    
    # Step 2: Test pathway mapping
    print("Step 2: Testing Gene → Pathway Mapping")
    print("-" * 80)
    print("Mapping genes to KEGG pathways (this may take a minute)...\n")
    
    try:
        pathways = map_features_to_kegg_pathways(
            features=set(mapped_genes.keys()),
            feature_type="gene",
        )
        
        print(f"✅ Found {len(pathways)} pathways\n")
        
        if pathways:
            print("Top 5 pathways:")
            for i, (pathway_id, pathway) in enumerate(list(pathways.items())[:5], 1):
                print(f"  {i}. {pathway.name} ({pathway_id})")
                print(f"     Features: {len(pathway.features)}")
                print(f"     Matched genes: {', '.join(list(pathway.features)[:5])}")
                if len(pathway.features) > 5:
                    print(f"     ... and {len(pathway.features) - 5} more")
                print()
        else:
            print("⚠️  No pathways found. This could be due to:")
            print("   - API rate limiting")
            print("   - Network issues")
            print("   - KEGG API changes")
            return
        
    except Exception as e:
        print(f"❌ Error mapping to pathways: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Step 3: Test pathway enrichment
    print("Step 3: Testing Pathway Enrichment Analysis")
    print("-" * 80)
    print("Performing enrichment analysis...\n")
    
    try:
        enrichment_results = perform_pathway_enrichment(
            input_features=set(mapped_genes.keys()),
            input_feature_types={"gene"},
            pathway_sources=["KEGG"],
            p_value_threshold=0.1,  # More lenient for testing
        )
        
        print(f"✅ Found {len(enrichment_results)} significantly enriched pathways\n")
        
        if enrichment_results:
            print("Top 5 enriched pathways:")
            for i, result in enumerate(enrichment_results[:5], 1):
                pathway = result.pathway
                print(f"\n  {i}. {pathway.name} ({pathway.source})")
                print(f"     Pathway ID: {pathway.pathway_id}")
                print(f"     P-value: {result.p_value:.4f}")
                print(f"     Adjusted P-value: {result.adjusted_p_value:.4f}")
                print(f"     Enrichment Ratio: {result.enrichment_ratio:.2f}x")
                print(f"     Matched Features: {result.input_features}/{result.pathway_size}")
                if result.matched_features:
                    print(f"     Genes: {', '.join(result.matched_features)}")
        else:
            print("⚠️  No significantly enriched pathways found")
            print("   This could be normal if the test genes don't cluster in pathways")
    
    except Exception as e:
        print(f"❌ Error in enrichment analysis: {e}")
        import traceback
        traceback.print_exc()
        return
    
    print("\n" + "=" * 80)
    print("✅ PATHWAY ENRICHMENT TEST COMPLETE!")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()

