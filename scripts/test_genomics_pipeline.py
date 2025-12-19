#!/usr/bin/env python3
"""
Test the full genomics pipeline end-to-end.

This script tests:
1. Search ENA for FASTQ files
2. Download a small FASTQ file (subset for testing)
3. Run Salmon quantification
4. Extract gene counts
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.genomics.pipeline import (
    get_ena_fastqs,
    download_fastq,
    quantify_with_salmon,
    extract_gene_counts_from_salmon,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def test_ena_search():
    """Test 1: Search ENA for FASTQ files."""
    print("\n" + "="*60)
    print("TEST 1: ENA Search")
    print("="*60)

    # Search for a small RNA-seq study
    keyword = "Homo sapiens"
    limit = 3

    logger.info("[PIPELINE-TEST] Searching ENA for: %s (limit: %d)", keyword, limit)

    try:
        runs = get_ena_fastqs(keyword, limit=limit)

        if not runs:
            print("❌ FAIL: No FASTQ files found")
            return None

        print(f"✅ PASS: Found {len(runs)} FASTQ run(s)")
        for i, run in enumerate(runs, 1):
            print(f"   {i}. Run: {run['Run']}, Sample: {run['Sample']}")

        return runs[0]  # Return first run for testing

    except Exception as e:
        logger.error("[PIPELINE-TEST] ENA search failed: %r", e)
        print(f"❌ FAIL: ENA search error: {e}")
        return None


def test_fastq_download(run_info):
    """Test 2: Download FASTQ file (subset for testing)."""
    print("\n" + "="*60)
    print("TEST 2: FASTQ Download (Subset)")
    print("="*60)

    if not run_info:
        print("⚠️  SKIP: No run info available")
        return None

    run_id = run_info['Run']

    # Create test directory
    test_dir = Path("./test_genomics")
    test_dir.mkdir(exist_ok=True)

    logger.info("[PIPELINE-TEST] Downloading FASTQ subset for %s", run_id)

    try:
        # Download subset (first 1000 lines = ~250 reads for testing)
        downloaded_path = download_fastq(
            run_info=run_info,
            output_dir=test_dir,
            confirm=True,  # Required to proceed with download
            subset=True,  # Download subset only
            max_lines=1000,  # First 1000 lines for testing
        )

        if downloaded_path and downloaded_path.exists():
            file_size_mb = downloaded_path.stat().st_size / (1024 * 1024)
            print(f"✅ PASS: Downloaded FASTQ subset: {downloaded_path}")
            print(f"   Size: {file_size_mb:.2f} MB")
            return downloaded_path
        else:
            print("❌ FAIL: FASTQ download failed")
            return None

    except Exception as e:
        logger.error("[PIPELINE-TEST] FASTQ download failed: %r", e)
        print(f"❌ FAIL: Download error: {e}")
        return None


def test_salmon_quantification(fastq_path):
    """Test 3: Run Salmon quantification."""
    print("\n" + "="*60)
    print("TEST 3: Salmon Quantification")
    print("="*60)

    if not fastq_path or not fastq_path.exists():
        print("⚠️  SKIP: No FASTQ file available")
        return None

    index_path = Path("./salmon_index")

    if not index_path.exists():
        print(f"❌ FAIL: Salmon index not found at {index_path}")
        print("   -> Please run: python scripts/setup_salmon_index.py")
        return None

    logger.info("[PIPELINE-TEST] Running Salmon quantification")

    try:
        quant_file = quantify_with_salmon(
            fastq_path=fastq_path,
            index_path=index_path,
            output_dir=Path("./test_genomics/quants"),
            library_type="A",  # Auto-detect
            validate_mappings=True,
        )

        if quant_file and quant_file.exists():
            print(f"✅ PASS: Quantification complete: {quant_file}")
            return quant_file
        else:
            print("❌ FAIL: Quantification failed or output file not found")
            return None

    except Exception as e:
        logger.error("[PIPELINE-TEST] Salmon quantification failed: %r", e)
        print(f"❌ FAIL: Quantification error: {e}")
        return None


def test_gene_count_extraction(quant_file):
    """Test 4: Extract gene counts from Salmon output."""
    print("\n" + "="*60)
    print("TEST 4: Gene Count Extraction")
    print("="*60)

    if not quant_file or not quant_file.exists():
        print("⚠️  SKIP: No quantification file available")
        return None

    logger.info("[PIPELINE-TEST] Extracting gene counts")

    try:
        gene_counts = extract_gene_counts_from_salmon(quant_file)

        if gene_counts:
            print(f"✅ PASS: Extracted {len(gene_counts)} gene/transcript counts")

            # Show top 10 by count
            sorted_counts = sorted(
                gene_counts.items(),
                key=lambda x: x[1],
                reverse=True,
            )[:10]

            print("\n   Top 10 genes/transcripts by TPM:")
            for i, (gene_id, tpm) in enumerate(sorted_counts, 1):
                print(f"   {i:2d}. {gene_id}: {tpm:.2f} TPM")

            return gene_counts
        else:
            print("❌ FAIL: No gene counts extracted")
            return None

    except Exception as e:
        logger.error("[PIPELINE-TEST] Gene count extraction failed: %r", e)
        print(f"❌ FAIL: Extraction error: {e}")
        return None


def main():
    """Run full pipeline test."""
    print("\n" + "╔" + "="*58 + "╗")
    print("║" + " "*15 + "GENOMICS PIPELINE TEST" + " "*20 + "║")
    print("╚" + "="*58 + "╝")

    results = []

    # Test 1: ENA Search
    run_info = test_ena_search()
    results.append(("ENA Search", run_info is not None))

    # Test 2: FASTQ Download
    fastq_path = test_fastq_download(run_info)
    results.append(("FASTQ Download", fastq_path is not None))

    # Test 3: Salmon Quantification
    quant_file = test_salmon_quantification(fastq_path)
    results.append(("Salmon Quantification", quant_file is not None))

    # Test 4: Gene Count Extraction
    gene_counts = test_gene_count_extraction(quant_file)
    results.append(("Gene Count Extraction", gene_counts is not None))

    # Summary
    print("\n" + "="*60)
    print("PIPELINE TEST SUMMARY")
    print("="*60)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {test_name}")

    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All pipeline tests passed! Genomics pipeline is working.")
        return 0
    elif passed > 0:
        print(f"\n⚠️  Partial success: {passed}/{total} tests passed")
        return 1
    else:
        print("\n❌ All tests failed. Please check logs for details.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

