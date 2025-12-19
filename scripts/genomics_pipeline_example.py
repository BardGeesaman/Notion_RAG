#!/usr/bin/env python3
"""
Example script for Genomics Pipeline Protocol.

Demonstrates:
1. Searching for FASTQ files in ENA
2. Downloading FASTQ files (with HTTP conversion)
3. Quantifying with Salmon/Kallisto
4. Extracting gene counts

Following the Genomics Pipeline Protocol from Gemini.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.genomics.pipeline import (
    download_fastq,
    get_ena_fastqs,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main():
    """Example genomics pipeline workflow."""
    print("\n" + "="*60)
    print("GENOMICS PIPELINE PROTOCOL EXAMPLE")
    print("="*60)
    print("\nFollowing Pipeline Protocol:")
    print("  • Search ENA for FASTQ files")
    print("  • Download FASTQ (FTP → HTTP conversion)")
    print("  • Quantify with Salmon/Kallisto")
    print("  • Extract gene counts")

    # Step 1: Search for FASTQ files
    print("\n" + "="*60)
    print("STEP 1: Search ENA for FASTQ files")
    print("="*60)

    keyword = "breast cancer"
    runs = get_ena_fastqs(keyword, limit=1)

    if not runs:
        print("❌ No FASTQ files found")
        return

    print(f"✅ Found {len(runs)} run(s) with FASTQ files:")
    for i, run_info in enumerate(runs, 1):
        print(f"\n   Run {i}:")
        print(f"   - Run ID: {run_info['Run']}")
        print(f"   - Sample: {run_info.get('Sample', 'N/A')}")
        print(f"   - HTTP URL: {run_info['URL'][:80]}...")
        print(f"   - Files: {len(run_info.get('Filenames', []))} FASTQ file(s)")

    # Select first run
    run_info = runs[0]
    run_id = run_info["Run"]

    print(f"\n📋 Using run: {run_id}")

    # Step 2: Download FASTQ (with confirmation)
    print("\n" + "="*60)
    print("STEP 2: Download FASTQ File")
    print("="*60)

    print("\n⚠️  FASTQ files can be very large (GB to TB).")
    print("   For this example, we'll download a subset (first 1000 lines) for testing.")
    print("   Set confirm=True and subset=False for full download.\n")

    fastq_path = download_fastq(
        run_info=run_info,
        output_dir=Path("./fastq_downloads"),
        confirm=True,  # User confirmation required
        subset=True,  # Download subset for testing
        max_lines=1000,
    )

    if not fastq_path:
        print("⚠️  FASTQ download skipped or failed")
        print("   (This is expected if confirm=False or download failed)")
        return

    print(f"✅ FASTQ file downloaded: {fastq_path}")

    # Step 3: Quantification with Salmon
    print("\n" + "="*60)
    print("STEP 3: Quantification with Salmon")
    print("="*60)

    print("\n⚠️  This step requires:")
    print("   1. Salmon installed (conda install -c bioconda salmon)")
    print("   2. Pre-built transcriptome index")
    print("\n   Skipping quantification in this example.")
    print("   Uncomment the code below when Salmon is installed.\n")

    # Uncomment when ready:
    # index_path = Path("./human_index")  # Path to your Salmon index
    # quant_file = quantify_with_salmon(
    #     fastq_path=fastq_path,
    #     index_path=index_path,
    #     output_dir=Path("./quants") / run_id,
    # )
    #
    # if quant_file:
    #     print(f"✅ Quantification complete: {quant_file}")
    #
    #     # Step 4: Extract gene counts
    #     print("\n" + "="*60)
    #     print("STEP 4: Extract Gene Counts")
    #     print("="*60)
    #
    #     gene_counts = extract_gene_counts_from_salmon(quant_file)
    #
    #     if gene_counts:
    #         print(f"✅ Extracted {len(gene_counts)} gene/transcript counts")
    #         print(f"\n   Sample counts:")
    #         for gene_id, count in list(gene_counts.items())[:10]:
    #             print(f"   - {gene_id}: {count:.2f} TPM")

    print("\n" + "="*60)
    print("PIPELINE EXAMPLE COMPLETE")
    print("="*60)
    print("\n📋 Summary:")
    print(f"   • Found: {len(runs)} run(s) with FASTQ files")
    print(f"   • Downloaded: {fastq_path.name if fastq_path else 'None'}")
    print(f"   • Quantification: Requires Salmon/Kallisto installation")

    print("\n📝 Next Steps:")
    print("   1. Install Salmon: conda install -c bioconda salmon")
    print("   2. Build transcriptome index")
    print("   3. Run quantification pipeline")
    print("   4. Extract gene counts and link to Postgres")


if __name__ == "__main__":
    main()

