#!/usr/bin/env python3
"""
Helper script to set up Salmon transcriptome index.

Downloads reference transcriptome (if needed) and builds Salmon index.
"""

import sys
import subprocess
from pathlib import Path
from typing import Optional

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.genomics.pipeline import check_salmon_installed
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def download_human_transcriptome(output_dir: Path = None) -> Optional[Path]:
    """
    Download human reference transcriptome from Ensembl.

    Downloads the Homo sapiens cDNA FASTA file.

    Args:
        output_dir: Directory to save transcriptome file

    Returns:
        Path to downloaded transcriptome FASTA file
    """
    if output_dir is None:
        output_dir = Path("./reference_data")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Ensembl GRCh38 release (latest stable)
    # cDNA FASTA file URL for human
    transcriptome_url = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"

    transcriptome_file = output_dir / "Homo_sapiens.GRCh38.cdna.all.fa.gz"

    if transcriptome_file.exists():
        logger.info("[INDEX-SETUP] Transcriptome file already exists: %s", transcriptome_file)
        return transcriptome_file

    logger.info("[INDEX-SETUP] Downloading human transcriptome from Ensembl...")
    logger.info("[INDEX-SETUP] This may take several minutes (file is ~100MB compressed)...")

    try:
        import requests
        from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT

        headers = {"User-Agent": REPOSITORY_USER_AGENT}
        response = requests.get(transcriptome_url, headers=headers, stream=True, timeout=300)
        response.raise_for_status()

        total_size = 0
        with open(transcriptome_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    total_size += len(chunk)

                    # Log progress every 10MB
                    if total_size % (10 * 1024 * 1024) == 0:
                        mb = total_size / (1024 * 1024)
                        logger.info("[INDEX-SETUP] Downloaded %.1f MB...", mb)

        logger.info(
            "[INDEX-SETUP] Successfully downloaded transcriptome: %.1f MB",
            total_size / (1024 * 1024),
        )

        return transcriptome_file

    except Exception as e:
        logger.error("[INDEX-SETUP] Error downloading transcriptome: %r", e)
        return None


def build_salmon_index(
    transcriptome_file: Path,
    index_dir: Path = None,
    kmer_size: int = 31,
) -> Optional[Path]:
    """
    Build Salmon transcriptome index.

    Following Subprocess Protocol (STRICT):
    - Uses subprocess.run() to execute Salmon binary
    - Does NOT import salmon as a Python package
    - Checks installation before running

    Args:
        transcriptome_file: Path to transcriptome FASTA file (can be gzipped)
        index_dir: Directory for Salmon index output
        kmer_size: K-mer size for indexing (default: 31)

    Returns:
        Path to index directory if successful, None otherwise
    """
    # Check if Salmon is installed (Subprocess Protocol)
    if not check_salmon_installed():
        return None

    if index_dir is None:
        index_dir = Path("./salmon_index")
    else:
        index_dir = Path(index_dir)

    # Check if index already exists (newer Salmon versions use info.json, older use hash.bin)
    if index_dir.exists():
        if (index_dir / "info.json").exists() or (index_dir / "hash.bin").exists():
            logger.info("[INDEX-SETUP] Salmon index already exists: %s", index_dir)
            return index_dir

    logger.info("[INDEX-SETUP] Building Salmon index...")
    logger.info("[INDEX-SETUP] This may take 5-10 minutes...")

    index_dir.mkdir(parents=True, exist_ok=True)

    # Check if transcriptome is gzipped
    transcriptome_input = str(transcriptome_file)

    # Build Salmon index command
    cmd = [
        "salmon",
        "index",
        "-t", transcriptome_input,
        "-i", str(index_dir),
        "-k", str(kmer_size),
        "--gencode",  # Indicates GENCODE/Ensembl format
    ]

    logger.info("[INDEX-SETUP] Running: %s", " ".join(cmd))

    try:
        subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )

        # Verify index was created (newer Salmon versions use info.json, older use hash.bin)
        if (index_dir / "info.json").exists() or (index_dir / "hash.bin").exists():
            logger.info("[INDEX-SETUP] ✅ Salmon index built successfully: %s", index_dir)
            return index_dir
        else:
            logger.error(
                "[INDEX-SETUP] Index directory exists but required files (info.json or hash.bin) not found"
            )
            return None

    except subprocess.CalledProcessError as e:
        logger.error(
            "[INDEX-SETUP] Error building index: %s\nStdout: %s\nStderr: %s",
            e,
            e.stdout,
            e.stderr,
        )
        return None
    except FileNotFoundError:
        logger.error(
            "[INDEX-SETUP] Error: 'salmon' command not found. "
            "Please ensure Salmon is installed and in PATH."
        )
        return None


def setup_human_salmon_index(
    reference_dir: Path = None,
    index_dir: Path = None,
) -> Optional[Path]:
    """
    Complete setup: Download transcriptome and build Salmon index.

    Args:
        reference_dir: Directory for reference data
        index_dir: Directory for Salmon index

    Returns:
        Path to index directory if successful
    """
    if reference_dir is None:
        reference_dir = Path("./reference_data")
    if index_dir is None:
        index_dir = Path("./salmon_index")

    logger.info("[INDEX-SETUP] Setting up human transcriptome index...")

    # Step 1: Download transcriptome
    transcriptome_file = download_human_transcriptome(reference_dir)
    if not transcriptome_file:
        logger.error("[INDEX-SETUP] Failed to download transcriptome")
        return None

    # Step 2: Build index
    index_path = build_salmon_index(transcriptome_file, index_dir)
    if not index_path:
        logger.error("[INDEX-SETUP] Failed to build index")
        return None

    logger.info("[INDEX-SETUP] ✅ Complete! Index ready at: %s", index_path)
    return index_path


def main():
    """Main setup function."""
    print("\n" + "="*60)
    print("SALMON TRANSCRIPTOME INDEX SETUP")
    print("="*60)

    print("\nThis will:")
    print("  1. Download human reference transcriptome from Ensembl")
    print("  2. Build Salmon index")
    print("\n⚠️  This may take 10-15 minutes and requires ~500MB disk space.")

    # Setup paths
    reference_dir = Path("./reference_data")
    index_dir = Path("./salmon_index")

    index_path = setup_human_salmon_index(
        reference_dir=reference_dir,
        index_dir=index_dir,
    )

    if index_path:
        print("\n✅ Index setup complete!")
        print(f"   Index location: {index_path}")
        print("\n📋 You can now use this index for quantification:")
        print(f"   index_path = Path('{index_path}')")
    else:
        print("\n❌ Index setup failed. Check logs for details.")


if __name__ == "__main__":
    main()

