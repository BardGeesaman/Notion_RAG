"""
Genomics pipeline for downloading FASTQ files and extracting features.

Following the Genomics Pipeline Protocol:
- Search and download FASTQ files from ENA
- Quantify using Salmon/Kallisto
- Extract gene count matrices

Note: This requires Salmon or Kallisto to be installed in the user's environment.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import requests

from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# ENA API endpoint (Master Protocol)
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/search"


def get_ena_fastqs(
    keyword: str,
    limit: int = 3,
    filters: Optional[Dict[str, str]] = None,
) -> List[Dict[str, str]]:
    """
    Step 1: Find FASTQ download links for a study/keyword.

    Following Pipeline Protocol:
    - Uses ENA repository (which handles query formatting correctly)
    - Requests critical fields: fastq_ftp, run_accession, sample_alias
    - Returns valid runs with HTTP download URLs

    Args:
        keyword: Search keyword (e.g., "Homo sapiens" or organism name)
        limit: Maximum number of results
        filters: Optional filters (organism, library_strategy, etc.)

    Returns:
        List of dictionaries with Run, Sample, and URL (HTTP-converted)
    """
    logger.info("[GENOMICS-PIPELINE] Searching ENA for FASTQ files: %s", keyword)

    # Use ENA repository for reliable searching
    from amprenta_rag.ingestion.repositories import ENARepository

    ena = ENARepository()

    # Search using repository (handles query formatting correctly)
    keywords_list = [keyword]
    if filters and "organism" in filters:
        keywords_list.append(filters["organism"])

    study_ids = ena.search_studies(
        keywords=keywords_list,
        filters=filters,
        max_results=limit,
    )

    if not study_ids:
        logger.warning("[GENOMICS-PIPELINE] No runs found for keyword: %s", keyword)
        return []

    # Fetch metadata and FASTQ links for each run
    valid_runs = []

    for run_id in study_ids[:limit]:
        try:
            # Fetch metadata to get FASTQ links
            metadata = ena.fetch_study_metadata(run_id)
            if not metadata or not metadata.raw_metadata:
                continue

            run_data = metadata.raw_metadata
            fastq_ftp = run_data.get("fastq_ftp", "")

            if not fastq_ftp:
                continue

            # Parse FTP links
            ftp_links = [link.strip() for link in fastq_ftp.split(";") if link.strip()]

            if not ftp_links:
                continue

            # Convert FTP to HTTP for requests library (Pipeline Protocol)
            # ENA FTP: ftp.sra.ebi.ac.uk/vol1/fastq/...
            # Convert to: http://ftp.sra.ebi.ac.uk/vol1/fastq/...
            ftp_link = ftp_links[0]
            http_url = f"http://{ftp_link}" if not ftp_link.startswith("http") else ftp_link

            sample_alias = run_data.get("sample_alias", "")

            valid_runs.append({
                "Run": run_id,
                "Sample": sample_alias if sample_alias else run_data.get("sample_accession", ""),
                "URL": http_url,
                "FTP": f"ftp://{ftp_link}",
                "Filenames": ftp_links,
            })

        except Exception as e:
            logger.warning("[GENOMICS-PIPELINE] Error processing run %s: %r", run_id, e)
            continue

    logger.info("[GENOMICS-PIPELINE] Found %d valid runs with FASTQ files", len(valid_runs))
    return valid_runs


def download_fastq(
    run_info: Dict[str, str],
    output_dir: Path = None,
    confirm: bool = False,
    subset: bool = False,
    max_lines: int = 1000,
) -> Optional[Path]:
    """
    Step 2: Download FASTQ file from ENA.

    Following Pipeline Protocol:
    - Converts FTP links to HTTP
    - Uses requests.get(stream=True) for large files
    - Warns about large file sizes
    - Optional: Download subset (first N lines) for testing

    Args:
        run_info: Dictionary with Run, Sample, URL keys
        output_dir: Directory to save FASTQ file
        confirm: Whether user confirmed download (required for large files)
        subset: If True, only download first N lines (for testing)
        max_lines: Number of lines to download if subset=True

    Returns:
        Path to downloaded FASTQ file, or None if download skipped/failed
    """
    if output_dir is None:
        output_dir = Path("./fastq_downloads")

    output_dir.mkdir(parents=True, exist_ok=True)

    run_id = run_info.get("Run", "")
    url = run_info.get("URL", "")

    if not url:
        logger.error("[GENOMICS-PIPELINE] No URL provided for run %s", run_id)
        return None

    # Determine output filename
    filename = url.split("/")[-1] if "/" in url else f"{run_id}.fastq.gz"
    fastq_path = output_dir / filename

    # Check if already exists
    if fastq_path.exists():
        logger.info("[GENOMICS-PIPELINE] File already exists: %s", fastq_path)
        return fastq_path

    # Pipeline Protocol: Warn about large files, require confirmation
    if not confirm:
        logger.warning(
            "[GENOMICS-PIPELINE] FASTQ files can be very large (GB to TB). "
            "Skipping download. Set confirm=True to proceed."
        )
        return None

    logger.info("[GENOMICS-PIPELINE] Downloading FASTQ for run %s from %s", run_id, url)

    try:
        headers = {"User-Agent": REPOSITORY_USER_AGENT}

        if subset:
            # Download subset for testing (Pipeline Protocol)
            # For gzipped FASTQ files, we need to download, decompress, subset, then optionally re-compress
            logger.info(
                "[GENOMICS-PIPELINE] Downloading subset (first %d complete FASTQ records) for testing",
                max_lines // 4,  # FASTQ records are 4 lines each
            )

            import gzip

            resp = requests.get(url, headers=headers, stream=True, timeout=300)
            resp.raise_for_status()

            # Download compressed data in chunks and decompress on-the-fly
            lines_downloaded = 0
            records_downloaded = 0
            target_records = max_lines // 4  # Number of complete FASTQ records

            # Use uncompressed output for subset (easier for testing)
            if fastq_path.suffix == ".gz":
                # Remove .gz extension, but keep base filename
                fastq_path_uncompressed = fastq_path.with_suffix("")
            else:
                fastq_path_uncompressed = fastq_path

            with gzip.open(resp.raw, "rt") as gz_stream:
                with open(fastq_path_uncompressed, "wt") as f:
                    for line in gz_stream:
                        if records_downloaded >= target_records:
                            break

                        f.write(line)
                        lines_downloaded += 1

                        # FASTQ records are 4 lines: header, sequence, +, quality
                        if lines_downloaded % 4 == 0:
                            records_downloaded += 1

            logger.info(
                "[GENOMICS-PIPELINE] Downloaded %d complete FASTQ records (%d lines) to %s",
                records_downloaded,
                lines_downloaded,
                fastq_path_uncompressed,
            )

            # Return uncompressed path for subset
            fastq_path = fastq_path_uncompressed
        else:
            # Full download (streaming)
            resp = requests.get(url, headers=headers, stream=True, timeout=300)
            resp.raise_for_status()

            total_size = 0
            with open(fastq_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        total_size += len(chunk)

                        # Log progress every 100MB
                        if total_size % (100 * 1024 * 1024) == 0:
                            mb = total_size / (1024 * 1024)
                            logger.info(
                                "[GENOMICS-PIPELINE] Downloaded %.1f MB...",
                                mb,
                            )

            logger.info(
                "[GENOMICS-PIPELINE] Successfully downloaded %.1f MB to %s",
                total_size / (1024 * 1024),
                fastq_path,
            )

        return fastq_path

    except Exception as e:
        logger.error("[GENOMICS-PIPELINE] Error downloading FASTQ: %r", e)
        # Clean up partial download
        if fastq_path.exists():
            fastq_path.unlink()
        return None


def check_salmon_installed() -> bool:
    """
    Verifies that the Salmon binary is accessible.

    Following Subprocess Protocol (STRICT):
    - Checks if 'salmon' command exists in PATH
    - Verifies it's the bioinformatics tool (not mail server)
    - Does NOT import salmon as a Python package
    - NEVER use: import salmon (this is a mail server utility on PyPI!)

    Returns:
        True if Salmon is installed and accessible, False otherwise
    """
    try:
        result = subprocess.run(
            ["salmon", "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )

        version_output = result.stdout.strip()
        logger.info("[GENOMICS-PIPELINE] ✅ Found Salmon: %s", version_output)
        return True

    except FileNotFoundError:
        logger.error(
            "[GENOMICS-PIPELINE] ❌ Error: 'salmon' command not found. "
            "Please install it via Conda: 'conda install -c bioconda salmon'"
        )
        return False
    except Exception as e:
        logger.error("[GENOMICS-PIPELINE] Error checking Salmon installation: %r", e)
        return False


def quantify_with_salmon(
    fastq_path: Path,
    index_path: Path,
    output_dir: Path = None,
    library_type: str = "A",
    validate_mappings: bool = True,
) -> Optional[Path]:
    """
    Step 3: Run Salmon quantification on FASTQ file.

    Following Subprocess Protocol (STRICT):
    - Uses subprocess.run() to execute Salmon binary
    - Does NOT import salmon as a Python package (it's a mail server utility!)
    - Checks installation before running
    - Uses proper command-line flags (--validateMappings, --quiet)
    - Command arguments as list of strings to avoid shell injection

    Args:
        fastq_path: Path to input FASTQ file
        index_path: Path to pre-built Salmon transcriptome index
        output_dir: Directory for quantification output
        library_type: Library type ('A' = auto-detect, 'U' = unstranded, etc.)
        validate_mappings: Use --validateMappings flag (recommended for accuracy)

    Returns:
        Path to quant.sf output file, or None if quantification failed
    """
    logger.info(
        "[GENOMICS-PIPELINE] Running Salmon quantification for %s",
        fastq_path,
    )

    # Check if Salmon is installed (Subprocess Protocol)
    if not check_salmon_installed():
        return None

    if not fastq_path.exists():
        logger.error("[GENOMICS-PIPELINE] FASTQ file not found: %s", fastq_path)
        return None

    if not index_path.exists():
        logger.error("[GENOMICS-PIPELINE] Salmon index not found: %s", index_path)
        return None

    if output_dir is None:
        output_dir = Path("./quants") / fastq_path.stem
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build Salmon command following Subprocess Protocol
    # Command arguments as list of strings to avoid shell injection
    cmd = [
        "salmon",
        "quant",
        "-i", str(index_path),  # Path to pre-built index
        "-l", library_type,  # Library Type: Automatic (A)
        "-r", str(fastq_path),  # Input FASTQ (single-end reads)
        # For paired-end, use: "-1", read1_path, "-2", read2_path
        "-o", str(output_dir),  # Output Directory
    ]

    # Add recommended flags
    if validate_mappings:
        cmd.append("--validateMappings")  # Recommended flag for accuracy

    cmd.append("--quiet")  # Reduce console spam

    logger.info("[GENOMICS-PIPELINE] Running command: %s", " ".join(cmd))

    try:
        # Run Salmon via subprocess (Subprocess Protocol - STRICT)
        subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )

        quant_file = output_dir / "quant.sf"

        if quant_file.exists():
            logger.info(
                "[GENOMICS-PIPELINE] ✅ Quantification complete! Results saved to %s",
                quant_file,
            )
            return quant_file
        else:
            logger.error(
                "[GENOMICS-PIPELINE] Quantification completed but quant.sf not found in %s",
                output_dir,
            )
            return None

    except subprocess.CalledProcessError as e:
        logger.error(
            "[GENOMICS-PIPELINE] ❌ Salmon failed with error code %d",
            e.returncode,
        )
        if e.stdout:
            logger.debug("[GENOMICS-PIPELINE] Stdout: %s", e.stdout)
        if e.stderr:
            logger.error("[GENOMICS-PIPELINE] Stderr: %s", e.stderr)
        return None
    except Exception as e:
        logger.error("[GENOMICS-PIPELINE] Unexpected error running Salmon: %r", e)
        return None


def extract_gene_counts_from_salmon(quant_file: Path) -> Dict[str, float]:
    """
    Extract gene counts from Salmon quantification output.

    Reads quant.sf file and extracts transcript/gene counts.

    Args:
        quant_file: Path to quant.sf file from Salmon

    Returns:
        Dictionary mapping transcript/gene IDs to TPM (or NumReads) values
    """
    import pandas as pd

    logger.info("[GENOMICS-PIPELINE] Extracting gene counts from %s", quant_file)

    if not quant_file.exists():
        logger.error("[GENOMICS-PIPELINE] Quantification file not found: %s", quant_file)
        return {}

    try:
        # Read Salmon quant.sf file (tab-separated)
        df = pd.read_csv(quant_file, sep="\t")

        # quant.sf columns: Name, Length, EffectiveLength, TPM, NumReads
        # Use TPM (Transcripts Per Million) as the count value
        gene_counts = {}

        for _, row in df.iterrows():
            transcript_id = row.get("Name", "")
            tpm = row.get("TPM", 0.0)

            if transcript_id and tpm > 0:
                gene_counts[transcript_id] = float(tpm)

        logger.info(
            "[GENOMICS-PIPELINE] Extracted %d gene/transcript counts",
            len(gene_counts),
        )

        return gene_counts

    except Exception as e:
        logger.error(
            "[GENOMICS-PIPELINE] Error extracting gene counts: %r",
            e,
        )
        return {}


def quantify_with_kallisto(
    fastq_path: Path,
    index_path: Path,
    output_dir: Path = None,
) -> Optional[Path]:
    """
    Alternative quantification using Kallisto.

    Following Pipeline Protocol:
    - Uses Kallisto pseudo-alignment tool
    - Calls via subprocess (assumes Kallisto is installed)
    - Generates abundance.tsv output file

    Args:
        fastq_path: Path to input FASTQ file
        index_path: Path to pre-built Kallisto transcriptome index
        output_dir: Directory for quantification output

    Returns:
        Path to abundance.tsv output file, or None if quantification failed
    """
    logger.info(
        "[GENOMICS-PIPELINE] Running Kallisto quantification for %s",
        fastq_path,
    )

    if not fastq_path.exists():
        logger.error("[GENOMICS-PIPELINE] FASTQ file not found: %s", fastq_path)
        return None

    if not index_path.exists():
        logger.error("[GENOMICS-PIPELINE] Kallisto index not found: %s", index_path)
        return None

    if output_dir is None:
        output_dir = Path("./kallisto_quants") / fastq_path.stem
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build Kallisto command
    cmd = [
        "kallisto",
        "quant",
        "-i", str(index_path),  # Index file
        "-o", str(output_dir),  # Output directory
        "--single",  # Single-end reads
        "-l", "200",  # Estimated average fragment length
        "-s", "20",  # Estimated standard deviation
        str(fastq_path),  # Input FASTQ file
    ]

    logger.info("[GENOMICS-PIPELINE] Running command: %s", " ".join(cmd))

    try:
        subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )

        quant_file = output_dir / "abundance.tsv"

        if quant_file.exists():
            logger.info(
                "[GENOMICS-PIPELINE] Kallisto quantification complete! Results in %s",
                quant_file,
            )
            return quant_file
        else:
            logger.error(
                "[GENOMICS-PIPELINE] Quantification completed but abundance.tsv not found in %s",
                output_dir,
            )
            return None

    except FileNotFoundError:
        logger.error(
            "[GENOMICS-PIPELINE] Error: 'kallisto' tool not found. "
            "Please install Kallisto via conda/brew: conda install -c bioconda kallisto"
        )
        return None
    except subprocess.CalledProcessError as e:
        logger.error(
            "[GENOMICS-PIPELINE] Kallisto quantification failed: %s\nStdout: %s\nStderr: %s",
            e,
            e.stdout,
            e.stderr,
        )
        return None

