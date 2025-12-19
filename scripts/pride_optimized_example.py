#!/usr/bin/env python3
"""
Optimized PRIDE extraction example following best practices.

Demonstrates the recommended approach:
- Priority 1: *.mzTab files
- Priority 2: protein_groups.txt (MaxQuant)
- Priority 3: *.xls/*.xlsx spreadsheets
- Uses pandas for data extraction
"""

import io
import sys
import time
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import requests

from amprenta_rag.ingestion.repositories.pride import PRIDERepository
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

PRIDE_BASE_URL = "https://www.ebi.ac.uk/pride/ws/archive/v2"


def search_pride_projects(keyword: str, limit: int = 5) -> list[str]:
    """
    Step 1: Find PXD IDs using PRIDE API v2.

    Args:
        keyword: Search term (e.g., "breast cancer")
        limit: Maximum number of results

    Returns:
        List of PRIDE project IDs (e.g., ["PXD012345"])
    """
    logger.info("[PRIDE] Searching PRIDE for: %s (limit=%d)", keyword, limit)

    try:
        url = f"{PRIDE_BASE_URL}/search/projects"
        params = {"keyword": keyword, "pageSize": limit}

        time.sleep(1.0)  # Rate limiting
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()

        # PRIDE v2 returns a wrapper with _embedded structure
        data = resp.json()
        projects = data.get("_embedded", {}).get("projects", [])

        pxd_ids = [p.get("accession", "") for p in projects if p.get("accession")]
        pxd_ids = [pid for pid in pxd_ids if pid and pid.startswith("PXD")]

        logger.info("[PRIDE] Found %d projects: %s", len(pxd_ids), pxd_ids[:5])
        return pxd_ids

    except Exception as e:
        logger.error("[PRIDE] Search error: %r", e)
        return []


def find_best_result_file(pxd_id: str) -> tuple[str, str] | None:
    """
    Step 2 & 3: Find the best 'Result' file with priority ordering.

    Priority:
    1. *.mzTab (Community standard)
    2. protein_groups.txt (MaxQuant standard)
    3. *.xls or *.xlsx (Supplied spreadsheets)

    Args:
        pxd_id: PRIDE project ID (e.g., "PXD012345")

    Returns:
        Tuple of (filename, download_url) or None if not found
    """
    logger.info("[PRIDE] Scanning files for %s...", pxd_id)

    try:
        repo = PRIDERepository()

        # Get all files
        data_files = repo.fetch_study_data_files(study_id=pxd_id)

        if not data_files:
            logger.warning("[PRIDE] No files found for %s", pxd_id)
            return None

        logger.info("[PRIDE] Found %d files for %s", len(data_files), pxd_id)

        # Priority 1: Look for *.mzTab files
        mztab_files = [
            f for f in data_files
            if f.filename.lower().endswith(".mztab")
        ]

        if mztab_files:
            target_file = mztab_files[0]
            logger.info("[PRIDE] Priority 1: Found mzTab file: %s", target_file.filename)
            return target_file.filename, target_file.download_url

        # Priority 2: Look for protein_groups.txt (MaxQuant)
        maxquant_files = [
            f for f in data_files
            if "protein_groups.txt" in f.filename.lower()
        ]

        if maxquant_files:
            target_file = maxquant_files[0]
            logger.info("[PRIDE] Priority 2: Found MaxQuant file: %s", target_file.filename)
            return target_file.filename, target_file.download_url

        # Priority 3: Look for Excel files (*.xls, *.xlsx)
        excel_files = [
            f for f in data_files
            if f.filename.lower().endswith((".xls", ".xlsx"))
        ]

        if excel_files:
            target_file = excel_files[0]
            logger.info("[PRIDE] Priority 3: Found Excel file: %s", target_file.filename)
            return target_file.filename, target_file.download_url

        # Fallback: Look for any TSV/CSV result file
        result_files = [
            f for f in data_files
            if (f.filename.lower().endswith((".tsv", ".csv", ".txt")) and
                any(kw in f.filename.lower() for kw in ["protein", "result", "identification"]))
        ]

        if result_files:
            target_file = result_files[0]
            logger.info("[PRIDE] Fallback: Found result file: %s", target_file.filename)
            return target_file.filename, target_file.download_url

        logger.warning("[PRIDE] No easy-to-parse result file (mzTab/txt/xls) found for %s", pxd_id)
        return None

    except Exception as e:
        logger.error("[PRIDE] Error finding result file for %s: %r", pxd_id, e)
        return None


def extract_proteins_from_mztab(download_url: str) -> pd.DataFrame:
    """
    Extract protein data from mzTab file.

    mzTab format:
    - Header lines start with 'PRH' (Protein Header)
    - Data lines start with 'PRT' (Protein)

    Args:
        download_url: URL to download mzTab file

    Returns:
        DataFrame with protein data
    """
    logger.info("[PRIDE] Downloading and parsing mzTab file...")

    # Convert FTP to HTTP if needed
    if download_url.startswith("ftp://"):
        download_url = download_url.replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")

    time.sleep(1.0)  # Rate limiting
    response = requests.get(download_url, timeout=300, stream=True)
    response.raise_for_status()

    # Read file content
    content = response.text if hasattr(response, 'text') else response.content.decode('utf-8', errors='ignore')
    lines = content.split('\n')

    # Find PRH (Protein Header) line - contains column names
    header_line = None
    for line in lines:
        if line.startswith('PRH'):
            header_line = line
            break

    if not header_line:
        logger.warning("[PRIDE] No PRH header line found in mzTab file")
        return pd.DataFrame()

    # Extract column names (skip 'PRH\t' prefix)
    cols = header_line.replace('PRH\t', '').split('\t')

    # Extract PRT (Protein) data lines
    data_lines = [line.replace('PRT\t', '') for line in lines if line.startswith('PRT')]

    if not data_lines:
        logger.warning("[PRIDE] No PRT data lines found in mzTab file")
        return pd.DataFrame()

    # Load into pandas DataFrame
    df = pd.read_csv(
        io.StringIO('\n'.join(data_lines)),
        sep='\t',
        names=cols,
        low_memory=False,
    )

    logger.info("[PRIDE] Extracted %d proteins from mzTab file", len(df))
    return df


def extract_proteins_from_maxquant(download_url: str) -> pd.DataFrame:
    """
    Extract protein data from MaxQuant protein_groups.txt file.

    Args:
        download_url: URL to download protein_groups.txt

    Returns:
        DataFrame with protein data
    """
    logger.info("[PRIDE] Downloading and parsing MaxQuant protein_groups.txt...")

    # Convert FTP to HTTP if needed
    if download_url.startswith("ftp://"):
        download_url = download_url.replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")

    time.sleep(1.0)  # Rate limiting
    response = requests.get(download_url, timeout=300)
    response.raise_for_status()

    # Load directly into pandas (tab-separated)
    df = pd.read_csv(
        io.StringIO(response.text),
        sep='\t',
        low_memory=False,
    )

    logger.info("[PRIDE] Extracted %d protein groups from MaxQuant file", len(df))
    return df


def extract_proteins_from_excel(download_url: str) -> pd.DataFrame:
    """
    Extract protein data from Excel file.

    Args:
        download_url: URL to download Excel file

    Returns:
        DataFrame with protein data
    """
    logger.info("[PRIDE] Downloading and parsing Excel file...")

    # Convert FTP to HTTP if needed
    if download_url.startswith("ftp://"):
        download_url = download_url.replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")

    time.sleep(1.0)  # Rate limiting
    response = requests.get(download_url, timeout=300)
    response.raise_for_status()

    # Load Excel file into pandas
    df = pd.read_excel(
        io.BytesIO(response.content),
        engine='openpyxl' if download_url.endswith('.xlsx') else None,
    )

    logger.info("[PRIDE] Extracted %d rows from Excel file", len(df))
    return df


def get_project_data(pxd_id: str) -> pd.DataFrame | None:
    """
    Step 2 & 3: Find the best 'Result' file and extract protein data.

    Args:
        pxd_id: PRIDE project ID

    Returns:
        DataFrame with protein data or None
    """
    result = find_best_result_file(pxd_id)

    if not result:
        return None

    filename, download_url = result

    try:
        # Parse based on file type
        if filename.lower().endswith('.mztab'):
            df = extract_proteins_from_mztab(download_url)
        elif 'protein_groups.txt' in filename.lower():
            df = extract_proteins_from_maxquant(download_url)
        elif filename.lower().endswith(('.xls', '.xlsx')):
            df = extract_proteins_from_excel(download_url)
        else:
            logger.warning("[PRIDE] Unknown file type: %s", filename)
            return None

        return df

    except Exception as e:
        logger.error("[PRIDE] Error extracting data from %s: %r", filename, e)
        return None


def extract_protein_accessions(df: pd.DataFrame) -> set[str]:
    """
    Extract protein accessions from DataFrame.

    Looks for common column names:
    - accessions, accession
    - protein, protein id
    - uniprot, uniprot id

    Args:
        df: DataFrame with protein data

    Returns:
        Set of normalized protein accessions
    """
    protein_set: set[str] = set()

    if df.empty:
        return protein_set

    # Find protein accession column
    accession_col = None
    for col in df.columns:
        col_lower = col.lower().strip()
        if any(keyword in col_lower for keyword in [
            "accession", "accessions", "protein", "protein id",
            "uniprot", "uniprot id", "proteinid"
        ]):
            accession_col = col
            break

    if not accession_col:
        # Default to first column
        accession_col = df.columns[0]
        logger.warning("[PRIDE] No accession column found, using first column: %s", accession_col)

    # Extract and normalize protein accessions
    for protein_id in df[accession_col]:
        if pd.isna(protein_id):
            continue

        protein_id_str = str(protein_id).strip().strip('"').strip("'")

        if protein_id_str and not protein_id_str.startswith('#'):
            # Handle formats like "sp|P12345|PROTEIN_NAME" -> extract "P12345"
            if '|' in protein_id_str:
                protein_id_str = protein_id_str.split('|')[1] if len(protein_id_str.split('|')) > 1 else protein_id_str.split('|')[0]

            # Take first word/token
            protein_id_str = protein_id_str.split()[0] if protein_id_str.split() else protein_id_str

            if protein_id_str and len(protein_id_str) > 2:
                protein_set.add(protein_id_str)

    logger.info("[PRIDE] Extracted %d unique protein accessions", len(protein_set))
    return protein_set


def main():
    """Example usage."""
    # 1. Search for studies
    keyword = "breast cancer"
    studies = search_pride_projects(keyword, limit=3)

    print(f"\n✅ Found {len(studies)} studies: {studies[:5]}")

    # 2. Extract data from first study
    if studies:
        pxd_id = studies[0]
        print(f"\n📊 Extracting data from {pxd_id}...")

        df = get_project_data(pxd_id)

        if df is not None and not df.empty:
            print("\n--- Preview: Protein Data ---")
            print(df.head())

            # Extract protein accessions
            print(f"\n🧬 Extracting protein accessions...")
            proteins = extract_protein_accessions(df)
            print(f"✅ Extracted {len(proteins)} unique protein accessions")
            if proteins:
                print(f"Sample proteins: {list(proteins)[:10]}")


if __name__ == "__main__":
    main()

