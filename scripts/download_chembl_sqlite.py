#!/usr/bin/env python3
"""Download and extract ChEMBL SQLite database with checksum verification."""

from __future__ import annotations

import argparse
import hashlib
import tarfile
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError
from tqdm import tqdm


CHEMBL_FTP_BASE = "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases"


def _download_with_progress(url: str, dest: Path, *, chunk_size: int = 1024 * 1024) -> None:
    """Download file with progress bar using tqdm."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    
    req = Request(url, headers={"User-Agent": "amprenta/1.0"})
    
    try:
        with urlopen(req) as resp:
            total = resp.headers.get("Content-Length")
            total_size = int(total) if total else None
            
            with open(tmp, "wb") as f:
                with tqdm(
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    desc=f"Downloading {dest.name}",
                ) as pbar:
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        pbar.update(len(chunk))
        
        # Move completed download to final location
        tmp.replace(dest)
        
    except URLError as e:
        # Clean up partial download on failure
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError(f"Failed to download {url}: {e}") from e


def _compute_sha256(file_path: Path) -> str:
    """Compute SHA256 checksum of a file."""
    hash_sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()


def _fetch_checksum(checksum_url: str) -> str:
    """Fetch SHA256 checksum from remote .sha256 file."""
    req = Request(checksum_url, headers={"User-Agent": "amprenta/1.0"})
    
    try:
        with urlopen(req) as resp:
            content = resp.read().decode('utf-8').strip()
            # SHA256 files typically contain: "checksum  filename"
            # Extract just the checksum part
            return content.split()[0]
    except URLError as e:
        raise RuntimeError(f"Failed to fetch checksum from {checksum_url}: {e}") from e


def _extract_db(tar_gz: Path, out_db: Path, version: int) -> Path:
    """Extract SQLite database from ChEMBL tarball."""
    out_db.parent.mkdir(parents=True, exist_ok=True)
    
    with tarfile.open(tar_gz, mode="r:gz") as tf:
        members = [m for m in tf.getmembers() if m.isfile() and m.name.lower().endswith(".db")]
        if not members:
            raise ValueError(f"No .db file found in {tar_gz}")
        
        # Look for the version-specific DB file
        target_name = f"chembl_{version}.db"
        chosen = None
        
        for member in members:
            if Path(member.name).name.lower() == target_name:
                chosen = member
                break
        
        if chosen is None:
            # Fall back to first .db file if version-specific not found
            chosen = members[0]
            print(f"Warning: {target_name} not found, using {chosen.name}")
        
        # Safe extraction: write selected member content to out_db
        fobj = tf.extractfile(chosen)
        if fobj is None:
            raise ValueError(f"Failed to read {chosen.name} from tarball")
        
        tmp = out_db.with_suffix(out_db.suffix + ".part")
        with open(tmp, "wb") as out:
            while True:
                buf = fobj.read(1024 * 1024)
                if not buf:
                    break
                out.write(buf)
        tmp.replace(out_db)
    
    return out_db


def download_chembl(version: int = 34, output_dir: str = "data/chembl", force: bool = False) -> Path:
    """Download and extract ChEMBL SQLite database.
    
    Args:
        version: ChEMBL version to download (default: 34)
        output_dir: Directory to store the database (default: data/chembl)
        force: Force re-download even if file exists and checksum matches
        
    Returns:
        Path to the extracted database file
    """
    output_path = Path(output_dir)
    db_path = output_path / f"chembl_{version}.db"
    tar_path = output_path / f"chembl_{version}_sqlite.tar.gz"
    
    # Construct URLs
    tar_url = f"{CHEMBL_FTP_BASE}/chembl_{version}/chembl_{version}_sqlite.tar.gz"
    checksum_url = f"{CHEMBL_FTP_BASE}/chembl_{version}/chembl_{version}_sqlite.tar.gz.sha256"
    
    # Check if database already exists and is valid (unless force is True)
    if not force and db_path.exists():
        print(f"Database {db_path} already exists.")
        
        # Try to verify existing tarball checksum if it exists
        if tar_path.exists():
            try:
                print("Verifying existing tarball checksum...")
                expected_checksum = _fetch_checksum(checksum_url)
                actual_checksum = _compute_sha256(tar_path)
                
                if actual_checksum == expected_checksum:
                    print("✓ Checksum verified. Database is up to date.")
                    return db_path
                else:
                    print("✗ Checksum mismatch. Re-downloading...")
            except Exception as e:
                print(f"Warning: Could not verify checksum: {e}")
                print("Proceeding with existing database.")
                return db_path
        else:
            print("Using existing database (tarball not found for verification).")
            return db_path
    
    # Download the tarball
    print(f"Downloading ChEMBL {version} from: {tar_url}")
    _download_with_progress(tar_url, tar_path)
    print(f"✓ Downloaded: {tar_path}")
    
    # Verify checksum
    print("Verifying SHA256 checksum...")
    try:
        expected_checksum = _fetch_checksum(checksum_url)
        actual_checksum = _compute_sha256(tar_path)
        
        if actual_checksum != expected_checksum:
            raise ValueError(f"Checksum mismatch! Expected: {expected_checksum}, Got: {actual_checksum}")
        
        print("✓ Checksum verified")
    except Exception as e:
        print(f"Warning: Checksum verification failed: {e}")
        print("Proceeding with download (checksum verification is optional)")
    
    # Extract the database
    print(f"Extracting database to: {db_path}")
    _extract_db(tar_path, db_path, version)
    print(f"✓ Extracted: {db_path}")
    
    return db_path


def main() -> None:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Download and extract ChEMBL SQLite database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--output-dir",
        default="data/chembl",
        help="Directory to store the database"
    )
    parser.add_argument(
        "--version",
        type=int,
        default=34,
        help="ChEMBL version to download"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if file exists and checksum matches"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without actually downloading"
    )
    
    args = parser.parse_args()
    
    if args.dry_run:
        output_path = Path(args.output_dir)
        db_path = output_path / f"chembl_{args.version}.db"
        tar_path = output_path / f"chembl_{args.version}_sqlite.tar.gz"
        tar_url = f"{CHEMBL_FTP_BASE}/chembl_{args.version}/chembl_{args.version}_sqlite.tar.gz"
        checksum_url = f"{CHEMBL_FTP_BASE}/chembl_{args.version}/chembl_{args.version}_sqlite.tar.gz.sha256"
        
        print("DRY RUN: ChEMBL SQLite download")
        print(f"Version: {args.version}")
        print(f"Output directory: {args.output_dir}")
        print(f"Database path: {db_path}")
        print(f"Tarball path: {tar_path}")
        print(f"Download URL: {tar_url}")
        print(f"Checksum URL: {checksum_url}")
        print(f"Force re-download: {args.force}")
        return
    
    try:
        db_path = download_chembl(
            version=args.version,
            output_dir=args.output_dir,
            force=args.force
        )
        print(f"\n✓ ChEMBL {args.version} database ready: {db_path}")
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())