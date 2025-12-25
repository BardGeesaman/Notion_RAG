from __future__ import annotations

import argparse
import tarfile
from pathlib import Path
from typing import Optional


DEFAULT_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_33_sqlite.tar.gz"


def _require_tqdm():
    try:
        from tqdm import tqdm  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("tqdm is required for progress bars (pip install tqdm)") from e
    return tqdm


def _download(url: str, dest: Path, *, chunk_size: int = 1024 * 1024) -> None:
    import urllib.request

    tqdm = _require_tqdm()

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")

    req = urllib.request.Request(url, headers={"User-Agent": "amprenta/1.0"})
    with urllib.request.urlopen(req) as resp:  # noqa: S310
        total = resp.headers.get("Content-Length")
        total_i: Optional[int] = None
        try:
            total_i = int(total) if total else None
        except Exception:  # noqa: BLE001
            total_i = None

        mode = "ab" if tmp.exists() else "wb"
        downloaded = tmp.stat().st_size if tmp.exists() else 0

        with open(tmp, mode) as f:
            if downloaded and resp.headers.get("Accept-Ranges") == "bytes":
                # Best-effort resume is not implemented for urllib; warn via progress total.
                pass

            with tqdm(
                total=total_i,
                initial=downloaded,
                unit="B",
                unit_scale=True,
                desc="Downloading chembl sqlite tar.gz",
            ) as pbar:
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    pbar.update(len(chunk))

    tmp.replace(dest)


def _extract_db(tar_gz: Path, out_db: Path) -> Path:
    out_db.parent.mkdir(parents=True, exist_ok=True)

    with tarfile.open(tar_gz, mode="r:gz") as tf:
        members = [m for m in tf.getmembers() if m.isfile() and m.name.lower().endswith(".db")]
        if not members:
            raise ValueError("No .db file found in tarball")

        # Prefer chembl_33.db if present; else take first.
        chosen = None
        for m in members:
            if Path(m.name).name.lower() == "chembl_33.db":
                chosen = m
                break
        if chosen is None:
            chosen = members[0]

        # Safe extraction: write selected member content to out_db without extracting arbitrary paths.
        fobj = tf.extractfile(chosen)
        if fobj is None:
            raise ValueError("Failed to read .db member from tarball")

        tmp = out_db.with_suffix(out_db.suffix + ".part")
        with open(tmp, "wb") as out:
            while True:
                buf = fobj.read(1024 * 1024)
                if not buf:
                    break
                out.write(buf)
        tmp.replace(out_db)

    return out_db


def main() -> None:
    ap = argparse.ArgumentParser(description="Download ChEMBL SQLite and extract chembl_33.db")
    ap.add_argument("--url", default=DEFAULT_URL, help="Tarball URL (default: latest chembl_33 sqlite tar.gz)")
    ap.add_argument("--out-db", default="data/chembl/chembl_33.db", help="Output DB path")
    ap.add_argument(
        "--cache-tar",
        default="data/chembl/chembl_33_sqlite.tar.gz",
        help="Where to store downloaded tarball",
    )
    ap.add_argument("--dry-run", action="store_true", help="Print planned actions without downloading")
    args = ap.parse_args()

    out_db = Path(args.out_db)
    tar_path = Path(args.cache_tar)

    if args.dry_run:
        print("DRY RUN: download chembl sqlite")
        print(f"URL: {args.url}")
        print(f"Tarball: {tar_path}")
        print(f"DB out: {out_db}")
        return

    print(f"Downloading: {args.url}")
    _download(args.url, tar_path)
    print(f"Saved tarball: {tar_path}")

    print(f"Extracting .db to: {out_db}")
    db_path = _extract_db(tar_path, out_db)
    print(f"Extracted DB: {db_path}")


if __name__ == "__main__":
    main()



