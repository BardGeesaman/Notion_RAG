from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
from typing import Dict, Iterable, Optional, Set, Tuple

import requests

from amprenta_rag.database.models import LINCSGene
from amprenta_rag.database.session import db_session


NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"


def _download(url: str, dest: Path, *, timeout: int = 60) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)


def _default_cache_dir() -> Path:
    return Path(os.environ.get("LINCS_CACHE_DIR", "data/lincs"))


def _ensure_gene_info(cache_dir: Path) -> Path:
    gene_info = cache_dir / "Homo_sapiens.gene_info.gz"
    if not gene_info.exists():
        _download(NCBI_GENE_INFO_URL, gene_info)
    return gene_info


def _load_landmark_ids(landmarks_file: Path) -> Set[int]:
    out: Set[int] = set()
    if not landmarks_file.exists():
        return out
    for ln in landmarks_file.read_text(encoding="utf-8", errors="replace").splitlines():
        t = ln.strip()
        if not t or t.startswith("#"):
            continue
        try:
            out.add(int(t))
        except Exception:
            continue
    return out


def _ensure_landmarks(cache_dir: Path, landmarks_file: Path) -> Tuple[Set[int], str]:
    """Return (landmark_entrez_ids, source_label)."""
    if landmarks_file.exists():
        return _load_landmark_ids(landmarks_file), f"file:{landmarks_file}"

    url = os.environ.get("LINCS_LANDMARKS_URL")
    if url:
        # Download into cache dir (offline-friendly once cached)
        cache_path = cache_dir / "landmark_entrez_ids.txt"
        if not cache_path.exists():
            _download(url, cache_path)
        return _load_landmark_ids(cache_path), f"url:{url}"

    return set(), "none"


def _parse_gene_info_for_entrez_ids(gene_info_gz: Path, target_entrez: Set[int], *, limit: Optional[int] = None) -> Dict[int, Tuple[str, str]]:
    """Parse gene_info.gz and return {entrez_id: (symbol, description)} for target IDs."""
    out: Dict[int, Tuple[str, str]] = {}
    with gzip.open(gene_info_gz, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            # columns: tax_id, GeneID, Symbol, LocusTag, Synonyms, dbXrefs, chromosome, map_location, description, ...
            try:
                entrez = int(parts[1])
            except Exception:
                continue
            if entrez not in target_entrez:
                continue
            symbol = (parts[2] or "").strip()
            desc = (parts[8] or "").strip()
            out[entrez] = (symbol, desc)
            if limit is not None and len(out) >= int(limit):
                break
    return out


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Populate LINCSGene metadata from NCBI gene_info.gz and mark landmark genes.")
    ap.add_argument("--cache-dir", type=str, default=str(_default_cache_dir()), help="Cache directory (default: LINCS_CACHE_DIR or data/lincs).")
    ap.add_argument(
        "--landmarks-file",
        type=str,
        default="data/lincs/landmark_entrez_ids.txt",
        help="Path to landmark Entrez IDs file (one per line).",
    )
    ap.add_argument("--dry-run", action="store_true", help="Do not write changes to DB.")
    ap.add_argument("--limit", type=int, default=None, help="Limit number of LINCSGene rows to update (debug).")
    ap.add_argument("--batch-size", type=int, default=1000, help="Commit every N updates.")

    args = ap.parse_args(list(argv) if argv is not None else None)

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    landmarks_file = Path(args.landmarks_file)

    gene_info_gz = _ensure_gene_info(cache_dir)
    landmark_ids, landmark_src = _ensure_landmarks(cache_dir, landmarks_file)

    with db_session() as db:
        q = db.query(LINCSGene).order_by(LINCSGene.entrez_id.asc())
        if args.limit:
            q = q.limit(int(args.limit))
        genes = q.all()

        target_entrez = {int(g.entrez_id) for g in genes if g.entrez_id is not None}
        meta = _parse_gene_info_for_entrez_ids(gene_info_gz, target_entrez)

        updated = 0
        missing = 0
        landmark_count = 0

        for idx, g in enumerate(genes, start=1):
            entrez = int(g.entrez_id)
            sym, desc = meta.get(entrez, ("", ""))
            if not sym and not desc:
                missing += 1

            is_landmark = bool(entrez in landmark_ids) if landmark_ids else False
            if is_landmark:
                landmark_count += 1

            changed = False
            if sym and g.gene_symbol != sym:
                g.gene_symbol = sym
                changed = True
            if desc and g.gene_title != desc:
                g.gene_title = desc
                changed = True
            if g.is_landmark != is_landmark:
                g.is_landmark = is_landmark
                changed = True

            if changed:
                updated += 1
                db.add(g)

            if (idx % int(args.batch_size)) == 0:
                if args.dry_run:
                    db.rollback()
                else:
                    db.commit()

        # final flush
        if args.dry_run:
            db.rollback()
        else:
            db.commit()

    print("LINCSGene metadata backfill complete")
    print(f"- cache_dir: {cache_dir}")
    print(f"- gene_info: {gene_info_gz}")
    print(f"- landmarks_source: {landmark_src}")
    if landmark_src == "none":
        print("  WARNING: no landmark list found; is_landmark was set to False for all genes")
    print(f"- total_rows_seen: {len(genes)}")
    print(f"- updated_rows: {updated}")
    print(f"- missing_metadata_rows: {missing}")
    print(f"- landmark_rows: {landmark_count}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


