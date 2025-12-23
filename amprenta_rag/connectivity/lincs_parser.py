"""Parsers for LINCS Level 5 GCT/GCTX files."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Optional


@dataclass(frozen=True)
class LINCSSignatureData:
    sig_id: str
    pert_iname: Optional[str]
    pert_id: Optional[str]
    cell_id: Optional[str]
    gene_expression: Dict[int, float]  # {entrez_id: z-score}


def _load_sidecar_meta(base_path: Path) -> Dict[str, Dict[str, str]]:
    """Load optional metadata sidecar if present.

    Supported:
    - {path}.meta.json: {sig_id: {pert_iname, pert_id, cell_id}}
    """
    meta_json = base_path.with_suffix(base_path.suffix + ".meta.json")
    if meta_json.exists():
        try:
            raw = json.loads(meta_json.read_text(encoding="utf-8"))
            if isinstance(raw, dict):
                out: Dict[str, Dict[str, str]] = {}
                for k, v in raw.items():
                    if isinstance(k, str) and isinstance(v, dict):
                        out[k] = {str(kk): str(vv) for kk, vv in v.items() if vv is not None}
                return out
        except Exception:
            return {}
    return {}


def parse_gct(gct_path: str) -> Iterator[LINCSSignatureData]:
    """Parse a GCT file and yield per-signature gene expression dicts.

    Handles GCT v1.2-like format:
    - line1: #1.2 (ignored)
    - line2: <n_rows>\t<n_cols>
    - line3: Name\tDescription\t<sig_id_1>\t<sig_id_2>...
    - subsequent lines: <entrez>\t<desc>\t<val1>\t<val2>...

    For metadata fields (pert_iname/pert_id/cell_id), we check an optional sidecar
    JSON: {gct_path}.meta.json.
    """
    p = Path(gct_path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    meta = _load_sidecar_meta(p)

    with p.open("r", encoding="utf-8", errors="replace") as f:
        header1 = f.readline()
        if not header1:
            return
        dims = f.readline()
        if not dims:
            return
        header = f.readline()
        if not header:
            return

        cols = header.rstrip("\n").split("\t")
        if len(cols) < 3:
            return
        sig_ids = cols[2:]

        # Accumulate expression per signature
        expr: Dict[str, Dict[int, float]] = {sid: {} for sid in sig_ids}

        for ln in f:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                entrez = int(parts[0])
            except Exception:
                continue
            vals = parts[2:]
            for sid, v in zip(sig_ids, vals):
                try:
                    expr[sid][entrez] = float(v)
                except Exception:
                    continue

        for sid in sig_ids:
            m = meta.get(sid, {})
            yield LINCSSignatureData(
                sig_id=sid,
                pert_iname=m.get("pert_iname"),
                pert_id=m.get("pert_id"),
                cell_id=m.get("cell_id"),
                gene_expression=expr.get(sid, {}),
            )


def parse_gctx(gctx_path: str) -> Iterator[LINCSSignatureData]:
    """Parse a GCTX file via cmapPy (if installed) and yield per-signature data."""
    p = Path(gctx_path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    try:
        from cmapPy.pandasGEXpress.parse import parse  # type: ignore[import-not-found]
    except Exception as e:
        raise ImportError('parse_gctx requires "cmapPy" to be installed') from e

    gctoo = parse.parse(str(p))
    data_df = getattr(gctoo, "data_df", None)
    col_meta = getattr(gctoo, "col_metadata_df", None)

    if data_df is None:
        return iter(())

    meta_by_sig: Dict[str, Dict[str, str]] = {}
    if col_meta is not None:
        try:
            for sid, row in col_meta.iterrows():
                meta_by_sig[str(sid)] = {
                    "pert_iname": str(row.get("pert_iname")) if row.get("pert_iname") is not None else "",
                    "pert_id": str(row.get("pert_id")) if row.get("pert_id") is not None else "",
                    "cell_id": str(row.get("cell_id")) if row.get("cell_id") is not None else "",
                }
        except Exception:
            meta_by_sig = {}

    # data_df is gene x signature; index should be entrez (often strings)
    for sid in list(data_df.columns):
        series = data_df[sid]
        gene_expression: Dict[int, float] = {}
        for idx, v in series.items():
            try:
                entrez = int(idx)
                gene_expression[entrez] = float(v)
            except Exception:
                continue
        m = meta_by_sig.get(str(sid), {})
        yield LINCSSignatureData(
            sig_id=str(sid),
            pert_iname=(m.get("pert_iname") or None),
            pert_id=(m.get("pert_id") or None),
            cell_id=(m.get("cell_id") or None),
            gene_expression=gene_expression,
        )


__all__ = ["LINCSSignatureData", "parse_gct", "parse_gctx"]


