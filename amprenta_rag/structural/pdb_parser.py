"""PDB metadata parsing utilities.

Extracts:
- Sequence from SEQRES records (per chain)
- Resolution from REMARK   2
- Experimental method from EXPDTA
"""

from __future__ import annotations

import re
from typing import Dict, Optional


_RES_RE = re.compile(r"RESOLUTION\.\s+([0-9.]+)\s+ANGSTROMS", re.I)


def parse_pdb_metadata(pdb_bytes: bytes) -> Dict[str, object]:
    """Parse a PDB file and extract metadata.

    Returns:
      {
        "sequences": {chain_id: "ACDE..."},
        "resolution": float|None,
        "method": str|None,
        "chain_ids": [..],
      }
    """
    text = (pdb_bytes or b"").decode("utf-8", errors="replace")
    lines = text.splitlines()

    sequences_3: Dict[str, list[str]] = {}
    method: Optional[str] = None
    resolution: Optional[float] = None

    for ln in lines:
        if ln.startswith("EXPDTA") and method is None:
            m = ln.replace("EXPDTA", "", 1).strip()
            method = m or None
        if ln.startswith("REMARK   2") and resolution is None:
            mm = _RES_RE.search(ln)
            if mm:
                try:
                    resolution = float(mm.group(1))
                except Exception:
                    resolution = None
        if ln.startswith("SEQRES"):
            parts = ln.split()
            # SEQRES serNum chainID numRes residues...
            if len(parts) >= 5:
                chain = parts[2]
                residues = parts[4:]
                sequences_3.setdefault(chain, []).extend(residues)

    sequences: Dict[str, str] = {}
    if sequences_3:
        try:
            from Bio.SeqUtils import seq1  # type: ignore[import-not-found]
        except Exception:
            seq1 = None  # type: ignore[assignment]

        for chain, res3 in sequences_3.items():
            if seq1 is None:
                # Fallback: use X for all residues
                sequences[chain] = "X" * len(res3)
            else:
                # seq1 supports whitespace-delimited 3-letter codes
                try:
                    sequences[chain] = str(seq1(" ".join(res3), custom_map={"SEC": "U", "PYL": "O"}))
                except Exception:
                    sequences[chain] = "X" * len(res3)

    chain_ids = sorted(sequences.keys()) if sequences else sorted(sequences_3.keys())
    return {
        "sequences": sequences,
        "resolution": resolution,
        "method": method,
        "chain_ids": chain_ids,
    }


def to_fasta(seqs: Dict[str, str], header_prefix: str = "chain") -> str:
    """Convert chain->sequence mapping to a FASTA string."""
    out = []
    for chain, seq in (seqs or {}).items():
        out.append(f">{header_prefix}_{chain}")
        out.append(seq or "")
    return "\n".join(out).strip() + ("\n" if out else "")


__all__ = ["parse_pdb_metadata", "to_fasta"]


