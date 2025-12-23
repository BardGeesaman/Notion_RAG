"""Phenotype -> gene mapping utilities (HPO)."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Set

from amprenta_rag.database.models import PhenotypeGeneAssociation
from amprenta_rag.database.session import db_session

from .hpo_parser import HPOGeneAssociationRecord


HPO_ID_RE = re.compile(r"HP:\d{7}")


class PhenotypeMapper:
    """Lookup service for phenotype->gene mapping."""

    def __init__(self, records: Optional[Iterable[HPOGeneAssociationRecord]] = None):
        self._index: Dict[str, Set[str]] = defaultdict(set)
        if records is not None:
            for r in records:
                self._index[r.hpo_id].add(r.gene_symbol)

    @classmethod
    def from_db(cls) -> "PhenotypeMapper":
        with db_session() as db:
            rows: List[PhenotypeGeneAssociation] = db.query(PhenotypeGeneAssociation).all()
        recs = [
            HPOGeneAssociationRecord(
                ncbi_gene_id=r.ncbi_gene_id,
                gene_symbol=r.gene_symbol,
                hpo_id=r.hpo_id,
                hpo_name=r.hpo_name,
                frequency=r.frequency,
                disease_id=r.disease_id,
                source=r.source,
            )
            for r in rows
        ]
        return cls(recs)

    def get_genes_for_hpo(self, hpo_id: str) -> List[str]:
        return sorted(self._index.get(hpo_id, set()))

    def expand_query(self, query: str, limit: int = 200) -> Dict[str, Any]:
        """Extract HPO IDs from a free-form query and return union of mapped genes."""
        hpo_ids = sorted(set(HPO_ID_RE.findall(query or "")))
        genes: Set[str] = set()
        for hid in hpo_ids:
            genes.update(self._index.get(hid, set()))
        out = {
            "hpo_ids": hpo_ids,
            "genes": sorted(list(genes))[:limit],
            "gene_count": len(genes),
        }
        return out


_mapper_cache: Optional["PhenotypeMapper"] = None


def get_mapper(force_reload: bool = False) -> PhenotypeMapper:
    """Get a cached mapper built from the DB (singleton cache)."""
    global _mapper_cache
    if force_reload or _mapper_cache is None:
        _mapper_cache = PhenotypeMapper.from_db()
    return _mapper_cache


__all__ = ["PhenotypeMapper", "HPO_ID_RE", "get_mapper"]


