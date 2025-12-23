"""HPOA genes_to_phenotype parser.

Parses `genes_to_phenotype.txt` from HPOA:
`http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt`

The file is tab-delimited with header comments starting with '#'.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional


@dataclass(frozen=True)
class HPOGeneAssociationRecord:
    ncbi_gene_id: Optional[str]
    gene_symbol: str
    hpo_id: str
    hpo_name: Optional[str] = None
    frequency: Optional[str] = None
    disease_id: Optional[str] = None
    source: Optional[str] = None


def iter_genes_to_phenotype_lines(lines: Iterable[str]) -> Iterator[HPOGeneAssociationRecord]:
    reader = csv.reader((ln for ln in lines if ln and not ln.startswith("#")), delimiter="\t")
    for row in reader:
        if not row:
            continue
        # Most common HPOA format (7 columns):
        # 0 ncbi_gene_id, 1 gene_symbol, 2 hpo_id, 3 hpo_name, 4 frequency, 5 disease_id, 6 source
        ncbi_gene_id = row[0].strip() if len(row) > 0 and row[0].strip() else None
        gene_symbol = row[1].strip() if len(row) > 1 else ""
        hpo_id = row[2].strip() if len(row) > 2 else ""
        hpo_name = row[3].strip() if len(row) > 3 and row[3].strip() else None
        frequency = row[4].strip() if len(row) > 4 and row[4].strip() else None
        disease_id = row[5].strip() if len(row) > 5 and row[5].strip() else None
        source = row[6].strip() if len(row) > 6 and row[6].strip() else None

        if not gene_symbol or not hpo_id:
            continue

        yield HPOGeneAssociationRecord(
            ncbi_gene_id=ncbi_gene_id,
            gene_symbol=gene_symbol,
            hpo_id=hpo_id,
            hpo_name=hpo_name,
            frequency=frequency,
            disease_id=disease_id,
            source=source,
        )


def parse_genes_to_phenotype(path: str | Path) -> List[HPOGeneAssociationRecord]:
    p = Path(path)
    text = p.read_text(encoding="utf-8", errors="replace").splitlines()
    return list(iter_genes_to_phenotype_lines(text))


__all__ = ["HPOGeneAssociationRecord", "parse_genes_to_phenotype", "iter_genes_to_phenotype_lines"]


