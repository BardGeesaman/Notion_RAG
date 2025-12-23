from __future__ import annotations

from pathlib import Path

from amprenta_rag.phenotype.hpo_parser import parse_genes_to_phenotype
from amprenta_rag.phenotype.phenotype_mapper import PhenotypeMapper


SAMPLE = """# comment
1\tGENE1\tHP:0000001\tTest phenotype\t\tOMIM:1\tHPO
2\tGENE2\tHP:0000001\tTest phenotype\t\tOMIM:1\tHPO
3\tGENE3\tHP:0000002\tOther phenotype\t\tOMIM:2\tHPO
"""


def test_parse_genes_to_phenotype(tmp_path: Path):
    p = tmp_path / "genes_to_phenotype.txt"
    p.write_text(SAMPLE, encoding="utf-8")
    recs = parse_genes_to_phenotype(p)
    assert len(recs) == 3
    assert recs[0].gene_symbol == "GENE1"
    assert recs[0].hpo_id == "HP:0000001"


def test_get_genes_for_hpo(tmp_path: Path):
    p = tmp_path / "genes_to_phenotype.txt"
    p.write_text(SAMPLE, encoding="utf-8")
    recs = parse_genes_to_phenotype(p)

    mapper = PhenotypeMapper(recs)
    genes = mapper.get_genes_for_hpo("HP:0000001")
    assert genes == ["GENE1", "GENE2"]


def test_expand_query_extracts_hpo_ids(tmp_path: Path):
    p = tmp_path / "genes_to_phenotype.txt"
    p.write_text(SAMPLE, encoding="utf-8")
    recs = parse_genes_to_phenotype(p)

    mapper = PhenotypeMapper(recs)
    out = mapper.expand_query("Patient has HP:0000002 and HP:0000001")
    assert out["hpo_ids"] == ["HP:0000001", "HP:0000002"]
    assert set(out["genes"]) == {"GENE1", "GENE2", "GENE3"}


