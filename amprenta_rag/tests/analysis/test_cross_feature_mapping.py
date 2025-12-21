from __future__ import annotations

from uuid import uuid4

from amprenta_rag.analysis import cross_feature_mapping as cfm


class FakeResult:
    def __init__(self, **kwargs):
        self.uniprot_id = kwargs.get("uniprot_id")
        self.gene_symbol = kwargs.get("gene_symbol")
        self.pathway_id = kwargs.get("pathway_id")


class FakeQuery:
    def __init__(self, rows):
        self._rows = rows

    def fetchall(self):
        return self._rows


class FakeDB:
    def __init__(self, rows):
        self._rows = rows

    def execute(self, _stmt):
        return FakeQuery(self._rows)


def test_map_gene_to_proteins():
    db = FakeDB([FakeResult(uniprot_id="P1"), FakeResult(uniprot_id="P2")])
    assert cfm.map_gene_to_proteins("GENE", db) == ["P1", "P2"]


def test_map_protein_to_genes():
    db = FakeDB([FakeResult(gene_symbol="G1"), FakeResult(gene_symbol="G2")])
    assert cfm.map_protein_to_genes("P", db) == ["G1", "G2"]


def test_get_cross_omics_feature_neighbors():
    fid = uuid4()
    db = FakeDB([FakeResult(pathway_id="PWY1"), FakeResult(pathway_id="PWY2")])
    neighbors = cfm.get_cross_omics_feature_neighbors(fid, db)
    assert neighbors["pathways"] == ["PWY1", "PWY2"]

