from __future__ import annotations

from typing import Any, Dict, List, Type
from uuid import UUID

import pytest

from amprenta_rag.database.models import Feature, GeneBurden, Variant, VariantAnnotation
from amprenta_rag.variants.gene_burden import compute_gene_burden


class _FakeQuery:
    def __init__(self, model: Type, db: "_FakeDB", rows: List[Any]):
        self.model = model
        self._db = db
        self._rows = rows

    def filter(self, *args, **kwargs):  # noqa: ANN001, ANN002, ANN003
        # This is a lightweight in-memory stub; ignore SQLAlchemy filter expressions.
        return self

    def all(self) -> List[Any]:
        return list(self._rows)

    def delete(self, synchronize_session: bool = False):  # noqa: FBT001, ARG002
        if self.model is GeneBurden:
            self._db._gene_burdens = []
            self._rows = []
        return 0


class _FakeDB:
    def __init__(self, variants: List[Variant], annotations: List[VariantAnnotation], features: List[Feature]):
        self._variants = variants
        self._annotations = annotations
        self._features = features
        self._gene_burdens: List[GeneBurden] = []

    def query(self, model: Type):  # noqa: ANN001
        if model is Variant:
            return _FakeQuery(model, self, self._variants)
        if model is VariantAnnotation:
            return _FakeQuery(model, self, self._annotations)
        if model is Feature:
            return _FakeQuery(model, self, self._features)
        if model is GeneBurden:
            return _FakeQuery(model, self, self._gene_burdens)
        raise AssertionError(f"Unexpected model: {model}")

    def add(self, obj: Any) -> None:
        if isinstance(obj, GeneBurden):
            self._gene_burdens.append(obj)

    def commit(self) -> None:
        return


def test_compute_gene_burden_counts_and_score():
    vs_id = UUID("00000000-0000-0000-0000-000000000001")

    v1 = Variant(variant_set_id=vs_id, chromosome="1", position=1, ref_allele="A", alt_allele="T", gene_symbol="TP53")
    v1.id = UUID("00000000-0000-0000-0000-000000000010")
    v2 = Variant(variant_set_id=vs_id, chromosome="1", position=2, ref_allele="A", alt_allele="G", gene_symbol="TP53")
    v2.id = UUID("00000000-0000-0000-0000-000000000011")
    v3 = Variant(variant_set_id=vs_id, chromosome="1", position=3, ref_allele="C", alt_allele="T", gene_symbol="BRCA1")
    v3.id = UUID("00000000-0000-0000-0000-000000000012")

    a1 = VariantAnnotation(variant_id=v1.id, annotation_source="clinvar", clinical_significance="Pathogenic")
    a2 = VariantAnnotation(variant_id=v2.id, annotation_source="clinvar", clinical_significance="Uncertain significance")
    a3 = VariantAnnotation(variant_id=v3.id, annotation_source="clinvar", clinical_significance="Benign")

    f_tp53 = Feature(name="TP53", feature_type="gene")
    f_tp53.id = UUID("00000000-0000-0000-0000-000000000100")
    f_brca1 = Feature(name="BRCA1", feature_type="gene")
    f_brca1.id = UUID("00000000-0000-0000-0000-000000000101")

    db = _FakeDB(variants=[v1, v2, v3], annotations=[a1, a2, a3], features=[f_tp53, f_brca1])
    out = compute_gene_burden(vs_id, db)

    assert len(out) == 2
    by_gene: Dict[str, GeneBurden] = {g.gene_symbol: g for g in out}

    tp53 = by_gene["TP53"]
    assert tp53.n_variants == 2
    assert tp53.n_pathogenic == 1
    assert tp53.n_vus == 1
    assert tp53.n_benign == 0
    assert tp53.burden_score == pytest.approx(11.0)  # 1*10 + 1
    assert tp53.feature_id == f_tp53.id

    brca1 = by_gene["BRCA1"]
    assert brca1.n_variants == 1
    assert brca1.n_pathogenic == 0
    assert brca1.n_vus == 0
    assert brca1.n_benign == 1
    assert brca1.burden_score == pytest.approx(0.0)
    assert brca1.feature_id == f_brca1.id


