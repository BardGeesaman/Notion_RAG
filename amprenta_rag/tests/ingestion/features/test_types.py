from __future__ import annotations

from amprenta_rag.ingestion.features.types import FeatureType


def test_feature_type_members():
    assert FeatureType.GENE.value == "gene"
    assert FeatureType.METABOLITE.value == "metabolite"
    assert "protein" in [m.value for m in FeatureType]

