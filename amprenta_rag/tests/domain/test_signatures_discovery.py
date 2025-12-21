from uuid import uuid4

from amprenta_rag.domain.signatures_discovery import (
    Component,
    DiscoveredSignature,
    DiscoveryDatasetSummary,
)


def test_discovery_dataset_summary_defaults():
    ds = DiscoveryDatasetSummary(
        dataset_id=uuid4(),
        omics_type="gene",
        features={"A", "B"},
    )
    assert ds.matrix is None
    assert ds.directions is None


def test_component_and_discovered_signature():
    comp = Component(feature="F1", weight=0.5, direction="up")
    sig = DiscoveredSignature(
        name="Sig",
        modality="gene",
        components=[comp],
        support=2,
        provenance={"source": "test"},
    )
    assert sig.components[0].feature == "F1"
    assert sig.support == 2

