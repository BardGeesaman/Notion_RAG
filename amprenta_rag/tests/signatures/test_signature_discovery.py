from __future__ import annotations

from uuid import uuid4
from amprenta_rag.domain.signatures_discovery import DiscoveryDatasetSummary
from amprenta_rag.signatures import signature_discovery as sd


def test_discover_signatures_no_datasets():
    res = sd.discover_signatures_from_datasets([])
    assert res == []


def test_discover_signatures_single_cluster():
    d1 = uuid4()
    d2 = uuid4()
    
    ds1 = DiscoveryDatasetSummary(
        dataset_id=d1,
        omics_type="gene",
        features={"A", "B"},
        disease="ALS"
    )
    ds2 = DiscoveryDatasetSummary(
        dataset_id=d2,
        omics_type="gene",
        features={"A", "B", "C"},
        disease="ALS"
    )
    
    # A and B appear in 2 datasets (support=2)
    # C appears in 1 (support=1) -> filtered out by default min_support=2
    
    res = sd.discover_signatures_from_datasets([ds1, ds2])
    
    assert len(res) == 1
    sig = res[0]
    feats = sorted([c.feature for c in sig.components])
    assert feats == ["A", "B"]
    assert sig.support == 2
    assert sig.modality == "gene"
    assert str(d1) in sig.provenance["dataset_ids"]
    assert str(d2) in sig.provenance["dataset_ids"]


def test_discover_signatures_min_support_filtering():
    d1 = uuid4()
    ds1 = DiscoveryDatasetSummary(
        dataset_id=d1, 
        omics_type="gene", 
        features={"A"}, 
        disease="d"
    )
    
    # A only in 1 dataset, min_support=2 -> should be empty
    res = sd.discover_signatures_from_datasets([ds1], min_support=2)
    assert res == []
    
    # min_support=1 -> should find it
    res = sd.discover_signatures_from_datasets([ds1], min_support=1)
    assert len(res) == 1
    assert res[0].components[0].feature == "A"


def test_discover_signatures_overlap_clustering():
    d1, d2, d3 = uuid4(), uuid4(), uuid4()
    
    # Cluster 1: {A, B} (co-occur in d1, d2)
    # Cluster 2: {X, Y} (co-occur in d2, d3)
    
    ds1 = DiscoveryDatasetSummary(
        dataset_id=d1, 
        omics_type="gene", 
        features={"A", "B"}, 
        disease="d"
    )
    ds2 = DiscoveryDatasetSummary(
        dataset_id=d2, 
        omics_type="gene", 
        features={"A", "B", "X", "Y"}, 
        disease="d"
    )
    ds3 = DiscoveryDatasetSummary(
        dataset_id=d3, 
        omics_type="gene", 
        features={"X", "Y"}, 
        disease="d"
    )
    
    res = sd.discover_signatures_from_datasets([ds1, ds2, ds3], min_support=2, min_overlap=0.9)
    
    # Expect 2 distinct signatures
    assert len(res) == 2
    
    sigs = sorted([sorted([c.feature for c in s.components]) for s in res])
    assert sigs == [["A", "B"], ["X", "Y"]]
