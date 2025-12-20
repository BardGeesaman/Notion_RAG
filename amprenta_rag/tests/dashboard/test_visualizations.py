from __future__ import annotations

"""
Tests for visualization helper functions used by Streamlit dashboard modules.

We exercise:
- Volcano plot data extraction on a realistic dataset (volcano._get_volcano_data)
- Heatmap matrix building and clustering for multiple datasets
- PCA feature matrix construction for program-linked datasets
- Signature network construction
- Handling of empty inputs (no data)

These tests avoid hitting real Postgres or Streamlite by constructing in-memory
SQLAlchemy model instances and lightweight dummy Session/Query objects.
"""

from typing import Dict, Iterable, List
from uuid import uuid4


from amprenta_rag.database.models import Dataset, Feature, Signature
from scripts.dashboard.pages.visualizations.heatmap import _build_matrix as build_heatmap_matrix
from scripts.dashboard.pages.visualizations.heatmap import _cluster_matrix
from scripts.dashboard.pages.visualizations.network import _build_network
from scripts.dashboard.pages.visualizations.pca import _build_matrix as build_pca_matrix
from scripts.dashboard.pages.visualizations.volcano import _get_volcano_data


class DummyQuery:
    """Simple stand-in for SQLAlchemy Query."""

    def __init__(self, result):
        self._result = result

    def filter(self, *args, **kwargs):
        # For these tests, we ignore the filter criteria and just return self
        return self

    def all(self):
        return self._result

    def first(self):
        if isinstance(self._result, list):
            return self._result[0] if self._result else None
        return self._result


class DummySession:
    """Simple stand-in for SQLAlchemy Session."""

    def __init__(self, mapping: Dict[object, List[object]]):
        """
        mapping: dict from model class to list of instances to return from query().
        """
        self._mapping = mapping

    def query(self, model):
        return DummyQuery(self._mapping.get(model, []))


def make_feature(name: str, ftype: str = "gene", fc: float | None = None, p: float | None = None) -> Feature:
    ext: Dict[str, float] = {}
    if fc is not None:
        ext["log2FC"] = fc
    if p is not None:
        ext["pvalue"] = p
    return Feature(id=uuid4(), name=name, feature_type=ftype, external_ids=ext)  # type: ignore[arg-type]


def make_dataset(name: str, omics: str, disease: Iterable[str], features: list[Feature]) -> Dataset:
    ds = Dataset(
        id=uuid4(),  # type: ignore[arg-type]
        name=name,
        omics_type=omics,
        disease=list(disease),  # type: ignore[arg-type]
    )
    ds.features = list(features)  # type: ignore[assignment]
    return ds


def test_volcano_get_data_with_realistic_dataset():
    """Volcano helper should produce a non-empty DataFrame with expected columns."""
    feat1 = make_feature("TP53", "gene", fc=1.5, p=0.001)
    feat2 = make_feature("TNF", "gene", fc=-0.5, p=0.2)  # below default fc/p thresholds but still valid
    ds = make_dataset("ALS DS", "transcriptomics", ["ALS"], [feat1, feat2])
    dummy_session = DummySession({Dataset: [ds]})

    df = _get_volcano_data(dummy_session, dataset_id=str(ds.id), feature_type="gene")

    assert not df.empty
    assert set(df.columns) == {"feature", "log2FC", "pvalue", "-log10p"}
    # Ensure our features are present
    assert set(df["feature"]) == {"TP53", "TNF"}
    # -log10p should be positive for p < 1
    assert (df["-log10p"] > 0).all()


def test_volcano_empty_when_no_matching_dataset():
    """If no dataset is found, _get_volcano_data should return an empty DataFrame."""
    dummy_session = DummySession({Dataset: []})
    df = _get_volcano_data(dummy_session, dataset_id="nonexistent", feature_type="all")
    assert df.empty


def test_heatmap_build_and_cluster_multiple_datasets():
    """Heatmap matrix and clustering should handle multiple datasets with overlapping features."""
    f1 = make_feature("F1", ftype="gene", fc=1.0)
    f2 = make_feature("F2", ftype="gene", fc=2.0)
    f3 = make_feature("F3", ftype="gene", fc=3.0)

    ds1 = make_dataset("DS1", "lipidomics", ["ALS"], [f1, f2])
    ds2 = make_dataset("DS2", "lipidomics", ["ALS"], [f2, f3])

    dummy_session = DummySession({Dataset: [ds1, ds2]})

    mat = build_heatmap_matrix(dummy_session, [str(ds1.id), str(ds2.id)])
    assert not mat.empty
    # Expect 3 features (F1, F2, F3) and 2 datasets
    assert mat.shape == (3, 2)
    # After clustering, shape should remain the same
    clustered = _cluster_matrix(mat.fillna(0.0))
    assert clustered.shape == mat.shape


def test_heatmap_empty_when_no_datasets():
    """Heatmap matrix builder should return empty DataFrame when no dataset IDs are provided."""
    dummy_session = DummySession({Dataset: []})
    mat = build_heatmap_matrix(dummy_session, [])
    assert mat.empty


def test_pca_matrix_with_program_datasets():
    """
    PCA feature matrix should be constructed correctly from program-linked datasets.
    We only test the data matrix construction; PCA fitting is covered indirectly.

    Note: _build_matrix expects dicts from _load_dataset_info, not ORM objects.
    """
    # Create dataset dicts as _load_dataset_info would produce
    ds1_id = str(uuid4())
    ds2_id = str(uuid4())

    datasets = [
        {
            "id": ds1_id,
            "name": "DS1",
            "omics_type": "transcriptomics",
            "disease": ["ALS"],
            "programs": ["ALS Program"],
            "features": {"G1": 1.0, "G2": -1.5},
        },
        {
            "id": ds2_id,
            "name": "DS2",
            "omics_type": "transcriptomics",
            "disease": ["ALS"],
            "programs": ["ALS Program"],
            "features": {"G1": 0.5},  # G2 missing -> will be 0.0
        },
    ]

    mat = build_pca_matrix(datasets)
    # Rows = datasets, columns = union of feature names
    assert list(mat.index) == [ds1_id, ds2_id]
    assert set(mat.columns) == {"G1", "G2"}
    # Values should be numeric
    assert mat.dtypes.unique().tolist() == [mat.dtypes[0]]

    # PCA should run without error even if number of features < number of components requested
    from sklearn.decomposition import PCA

    n_components = min(2, mat.shape[1])
    model = PCA(n_components=n_components)
    coords = model.fit_transform(mat)
    assert coords.shape[0] == mat.shape[0]


def test_network_build_with_signature():
    """Network helper should produce nodes and edges for a signature with features."""
    f1 = Feature(id=uuid4(), name="F1", feature_type="gene")
    f2 = Feature(id=uuid4(), name="F2", feature_type="gene")

    sig = Signature(id=uuid4(), name="Sig1")
    sig.features = [f1, f2]

    dummy_session = DummySession({Signature: [sig]})

    nodes, edges = _build_network(dummy_session, signature_id=str(sig.id))

    assert set(nodes) == {"F1", "F2"}
    # With two nodes, there should be exactly one edge connecting them
    assert len(edges) == 1
    assert set(edges[0]) == {"F1", "F2"}


def test_network_empty_when_no_features():
    """Network helper should return empty lists when signature has no features."""
    sig = Signature(id=uuid4(), name="EmptySig")
    sig.features = []
    dummy_session = DummySession({Signature: [sig]})

    nodes, edges = _build_network(dummy_session, signature_id=str(sig.id))
    assert nodes == []
    assert edges == []


