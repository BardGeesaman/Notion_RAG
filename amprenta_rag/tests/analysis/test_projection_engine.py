"""
Tests for projection engine (UMAP, t-SNE, PCA).

Tests dimensionality reduction with various algorithms and caching.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from amprenta_rag.analysis.projection_engine import ProjectorEngine, ProjectionParams


class TestProjectorEngine:
    """Tests for ProjectorEngine class."""

    def test_compute_pca_2d(self):
        """Test 2D PCA projection."""
        engine = ProjectorEngine()
        
        # Create random data
        np.random.seed(42)
        data = np.random.randn(50, 10)
        
        result = engine.compute_pca(data, n_components=2)
        
        assert result.algorithm_used == "pca"
        assert result.n_samples == 50
        assert result.n_features == 10
        assert len(result.coordinates) == 50
        assert len(result.coordinates[0]) == 2
        assert result.cached is False

    def test_compute_pca_3d(self):
        """Test 3D PCA projection."""
        engine = ProjectorEngine()
        
        np.random.seed(42)
        data = np.random.randn(30, 20)
        
        result = engine.compute_pca(data, n_components=3)
        
        assert len(result.coordinates[0]) == 3
        assert result.params["n_components"] == 3

    def test_compute_tsne_2d(self):
        """Test 2D t-SNE projection."""
        engine = ProjectorEngine()
        
        np.random.seed(42)
        data = np.random.randn(40, 15)
        
        result = engine.compute_tsne(data, n_components=2, perplexity=10, random_state=42)
        
        assert result.algorithm_used == "tsne"
        assert result.n_samples == 40
        assert len(result.coordinates) == 40
        assert len(result.coordinates[0]) == 2

    def test_compute_umap_2d(self):
        """Test 2D UMAP projection."""
        engine = ProjectorEngine()
        
        np.random.seed(42)
        data = np.random.randn(60, 25)
        
        result = engine.compute_umap(data, n_components=2, n_neighbors=10, random_state=42)
        
        assert result.algorithm_used == "umap"
        assert result.n_samples == 60
        assert len(result.coordinates) == 60
        assert len(result.coordinates[0]) == 2

    def test_input_dataframe(self):
        """Test projection with pandas DataFrame input."""
        engine = ProjectorEngine()
        
        # Create DataFrame
        df = pd.DataFrame(np.random.randn(30, 5), columns=["A", "B", "C", "D", "E"])
        
        result = engine.compute_pca(df, n_components=2)
        
        assert result.n_samples == 30
        assert result.n_features == 5

    def test_caching_behavior(self):
        """Test that results are cached."""
        engine = ProjectorEngine()
        
        np.random.seed(42)
        data = np.random.randn(25, 8)
        
        # First call - not cached
        result1 = engine.compute_pca(data, n_components=2)
        assert result1.cached is False
        
        # Second call - should be cached
        result2 = engine.compute_pca(data, n_components=2)
        assert result2.cached is True
        
        # Same data, different params - not cached
        result3 = engine.compute_pca(data, n_components=3)
        assert result3.cached is False

    def test_clear_cache(self):
        """Test cache clearing."""
        engine = ProjectorEngine()
        
        np.random.seed(42)
        data = np.random.randn(20, 5)
        
        # Populate cache
        result1 = engine.compute_pca(data, n_components=2)
        assert result1.cached is False
        
        result2 = engine.compute_pca(data, n_components=2)
        assert result2.cached is True
        
        # Clear cache
        engine.clear_cache()
        
        # Should not be cached after clear
        result3 = engine.compute_pca(data, n_components=2)
        assert result3.cached is False

    def test_different_algorithms_use_separate_cache(self):
        """Test that different algorithms cache separately."""
        engine = ProjectorEngine()
        
        np.random.seed(42)
        data = np.random.randn(30, 10)
        
        # Compute with PCA
        pca_result = engine.compute_pca(data, n_components=2)
        assert pca_result.cached is False
        
        # Compute with t-SNE (same data)
        tsne_result = engine.compute_tsne(data, n_components=2, perplexity=10, random_state=42)
        assert tsne_result.cached is False
        
        # Both should now be cached separately
        pca_result2 = engine.compute_pca(data, n_components=2)
        assert pca_result2.cached is True
        
        tsne_result2 = engine.compute_tsne(data, n_components=2, perplexity=10, random_state=42)
        assert tsne_result2.cached is True


class TestProjectionParams:
    """Tests for ProjectionParams validation."""

    def test_valid_params(self):
        """Test valid parameter creation."""
        params = ProjectionParams(
            algorithm="umap",
            n_components=2,
            n_neighbors=15,
            perplexity=30,
        )
        
        assert params.algorithm == "umap"
        assert params.n_components == 2

    def test_n_components_validation(self):
        """Test n_components must be 2 or 3."""
        # Valid
        params2 = ProjectionParams(algorithm="pca", n_components=2)
        params3 = ProjectionParams(algorithm="pca", n_components=3)
        
        assert params2.n_components == 2
        assert params3.n_components == 3
        
        # Invalid values should fail validation
        with pytest.raises(Exception):  # Pydantic validation error
            ProjectionParams(algorithm="pca", n_components=4)

