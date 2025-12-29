"""
High-dimensional projection engine for dimensionality reduction.

Provides UMAP, t-SNE, and PCA projections with caching.
"""

from __future__ import annotations

import hashlib
import pickle
from typing import Any, Dict, Optional, Union

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Simple in-memory cache for projections
_PROJECTION_CACHE: Dict[str, Any] = {}


class ProjectionParams(BaseModel):
    """Parameters for dimensionality reduction."""

    algorithm: str = Field(..., description="umap, tsne, or pca")
    n_components: int = Field(2, ge=2, le=3, description="2D or 3D projection")
    n_neighbors: int = Field(15, ge=2, le=100, description="UMAP n_neighbors")
    min_dist: float = Field(0.1, ge=0.0, le=1.0, description="UMAP min_dist")
    perplexity: int = Field(30, ge=5, le=50, description="t-SNE perplexity")
    random_state: int = Field(42, description="Random seed for reproducibility")


class ProjectionResult(BaseModel):
    """Result of dimensionality reduction."""

    coordinates: Any  # numpy array - stored as list for JSON serialization
    algorithm_used: str
    params: Dict[str, Any]
    cached: bool = False
    n_samples: int
    n_features: int

    model_config = ConfigDict(arbitrary_types_allowed=True)


class ProjectorEngine:
    """
    High-dimensional data projection engine.

    Provides UMAP, t-SNE, and PCA for dimensionality reduction with caching.
    """

    def __init__(self):
        """Initialize projection engine."""
        logger.info("[PROJECTOR] Initialized projection engine")

    def compute_umap(
        self,
        data: Union[np.ndarray, pd.DataFrame],
        n_components: int = 2,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        random_state: int = 42,
    ) -> ProjectionResult:
        """
        Compute UMAP projection.

        Args:
            data: Input data (n_samples, n_features)
            n_components: Number of dimensions (2 or 3)
            n_neighbors: UMAP n_neighbors parameter
            min_dist: UMAP min_dist parameter
            random_state: Random seed

        Returns:
            ProjectionResult with UMAP coordinates
        """
        import umap
        
        X = self._to_numpy(data)
        params = {
            "n_components": n_components,
            "n_neighbors": n_neighbors,
            "min_dist": min_dist,
            "random_state": random_state,
        }
        
        # Check cache
        cache_key = self._get_cache_key(X, "umap", params)
        if cache_key in _PROJECTION_CACHE:
            logger.info("[PROJECTOR] UMAP cache hit")
            result = _PROJECTION_CACHE[cache_key]
            result.cached = True
            return result
        
        # Compute UMAP
        logger.info("[PROJECTOR] Computing UMAP (n_components=%d, n_neighbors=%d)", n_components, n_neighbors)
        reducer = umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=random_state,
        )
        
        coords = reducer.fit_transform(X)
        
        result = ProjectionResult(
            coordinates=coords.tolist(),
            algorithm_used="umap",
            params=params,
            cached=False,
            n_samples=X.shape[0],
            n_features=X.shape[1],
        )
        
        # Cache result
        _PROJECTION_CACHE[cache_key] = result
        
        logger.info("[PROJECTOR] UMAP complete: %s → %s", X.shape, coords.shape)
        return result

    def compute_tsne(
        self,
        data: Union[np.ndarray, pd.DataFrame],
        n_components: int = 2,
        perplexity: int = 30,
        random_state: int = 42,
    ) -> ProjectionResult:
        """
        Compute t-SNE projection.

        Args:
            data: Input data (n_samples, n_features)
            n_components: Number of dimensions (2 or 3)
            perplexity: t-SNE perplexity parameter
            random_state: Random seed

        Returns:
            ProjectionResult with t-SNE coordinates
        """
        X = self._to_numpy(data)
        params = {
            "n_components": n_components,
            "perplexity": perplexity,
            "random_state": random_state,
        }
        
        # Check cache
        cache_key = self._get_cache_key(X, "tsne", params)
        if cache_key in _PROJECTION_CACHE:
            logger.info("[PROJECTOR] t-SNE cache hit")
            result = _PROJECTION_CACHE[cache_key]
            result.cached = True
            return result
        
        # Compute t-SNE
        logger.info("[PROJECTOR] Computing t-SNE (n_components=%d, perplexity=%d)", n_components, perplexity)
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            random_state=random_state,
        )
        
        coords = tsne.fit_transform(X)
        
        result = ProjectionResult(
            coordinates=coords.tolist(),
            algorithm_used="tsne",
            params=params,
            cached=False,
            n_samples=X.shape[0],
            n_features=X.shape[1],
        )
        
        # Cache result
        _PROJECTION_CACHE[cache_key] = result
        
        logger.info("[PROJECTOR] t-SNE complete: %s → %s", X.shape, coords.shape)
        return result

    def compute_pca(
        self,
        data: Union[np.ndarray, pd.DataFrame],
        n_components: int = 2,
    ) -> ProjectionResult:
        """
        Compute PCA projection.

        Args:
            data: Input data (n_samples, n_features)
            n_components: Number of dimensions (2 or 3)

        Returns:
            ProjectionResult with PCA coordinates
        """
        X = self._to_numpy(data)
        params = {"n_components": n_components}
        
        # Check cache
        cache_key = self._get_cache_key(X, "pca", params)
        if cache_key in _PROJECTION_CACHE:
            logger.info("[PROJECTOR] PCA cache hit")
            result = _PROJECTION_CACHE[cache_key]
            result.cached = True
            return result
        
        # Compute PCA
        logger.info("[PROJECTOR] Computing PCA (n_components=%d)", n_components)
        pca = PCA(n_components=n_components)
        coords = pca.fit_transform(X)
        
        result = ProjectionResult(
            coordinates=coords.tolist(),
            algorithm_used="pca",
            params=params,
            cached=False,
            n_samples=X.shape[0],
            n_features=X.shape[1],
        )
        
        # Cache result
        _PROJECTION_CACHE[cache_key] = result
        
        logger.info("[PROJECTOR] PCA complete: %s → %s", X.shape, coords.shape)
        return result

    def _to_numpy(self, data: Union[np.ndarray, pd.DataFrame]) -> np.ndarray:
        """Convert input to numpy array."""
        if isinstance(data, pd.DataFrame):
            return data.values
        return data

    def _get_cache_key(self, X: np.ndarray, algorithm: str, params: Dict) -> str:
        """Generate cache key from data and parameters."""
        # Hash data
        data_hash = hashlib.md5(X.tobytes()).hexdigest()
        
        # Hash parameters
        params_str = str(sorted(params.items()))
        params_hash = hashlib.md5(params_str.encode()).hexdigest()
        
        return f"{algorithm}_{data_hash}_{params_hash}"

    def clear_cache(self) -> None:
        """Clear projection cache."""
        _PROJECTION_CACHE.clear()
        logger.info("[PROJECTOR] Cache cleared")

