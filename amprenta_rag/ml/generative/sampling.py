"""Latent space sampling strategies for molecular generation."""

from __future__ import annotations

from typing import List

import torch
import torch.nn.functional as F

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class LatentSampler:
    """Latent space sampling strategies for VAE-based molecular generation."""
    
    def __init__(self, latent_dim: int = 256, device: str = "cpu"):
        """Initialize latent sampler.
        
        Args:
            latent_dim: Dimensionality of latent space
            device: Device for tensor operations
        """
        self.latent_dim = latent_dim
        self.device = device
        
        logger.info(f"Initialized LatentSampler: latent_dim={latent_dim}, device={device}")
    
    def random_sample(self, n: int = 10) -> List[torch.Tensor]:
        """Sample random points from standard normal distribution.
        
        Args:
            n: Number of samples to generate
            
        Returns:
            List of latent vectors (each tensor shape: (latent_dim,))
        """
        samples = []
        
        for _ in range(n):
            z = torch.randn(self.latent_dim, device=self.device)
            samples.append(z)
        
        logger.debug(f"Generated {n} random samples from latent space")
        return samples
    
    def interpolate(
        self, 
        z1: torch.Tensor, 
        z2: torch.Tensor, 
        steps: int = 10
    ) -> List[torch.Tensor]:
        """Linear interpolation between two latent vectors.
        
        Args:
            z1: First latent vector (latent_dim,)
            z2: Second latent vector (latent_dim,)
            steps: Number of interpolation steps (including endpoints)
            
        Returns:
            List of interpolated latent vectors
        """
        if z1.shape != z2.shape:
            raise ValueError(f"Latent vectors must have same shape: {z1.shape} vs {z2.shape}")
        
        if steps < 2:
            raise ValueError("Steps must be at least 2 for interpolation")
        
        # Create interpolation weights
        alphas = torch.linspace(0, 1, steps, device=self.device)
        
        interpolated = []
        for alpha in alphas:
            z_interp = (1 - alpha) * z1 + alpha * z2
            interpolated.append(z_interp)
        
        logger.debug(f"Generated {steps} interpolation steps between latent vectors")
        return interpolated
    
    def perturb(
        self, 
        z: torch.Tensor, 
        noise_scale: float = 0.1, 
        n_perturbations: int = 1
    ) -> List[torch.Tensor]:
        """Add Gaussian noise to a latent vector.
        
        Args:
            z: Base latent vector (latent_dim,)
            noise_scale: Standard deviation of Gaussian noise
            n_perturbations: Number of perturbed versions to generate
            
        Returns:
            List of perturbed latent vectors
        """
        perturbed = []
        
        for _ in range(n_perturbations):
            noise = torch.randn_like(z, device=self.device) * noise_scale
            z_perturbed = z + noise
            perturbed.append(z_perturbed)
        
        logger.debug(f"Generated {n_perturbations} perturbed versions with noise_scale={noise_scale}")
        return perturbed
    
    def spherical_interpolation(
        self, 
        z1: torch.Tensor, 
        z2: torch.Tensor, 
        steps: int = 10
    ) -> List[torch.Tensor]:
        """Spherical linear interpolation (SLERP) between two latent vectors.
        
        More appropriate for high-dimensional spaces than linear interpolation.
        
        Args:
            z1: First latent vector (latent_dim,)
            z2: Second latent vector (latent_dim,)
            steps: Number of interpolation steps
            
        Returns:
            List of spherically interpolated latent vectors
        """
        if z1.shape != z2.shape:
            raise ValueError(f"Latent vectors must have same shape: {z1.shape} vs {z2.shape}")
        
        if steps < 2:
            raise ValueError("Steps must be at least 2 for interpolation")
        
        # Normalize vectors
        z1_norm = F.normalize(z1, dim=0)
        z2_norm = F.normalize(z2, dim=0)
        
        # Calculate angle between vectors
        dot_product = torch.dot(z1_norm, z2_norm)
        dot_product = torch.clamp(dot_product, -1.0, 1.0)  # Numerical stability
        omega = torch.acos(dot_product)
        
        # Handle parallel vectors
        if torch.abs(omega) < 1e-6:
            logger.warning("Vectors are nearly parallel, falling back to linear interpolation")
            return self.interpolate(z1, z2, steps)
        
        sin_omega = torch.sin(omega)
        
        interpolated = []
        alphas = torch.linspace(0, 1, steps, device=self.device)
        
        for alpha in alphas:
            # SLERP formula
            coeff1 = torch.sin((1 - alpha) * omega) / sin_omega
            coeff2 = torch.sin(alpha * omega) / sin_omega
            
            z_interp = coeff1 * z1_norm + coeff2 * z2_norm
            
            # Scale back to original magnitude
            magnitude = (1 - alpha) * torch.norm(z1) + alpha * torch.norm(z2)
            z_interp = z_interp * magnitude
            
            interpolated.append(z_interp)
        
        logger.debug(f"Generated {steps} spherical interpolation steps")
        return interpolated
    
    def random_walk(
        self, 
        z_start: torch.Tensor, 
        steps: int = 10, 
        step_size: float = 0.1
    ) -> List[torch.Tensor]:
        """Perform random walk in latent space.
        
        Args:
            z_start: Starting latent vector (latent_dim,)
            steps: Number of walk steps
            step_size: Size of each random step
            
        Returns:
            List of latent vectors along the walk path
        """
        walk_path = [z_start.clone()]
        z_current = z_start.clone()
        
        for _ in range(steps):
            # Random step direction
            step = torch.randn_like(z_current, device=self.device) * step_size
            z_current = z_current + step
            walk_path.append(z_current.clone())
        
        logger.debug(f"Generated random walk with {steps} steps, step_size={step_size}")
        return walk_path
