"""Property-guided optimization in latent space for molecular generation."""

from __future__ import annotations

import random
from typing import Dict, List, Optional, Tuple

import torch

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.generative.vae import MoleculeVAE
from amprenta_rag.ml.generative.tokenizer import SMILESTokenizer
from amprenta_rag.ml.generative.sampling import LatentSampler
from amprenta_rag.ml.generative.constraints import ConstraintSet

logger = get_logger(__name__)


class PropertyOptimizer:
    """Property-guided latent space optimization for molecular generation."""
    
    def __init__(
        self,
        vae: MoleculeVAE,
        tokenizer: SMILESTokenizer,
        admet_predictor=None,  # Optional[ADMETPredictor]
        qsar_predictor=None,   # Optional[QSARPredictor]
        alert_checker=None,    # Optional[StructuralAlertChecker]
    ):
        """Initialize property optimizer.
        
        Args:
            vae: Trained MoleculeVAE model
            tokenizer: SMILES tokenizer
            admet_predictor: ADMET predictor for property evaluation
            qsar_predictor: QSAR predictor for target activity
            alert_checker: Structural alert checker
        """
        self.vae = vae
        self.tokenizer = tokenizer
        self.admet_predictor = admet_predictor
        self.qsar_predictor = qsar_predictor
        self.alert_checker = alert_checker
        self.sampler = LatentSampler(latent_dim=vae.latent_dim)
        
        logger.info("PropertyOptimizer initialized")
    
    def optimize(
        self,
        seed_smiles: str,
        constraints: ConstraintSet,
        n_iterations: int = 100,
        n_samples_per_iter: int = 10,
        learning_rate: float = 0.1,
        temperature: float = 1.0,
    ) -> List[Dict]:
        """Optimize molecules in latent space guided by property constraints.
        
        Uses gradient-free optimization with random sampling and selection.
        
        Args:
            seed_smiles: Starting molecule SMILES
            constraints: Property constraints to optimize against
            n_iterations: Number of optimization iterations
            n_samples_per_iter: Candidates to generate per iteration
            learning_rate: Step size for latent space perturbation
            temperature: Sampling temperature (higher = more exploration)
            
        Returns:
            List of optimized molecules as dicts with keys:
            - smiles: Generated SMILES
            - properties: Predicted properties dict
            - score: Overall constraint satisfaction score
            - iteration: Iteration when found
        """
        logger.info(f"Starting optimization from seed: {seed_smiles}")
        
        # Encode seed molecule to latent space
        seed_tokens = self.tokenizer.tokenize(seed_smiles)
        seed_tensor = torch.tensor([seed_tokens])
        
        self.vae.eval()
        with torch.no_grad():
            mu, logvar = self.vae.encode(seed_tensor)
            # Use mean for deterministic starting point
            z_current = mu[0]  # (latent_dim,)
        
        best_molecules = []
        best_score = 0.0
        
        for iteration in range(n_iterations):
            # Generate candidate molecules by perturbing current latent vector
            candidates = self.sampler.perturb(
                z_current, 
                noise_scale=learning_rate * temperature,
                n_perturbations=n_samples_per_iter
            )
            
            # Add current position as a candidate
            candidates.append(z_current)
            
            # Evaluate all candidates
            candidate_results = []
            for z_candidate in candidates:
                result = self._evaluate_latent_vector(z_candidate, constraints, iteration)
                if result is not None:
                    candidate_results.append(result)
            
            if not candidate_results:
                logger.warning(f"No valid molecules generated in iteration {iteration}")
                continue
            
            # Select best candidate for next iteration
            best_candidate = max(candidate_results, key=lambda x: x["score"])
            
            # Update current position if improvement found
            if best_candidate["score"] > best_score:
                best_score = best_candidate["score"]
                # Move towards best candidate
                z_current = best_candidate["latent_vector"]
                best_molecules.append(best_candidate)
                
                logger.debug(f"Iteration {iteration}: new best score {best_score:.3f}")
            else:
                # Small random perturbation to escape local minima
                z_current = z_current + torch.randn_like(z_current) * (learning_rate * 0.1)
            
            # Reduce temperature over time (simulated annealing)
            temperature = max(0.1, temperature * 0.995)
        
        # Sort by score and return top results
        best_molecules.sort(key=lambda x: x["score"], reverse=True)
        
        # Remove latent_vector from results (not needed in output)
        for result in best_molecules:
            result.pop("latent_vector", None)
        
        logger.info(f"Optimization complete: {len(best_molecules)} improved molecules found")
        return best_molecules
    
    def _evaluate_latent_vector(
        self, 
        z: torch.Tensor, 
        constraints: ConstraintSet, 
        iteration: int
    ) -> Optional[Dict]:
        """Evaluate a latent vector by decoding and scoring.
        
        Args:
            z: Latent vector to evaluate
            constraints: Property constraints
            iteration: Current iteration number
            
        Returns:
            Result dict or None if invalid molecule
        """
        try:
            # Decode to SMILES
            z_batch = z.unsqueeze(0)  # (1, latent_dim)
            smiles_list = self.vae.sample(z_batch, self.tokenizer, max_len=100)
            
            if not smiles_list or not smiles_list[0]:
                return None
            
            smiles = smiles_list[0]
            
            # Score the molecule
            result = self.score_molecule(smiles, constraints)
            if result is None:
                return None
            
            # Add metadata
            result["iteration"] = iteration
            result["latent_vector"] = z  # Keep for optimization
            
            return result
        
        except Exception as e:
            logger.debug(f"Failed to evaluate latent vector: {e}")
            return None
    
    def score_molecule(self, smiles: str, constraints: ConstraintSet) -> Optional[Dict]:
        """Score a single molecule against property constraints.
        
        Args:
            smiles: Molecule SMILES
            constraints: Property constraints
            
        Returns:
            Dict with smiles, properties, score or None if invalid
        """
        try:
            # Predict properties
            properties = {}
            
            # ADMET predictions
            if self.admet_predictor is not None:
                try:
                    admet_props = self.admet_predictor.predict(smiles)
                    if admet_props:
                        properties.update(admet_props)
                except Exception as e:
                    logger.debug(f"ADMET prediction failed for {smiles}: {e}")
            
            # QSAR predictions
            if self.qsar_predictor is not None:
                try:
                    qsar_props = self.qsar_predictor.predict(smiles)
                    if qsar_props:
                        properties.update(qsar_props)
                except Exception as e:
                    logger.debug(f"QSAR prediction failed for {smiles}: {e}")
            
            # Structural alerts
            if self.alert_checker is not None:
                try:
                    alerts = self.alert_checker.check_alerts(smiles)
                    # Convert alerts to numeric score (fewer alerts = better)
                    properties["alert_count"] = len(alerts) if alerts else 0
                    properties["alert_score"] = max(0.0, 1.0 - len(alerts) * 0.1) if alerts else 1.0
                except Exception as e:
                    logger.debug(f"Alert checking failed for {smiles}: {e}")
            
            # Basic molecular properties (if no predictors available)
            if not properties:
                properties = self._compute_basic_properties(smiles)
            
            # Compute constraint satisfaction score
            score = constraints.compute_score(properties)
            
            return {
                "smiles": smiles,
                "properties": properties,
                "score": score,
            }
        
        except Exception as e:
            logger.debug(f"Failed to score molecule {smiles}: {e}")
            return None
    
    def _compute_basic_properties(self, smiles: str) -> Dict[str, float]:
        """Compute basic molecular properties using RDKit.
        
        Args:
            smiles: Molecule SMILES
            
        Returns:
            Dictionary of basic properties
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {}
            
            return {
                "mw": Descriptors.MolWt(mol),
                "logp": Descriptors.MolLogP(mol),
                "tpsa": Descriptors.TPSA(mol),
                "hbd": Descriptors.NumHDonors(mol),
                "hba": Descriptors.NumHAcceptors(mol),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "aromatic_rings": Descriptors.NumAromaticRings(mol),
            }
        
        except ImportError:
            logger.warning("RDKit not available for basic property calculation")
            return {}
        except Exception as e:
            logger.debug(f"Failed to compute basic properties for {smiles}: {e}")
            return {}
    
    def batch_optimize(
        self,
        seed_smiles_list: List[str],
        constraints: ConstraintSet,
        n_iterations: int = 50,
        **kwargs
    ) -> List[List[Dict]]:
        """Optimize multiple seed molecules in parallel.
        
        Args:
            seed_smiles_list: List of starting molecules
            constraints: Property constraints
            n_iterations: Number of iterations per seed
            **kwargs: Additional arguments for optimize()
            
        Returns:
            List of optimization results (one list per seed)
        """
        results = []
        
        for i, seed_smiles in enumerate(seed_smiles_list):
            logger.info(f"Optimizing seed {i+1}/{len(seed_smiles_list)}: {seed_smiles}")
            
            seed_results = self.optimize(
                seed_smiles=seed_smiles,
                constraints=constraints,
                n_iterations=n_iterations,
                **kwargs
            )
            
            results.append(seed_results)
        
        return results
    
    def random_search(
        self,
        constraints: ConstraintSet,
        n_samples: int = 1000,
        batch_size: int = 50,
    ) -> List[Dict]:
        """Random search in latent space (baseline comparison).
        
        Args:
            constraints: Property constraints
            n_samples: Total number of random samples
            batch_size: Samples to generate per batch
            
        Returns:
            List of scored molecules
        """
        logger.info(f"Starting random search with {n_samples} samples")
        
        all_results = []
        
        for batch_start in range(0, n_samples, batch_size):
            batch_end = min(batch_start + batch_size, n_samples)
            batch_samples = self.sampler.random_sample(batch_end - batch_start)
            
            for z in batch_samples:
                result = self._evaluate_latent_vector(z, constraints, iteration=-1)
                if result is not None:
                    result.pop("latent_vector", None)
                    all_results.append(result)
        
        # Sort by score
        all_results.sort(key=lambda x: x["score"], reverse=True)
        
        logger.info(f"Random search complete: {len(all_results)} valid molecules")
        return all_results
