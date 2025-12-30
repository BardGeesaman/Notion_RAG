"""Generative Chemistry Service for orchestrating molecular generation."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Dict, List, Optional

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.generative.vae import MoleculeVAE
from amprenta_rag.ml.generative.tokenizer import SMILESTokenizer
from amprenta_rag.ml.generative.sampling import LatentSampler
from amprenta_rag.ml.generative.optimizer import PropertyOptimizer
from amprenta_rag.ml.generative.constraints import ConstraintSet, PropertyConstraint
from amprenta_rag.ml.generative.filters import NoveltyChecker, DiversityFilter
from amprenta_rag.ml.generative.scaffolds import ScaffoldExtractor
from amprenta_rag.ml.generative.model_io import load_model

logger = get_logger(__name__)


class GenerativeModelNotFoundError(Exception):
    """Raised when generative model is not found or not loaded."""
    pass


class GenerativeChemistryService:
    """Service for generative chemistry operations."""
    
    def __init__(
        self,
        model_path: Optional[str] = None,
        admet_predictor=None,
        qsar_predictor=None,
        alert_checker=None,
    ):
        """Initialize generative chemistry service.
        
        Args:
            model_path: Path to trained VAE model directory
            admet_predictor: ADMET predictor instance
            qsar_predictor: QSAR predictor instance
            alert_checker: Structural alert checker instance
        """
        # Set default model path if not provided
        if model_path is None:
            default_path = Path("models/generative/demo_vae")
            if default_path.exists():
                model_path = str(default_path)
        
        self.model_path = model_path
        self.admet_predictor = admet_predictor
        self.qsar_predictor = qsar_predictor
        self.alert_checker = alert_checker
        
        # Model components (loaded lazily)
        self._vae: Optional[MoleculeVAE] = None
        self._tokenizer: Optional[SMILESTokenizer] = None
        self._sampler: Optional[LatentSampler] = None
        self._optimizer: Optional[PropertyOptimizer] = None
        
        self._model_loaded = False
        self._model_info = {}
        
        logger.info("GenerativeChemistryService initialized")
    
    def _load_model(self) -> None:
        """Load VAE model and tokenizer."""
        if self._model_loaded:
            return
        
        if not self.model_path:
            raise GenerativeModelNotFoundError("No model path specified")
        
        model_path = Path(self.model_path)
        if not model_path.exists():
            raise GenerativeModelNotFoundError(f"Model not found at {model_path}")
        
        try:
            # Load model and tokenizer
            self._vae, self._tokenizer = load_model(model_path)
            
            # Initialize components
            self._sampler = LatentSampler(latent_dim=self._vae.latent_dim)
            self._optimizer = PropertyOptimizer(
                vae=self._vae,
                tokenizer=self._tokenizer,
                admet_predictor=self.admet_predictor,
                qsar_predictor=self.qsar_predictor,
                alert_checker=self.alert_checker,
            )
            
            # Store model info
            self._model_info = {
                "name": model_path.name,
                "version": "1.0.0",
                "latent_dim": self._vae.latent_dim,
                "vocab_size": self._vae.vocab_size,
                "status": "loaded",
                "path": str(model_path),
            }
            
            self._model_loaded = True
            logger.info(f"Successfully loaded generative model from {model_path}")
        
        except Exception as e:
            logger.error(f"Failed to load model from {model_path}: {e}")
            raise GenerativeModelNotFoundError(f"Failed to load model: {e}")
    
    def get_model_info(self) -> Dict:
        """Get information about the loaded model.
        
        Returns:
            Dictionary with model information
        """
        if not self._model_loaded:
            try:
                self._load_model()
            except GenerativeModelNotFoundError:
                return {
                    "name": "default",
                    "version": "unknown",
                    "latent_dim": 0,
                    "vocab_size": 0,
                    "status": "not_found",
                }
        
        return self._model_info.copy()
    
    def sample(
        self,
        n_samples: int = 10,
        temperature: float = 1.0,
        max_length: int = 100,
    ) -> List[Dict]:
        """Generate random molecules from latent space.
        
        Args:
            n_samples: Number of molecules to generate
            temperature: Sampling temperature
            max_length: Maximum SMILES length
            
        Returns:
            List of generated molecules with properties
        """
        self._load_model()
        
        logger.info(f"Generating {n_samples} random molecules")
        
        # Generate random latent vectors
        samples = self._sampler.random_sample(n_samples)
        
        # Decode to SMILES
        molecules = []
        for i, z in enumerate(samples):
            try:
                z_batch = z.unsqueeze(0)
                smiles_list = self._vae.sample(z_batch, self._tokenizer, max_len=max_length)
                
                if smiles_list and smiles_list[0]:
                    smiles = smiles_list[0]
                    
                    # Basic validation - check if valid SMILES
                    if self._is_valid_smiles(smiles):
                        molecules.append({
                            "smiles": smiles,
                            "properties": {},
                            "score": None,
                        })
            
            except Exception as e:
                logger.debug(f"Failed to generate molecule {i}: {e}")
                continue
        
        logger.info(f"Successfully generated {len(molecules)}/{n_samples} valid molecules")
        
        return {
            "molecules": molecules,
            "count": len(molecules)
        }
    
    def interpolate(
        self,
        smiles_start: str,
        smiles_end: str,
        steps: int = 10,
        interpolation_type: str = "linear",
    ) -> List[Dict]:
        """Interpolate between two molecules in latent space.
        
        Args:
            smiles_start: Starting molecule SMILES
            smiles_end: Ending molecule SMILES
            steps: Number of interpolation steps
            interpolation_type: Type of interpolation ("linear" or "spherical")
            
        Returns:
            List of interpolated molecules
        """
        self._load_model()
        
        logger.info(f"Interpolating between {smiles_start} and {smiles_end}")
        
        # Encode both molecules to latent space
        start_tokens = self._tokenizer.tokenize(smiles_start)
        end_tokens = self._tokenizer.tokenize(smiles_end)
        
        import torch
        start_tensor = torch.tensor([start_tokens])
        end_tensor = torch.tensor([end_tokens])
        
        self._vae.eval()
        with torch.no_grad():
            mu_start, _ = self._vae.encode(start_tensor)
            mu_end, _ = self._vae.encode(end_tensor)
            
            z_start = mu_start[0]  # Use mean, not sample
            z_end = mu_end[0]
        
        # Interpolate in latent space
        if interpolation_type == "spherical":
            interpolated_vectors = self._sampler.spherical_interpolation(z_start, z_end, steps)
        else:
            interpolated_vectors = self._sampler.interpolate(z_start, z_end, steps)
        
        # Decode interpolated vectors
        molecules = []
        for step, z in enumerate(interpolated_vectors):
            try:
                z_batch = z.unsqueeze(0)
                smiles_list = self._vae.sample(z_batch, self._tokenizer, max_len=100)
                
                if smiles_list and smiles_list[0]:
                    smiles = smiles_list[0]
                    
                    if self._is_valid_smiles(smiles):
                        molecules.append({
                            "smiles": smiles,
                            "properties": {},
                            "score": None,
                            "step": step,
                        })
            
            except Exception as e:
                logger.debug(f"Failed to decode interpolation step {step}: {e}")
                continue
        
        logger.info(f"Generated {len(molecules)} interpolated molecules")
        return molecules
    
    def optimize(
        self,
        seed_smiles: str,
        constraints: List[Dict],
        n_iterations: int = 100,
        n_samples_per_iter: int = 10,
        learning_rate: float = 0.1,
        temperature: float = 1.0,
    ) -> Dict:
        """Optimize molecules for desired properties.
        
        Args:
            seed_smiles: Starting molecule
            constraints: List of property constraints
            n_iterations: Number of optimization iterations
            n_samples_per_iter: Samples per iteration
            learning_rate: Learning rate
            temperature: Sampling temperature
            
        Returns:
            Optimization results
        """
        self._load_model()
        
        logger.info(f"Optimizing from seed: {seed_smiles}")
        
        # Convert constraints to ConstraintSet
        constraint_objects = []
        for c in constraints:
            constraint = PropertyConstraint(
                name=c["name"],
                min_value=c.get("min_value"),
                max_value=c.get("max_value"),
                target_value=c.get("target_value"),
                weight=c.get("weight", 1.0),
            )
            constraint_objects.append(constraint)
        
        constraint_set = ConstraintSet(constraint_objects)
        
        # Get seed properties
        seed_properties = self._get_molecule_properties(seed_smiles)
        
        # Run optimization
        results = self._optimizer.optimize(
            seed_smiles=seed_smiles,
            constraints=constraint_set,
            n_iterations=n_iterations,
            n_samples_per_iter=n_samples_per_iter,
            learning_rate=learning_rate,
            temperature=temperature,
        )
        
        # Find best score
        best_score = max([r["score"] for r in results]) if results else 0.0
        
        logger.info(f"Optimization complete: {len(results)} improved molecules, best score: {best_score:.3f}")
        
        return {
            "optimized": results,
            "seed_properties": seed_properties,
            "best_score": best_score,
            "iterations_completed": n_iterations,
        }
    
    def scaffold_hop(
        self,
        smiles: str,
        n_analogs: int = 20,
        preserve_scaffold: bool = True,
        similarity_threshold: float = 0.7,
    ) -> Dict:
        """Generate scaffold analogs or novel scaffolds.
        
        Args:
            smiles: Input molecule SMILES
            n_analogs: Number of analogs to generate
            preserve_scaffold: Whether to preserve scaffold
            similarity_threshold: Similarity threshold for filtering
            
        Returns:
            Scaffold hopping results
        """
        self._load_model()
        
        logger.info(f"Scaffold hopping for: {smiles}")
        
        # Extract scaffold
        scaffold = ScaffoldExtractor.get_scaffold(smiles)
        
        if preserve_scaffold:
            # Generate molecules with same scaffold (not implemented in this version)
            # This would require more sophisticated scaffold-constrained generation
            analogs = []
            logger.warning("Scaffold-preserving generation not yet implemented")
        else:
            # Generate diverse molecules and filter by novelty
            generated = self.sample(n_samples=n_analogs * 3)  # Generate more to filter
            
            # Filter for novelty
            novelty_checker = NoveltyChecker([smiles], similarity_threshold=similarity_threshold)
            novel_smiles = [m["smiles"] for m in generated]
            novel_smiles = novelty_checker.filter_novel(novel_smiles)
            
            # Filter for diversity
            diversity_filter = DiversityFilter(similarity_threshold=similarity_threshold)
            diverse_smiles = diversity_filter.filter_diverse(novel_smiles, max_count=n_analogs)
            
            # Convert back to molecule format
            analogs = []
            for smiles_analog in diverse_smiles:
                properties = self._get_molecule_properties(smiles_analog)
                analogs.append({
                    "smiles": smiles_analog,
                    "properties": properties,
                    "score": None,
                })
        
        logger.info(f"Generated {len(analogs)} scaffold analogs")
        
        return {
            "scaffold": scaffold,
            "analogs": analogs,
            "n_generated": len(analogs),
        }
    
    def _get_molecule_properties(self, smiles: str) -> Dict[str, float]:
        """Get properties for a molecule using available predictors.
        
        Args:
            smiles: Molecule SMILES
            
        Returns:
            Dictionary of predicted properties
        """
        properties = {}
        
        # ADMET predictions
        if self.admet_predictor:
            try:
                admet_props = self.admet_predictor.predict(smiles)
                if admet_props:
                    properties.update(admet_props)
            except Exception as e:
                logger.debug(f"ADMET prediction failed for {smiles}: {e}")
        
        # QSAR predictions
        if self.qsar_predictor:
            try:
                qsar_props = self.qsar_predictor.predict(smiles)
                if qsar_props:
                    properties.update(qsar_props)
            except Exception as e:
                logger.debug(f"QSAR prediction failed for {smiles}: {e}")
        
        # Structural alerts
        if self.alert_checker:
            try:
                alerts = self.alert_checker.check_alerts(smiles)
                properties["alert_count"] = len(alerts) if alerts else 0
            except Exception as e:
                logger.debug(f"Alert checking failed for {smiles}: {e}")
        
        return properties
    
    def _is_valid_smiles(self, smiles: str) -> bool:
        """Check if SMILES string is valid.
        
        Args:
            smiles: SMILES string to validate
            
        Returns:
            True if valid SMILES
        """
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except ImportError:
            # If RDKit not available, do basic validation
            return bool(smiles and smiles.strip() and len(smiles) > 0)
        except Exception:
            return False
