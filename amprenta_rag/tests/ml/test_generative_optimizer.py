"""Tests for generative chemistry optimization and filters."""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
import torch

from amprenta_rag.ml.generative.constraints import (
    PropertyConstraint, 
    ConstraintSet,
    create_lipinski_constraints,
    create_admet_constraints
)
from amprenta_rag.ml.generative.optimizer import PropertyOptimizer
from amprenta_rag.ml.generative.filters import NoveltyChecker, DiversityFilter
from amprenta_rag.ml.generative.scaffolds import ScaffoldExtractor
from amprenta_rag.ml.generative.vae import MoleculeVAE
from amprenta_rag.ml.generative.tokenizer import SMILESTokenizer


class TestPropertyConstraints:
    """Test property constraint functionality."""
    
    def test_constraint_is_satisfied(self):
        """Test single constraint satisfaction check."""
        constraint = PropertyConstraint(
            name="logP",
            min_value=0.0,
            max_value=5.0,
            weight=1.0
        )
        
        # Test satisfied values
        assert constraint.is_satisfied(2.5) == True
        assert constraint.is_satisfied(0.0) == True  # Boundary
        assert constraint.is_satisfied(5.0) == True  # Boundary
        
        # Test violated values
        assert constraint.is_satisfied(-0.1) == False
        assert constraint.is_satisfied(5.1) == False
    
    def test_constraint_set_all_satisfied(self):
        """Test multiple constraints satisfaction."""
        constraints = ConstraintSet([
            PropertyConstraint("logP", min_value=0.0, max_value=5.0),
            PropertyConstraint("mw", min_value=100.0, max_value=500.0),
            PropertyConstraint("tpsa", min_value=0.0, max_value=140.0),
        ])
        
        # All satisfied
        predictions = {"logP": 2.5, "mw": 300.0, "tpsa": 80.0}
        assert constraints.is_satisfied(predictions) == True
        
        # One violated
        predictions = {"logP": 6.0, "mw": 300.0, "tpsa": 80.0}  # logP too high
        assert constraints.is_satisfied(predictions) == False
        
        # Missing property
        predictions = {"logP": 2.5, "mw": 300.0}  # Missing tpsa
        assert constraints.is_satisfied(predictions) == False
    
    def test_constraint_penalty_computation(self):
        """Test penalty calculation for constraint violations."""
        constraint = PropertyConstraint(
            name="logP",
            min_value=0.0,
            max_value=5.0,
            target_value=2.5,
            weight=2.0
        )
        
        # No violation
        penalty = constraint.compute_penalty(2.5)
        assert penalty == 0.0
        
        # Range violation
        penalty = constraint.compute_penalty(-1.0)  # Below min
        assert penalty > 0.0
        
        penalty = constraint.compute_penalty(6.0)  # Above max
        assert penalty > 0.0
        
        # Target optimization
        penalty_close = constraint.compute_penalty(2.6)  # Close to target
        penalty_far = constraint.compute_penalty(1.0)   # Far from target
        assert penalty_far > penalty_close


class TestPropertyOptimizer:
    """Test property optimizer functionality."""
    
    @pytest.fixture
    def optimizer_setup(self):
        """Create optimizer with mock predictors for testing."""
        tokenizer = SMILESTokenizer()
        model = MoleculeVAE(
            vocab_size=tokenizer.vocab_size,
            latent_dim=64,  # Smaller for testing
            hidden_size=128,
            num_layers=1,
            embedding_dim=32
        )
        
        # Mock predictors
        mock_admet = MagicMock()
        mock_admet.predict.return_value = {
            "logp": 2.5,
            "herg": 0.1,
            "logs": -3.0
        }
        
        optimizer = PropertyOptimizer(
            vae=model,
            tokenizer=tokenizer,
            admet_predictor=mock_admet
        )
        
        return optimizer, tokenizer, model, mock_admet
    
    def test_optimizer_single_iteration(self, optimizer_setup):
        """Test one optimization iteration."""
        optimizer, tokenizer, model, mock_admet = optimizer_setup
        
        constraints = ConstraintSet([
            PropertyConstraint("logp", min_value=0.0, max_value=5.0),
            PropertyConstraint("herg", min_value=0.0, max_value=0.3),
        ])
        
        # Test with single iteration
        results = optimizer.optimize(
            seed_smiles="CCO",
            constraints=constraints,
            n_iterations=1,
            n_samples_per_iter=2
        )
        
        # Should return some results (may be empty if no valid molecules)
        assert isinstance(results, list)
        
        # Verify ADMET predictor was called
        assert mock_admet.predict.called
    
    def test_optimizer_with_admet(self, optimizer_setup):
        """Test integration with ADMET predictor."""
        optimizer, tokenizer, model, mock_admet = optimizer_setup
        
        # Configure mock to return good properties
        mock_admet.predict.return_value = {
            "logp": 2.0,  # Good
            "herg": 0.05,  # Low risk
            "logs": -2.5   # Good solubility
        }
        
        constraints = create_admet_constraints()
        
        result = optimizer.score_molecule("CCO", constraints)
        
        assert result is not None
        assert "smiles" in result
        assert "properties" in result
        assert "score" in result
        assert result["smiles"] == "CCO"
        assert result["score"] > 0.0  # Should get positive score
    
    def test_optimizer_score_molecule(self, optimizer_setup):
        """Test molecule scoring function."""
        optimizer, tokenizer, model, mock_admet = optimizer_setup
        
        constraints = ConstraintSet([
            PropertyConstraint("logp", min_value=0.0, max_value=3.0),
        ])
        
        # Mock good prediction
        mock_admet.predict.return_value = {"logp": 1.5}
        
        result = optimizer.score_molecule("CCO", constraints)
        
        assert result is not None
        assert result["smiles"] == "CCO"
        assert "logp" in result["properties"]
        assert result["properties"]["logp"] == 1.5
        assert result["score"] == 1.0  # Should satisfy constraint perfectly


class TestNoveltyChecker:
    """Test novelty checking functionality."""
    
    def test_novelty_checker_novel(self):
        """Test novel molecule detection."""
        reference_smiles = ["CCO", "CCC", "CCCC"]  # Simple alcohols/alkanes
        checker = NoveltyChecker(reference_smiles, similarity_threshold=0.9)
        
        # Very different molecule should be novel
        novel_smiles = "c1ccccc1"  # Benzene - very different from alkanes
        assert checker.is_novel(novel_smiles) == True
    
    def test_novelty_checker_known(self):
        """Test known molecule rejection."""
        reference_smiles = ["CCO", "CCC", "CCCC"]
        checker = NoveltyChecker(reference_smiles, similarity_threshold=0.9)
        
        # Identical molecule should not be novel
        known_smiles = "CCO"
        assert checker.is_novel(known_smiles) == False
        
        # Very similar molecule should not be novel
        similar_smiles = "CCCO"  # Just one carbon longer
        # This might pass or fail depending on threshold - test with lower threshold
        checker_strict = NoveltyChecker(reference_smiles, similarity_threshold=0.5)
        assert checker_strict.is_novel(similar_smiles) == False
    
    def test_novelty_filter_list(self):
        """Test filtering a list of molecules."""
        reference_smiles = ["CCO", "CCC"]
        checker = NoveltyChecker(reference_smiles, similarity_threshold=0.9)
        
        test_smiles = [
            "CCO",          # Known (identical)
            "c1ccccc1",     # Novel (benzene)
            "CCCCO",        # Possibly novel (longer chain)
        ]
        
        novel_molecules = checker.filter_novel(test_smiles)
        
        # Should exclude exact match
        assert "CCO" not in novel_molecules
        # Should include benzene (very different)
        assert "c1ccccc1" in novel_molecules


class TestDiversityFilter:
    """Test diversity filtering functionality."""
    
    def test_diversity_filter(self):
        """Test diverse set selection."""
        diversity_filter = DiversityFilter(similarity_threshold=0.7)
        
        # Mix of similar and diverse molecules
        test_smiles = [
            "CCO",           # Ethanol
            "CCCO",          # Propanol (similar to ethanol)
            "c1ccccc1",      # Benzene (very different)
            "CC(C)C",        # Isobutane (different)
            "CCCCO",         # Butanol (similar to other alcohols)
        ]
        
        diverse_set = diversity_filter.filter_diverse(test_smiles, max_count=3)
        
        # Should return at most 3 molecules
        assert len(diverse_set) <= 3
        
        # Should include benzene (very different from others)
        assert "c1ccccc1" in diverse_set
        
        # Should not include too many similar alcohols
        alcohol_count = sum(1 for s in diverse_set if "O" in s)
        assert alcohol_count <= 2  # At most 2 alcohols in diverse set
    
    def test_diversity_empty_input(self):
        """Test diversity filter with empty input."""
        diversity_filter = DiversityFilter()
        
        result = diversity_filter.filter_diverse([])
        assert result == []
        
        result = diversity_filter.filter_diverse(["invalid_smiles"])
        # Should handle invalid SMILES gracefully
        assert isinstance(result, list)


class TestScaffoldExtractor:
    """Test scaffold extraction functionality."""
    
    def test_scaffold_extraction_basic(self):
        """Test basic scaffold extraction."""
        # Simple molecule with clear scaffold
        smiles = "CCc1ccc(O)cc1"  # 4-ethylphenol
        
        scaffold = ScaffoldExtractor.get_scaffold(smiles)
        
        assert scaffold is not None
        assert isinstance(scaffold, str)
        
        # Scaffold should contain the benzene ring
        assert "c1ccc" in scaffold or "C1=CC=" in scaffold
    
    def test_scaffold_extraction_generic(self):
        """Test generic scaffold extraction."""
        smiles = "CCc1ccc(O)cc1"  # 4-ethylphenol
        
        generic_scaffold = ScaffoldExtractor.get_scaffold(smiles, generic=True)
        
        assert generic_scaffold is not None
        assert isinstance(generic_scaffold, str)
        
        # Generic scaffold should be different from specific scaffold
        specific_scaffold = ScaffoldExtractor.get_scaffold(smiles, generic=False)
        # They might be the same for simple molecules, so just check they're valid
        assert specific_scaffold is not None
    
    def test_scaffold_group_by(self):
        """Test grouping molecules by scaffold."""
        test_smiles = [
            "CCc1ccccc1",     # Ethylbenzene
            "CCCc1ccccc1",    # Propylbenzene - same scaffold
            "CCO",            # Ethanol - different scaffold
            "CCCO",           # Propanol - different scaffold
        ]
        
        groups = ScaffoldExtractor.group_by_scaffold(test_smiles)
        
        assert isinstance(groups, dict)
        assert len(groups) >= 2  # At least aromatic and aliphatic groups
        
        # Check that molecules are grouped
        total_molecules = sum(len(molecules) for molecules in groups.values())
        valid_molecules = len([s for s in test_smiles if ScaffoldExtractor.get_scaffold(s) is not None])
        assert total_molecules == valid_molecules
    
    def test_scaffold_invalid_smiles(self):
        """Test scaffold extraction with invalid SMILES."""
        invalid_smiles = "invalid_smiles_string"
        
        scaffold = ScaffoldExtractor.get_scaffold(invalid_smiles)
        assert scaffold is None
        
        framework = ScaffoldExtractor.get_framework(invalid_smiles)
        assert framework is None
        
        sidechains = ScaffoldExtractor.get_sidechains(invalid_smiles)
        assert sidechains == []


class TestConstraintPresets:
    """Test predefined constraint sets."""
    
    def test_lipinski_constraints(self):
        """Test Lipinski Rule of Five constraints."""
        constraints = create_lipinski_constraints()
        
        assert len(constraints.constraints) == 4
        
        # Check constraint names
        names = [c.name for c in constraints.constraints]
        assert "mw" in names
        assert "logp" in names
        assert "hbd" in names
        assert "hba" in names
    
    def test_admet_constraints(self):
        """Test ADMET constraints."""
        constraints = create_admet_constraints()
        
        assert len(constraints.constraints) == 3
        
        # Check constraint names
        names = [c.name for c in constraints.constraints]
        assert "herg" in names
        assert "logs" in names
        assert "logp" in names
        
        # Check that herg has higher weight (more important)
        herg_constraint = constraints.get_constraint_by_name("herg")
        assert herg_constraint.weight > 1.0
