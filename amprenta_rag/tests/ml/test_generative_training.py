"""Tests for generative chemistry training functionality."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest
import torch

from amprenta_rag.ml.generative.tokenizer import SMILESTokenizer
from amprenta_rag.ml.generative.vae import MoleculeVAE
from scripts.train_generative_model import SMILESDataset, build_vocabulary, train_epoch
from scripts.seed_generative_demo import DEMO_SMILES


class TestGenerativeTraining:
    """Test generative model training functionality."""
    
    def test_tokenizer_build_vocab(self):
        """Test vocabulary building from SMILES data."""
        test_smiles = [
            "CCO",
            "CCC", 
            "c1ccccc1",
            "CC(=O)O",
        ]
        
        tokenizer = build_vocabulary(test_smiles, min_frequency=1)
        
        # Check vocabulary contains expected tokens
        assert tokenizer.vocab_size > 4  # At least special tokens
        assert "<PAD>" in tokenizer.vocab
        assert "<SOS>" in tokenizer.vocab
        assert "<EOS>" in tokenizer.vocab
        assert "<UNK>" in tokenizer.vocab
        
        # Check common tokens are present
        assert "C" in tokenizer.vocab
        assert "O" in tokenizer.vocab
        assert "=" in tokenizer.vocab
        assert "(" in tokenizer.vocab
        assert ")" in tokenizer.vocab
        
        # Test tokenization works
        tokens = tokenizer.tokenize("CCO")
        assert len(tokens) > 0
        assert all(isinstance(token, int) for token in tokens)
    
    def test_dataset_creation(self):
        """Test SMILES dataset creation."""
        tokenizer = SMILESTokenizer()
        
        test_smiles = [
            "CCO",
            "CCC",
            "c1ccccc1",
            "invalid_smiles_that_should_be_skipped_123456789",  # Invalid
        ]
        
        dataset = SMILESDataset(test_smiles, tokenizer, max_length=50)
        
        # Should have valid SMILES only
        assert len(dataset) >= 3  # At least 3 valid SMILES
        
        # Test dataset items
        item = dataset[0]
        assert isinstance(item, torch.Tensor)
        assert item.dtype == torch.long
        assert item.shape == (50,)  # max_length
        
        # Check padding
        assert tokenizer.pad_id in item  # Should have padding tokens
    
    def test_training_one_epoch(self):
        """Test training for one epoch."""
        # Create small model and dataset
        tokenizer = SMILESTokenizer()
        model = MoleculeVAE(
            vocab_size=tokenizer.vocab_size,
            latent_dim=32,  # Small for testing
            hidden_size=64,
            num_layers=1,
            embedding_dim=32
        )
        
        # Create dataset with demo SMILES
        dataset = SMILESDataset(DEMO_SMILES[:10], tokenizer, max_length=50)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=4, shuffle=False)
        
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
        
        # Train one epoch
        metrics = train_epoch(model, dataloader, optimizer, epoch=0, total_epochs=1)
        
        # Check metrics
        assert "epoch" in metrics
        assert "total_loss" in metrics
        assert "recon_loss" in metrics
        assert "kl_loss" in metrics
        assert "beta" in metrics
        
        assert metrics["epoch"] == 1
        assert metrics["total_loss"] > 0
        assert metrics["recon_loss"] > 0
        assert metrics["kl_loss"] >= 0  # Can be zero initially
        assert 0 <= metrics["beta"] <= 1
        
        # Loss should be finite
        assert torch.isfinite(torch.tensor(metrics["total_loss"]))
        assert torch.isfinite(torch.tensor(metrics["recon_loss"]))
        assert torch.isfinite(torch.tensor(metrics["kl_loss"]))
    
    def test_demo_model_inference(self):
        """Test demo model creation and inference."""
        # Test with minimal demo data
        demo_smiles = DEMO_SMILES[:5]  # Use first 5 for speed
        
        tokenizer = SMILESTokenizer()
        model = MoleculeVAE(
            vocab_size=tokenizer.vocab_size,
            latent_dim=16,  # Very small for testing
            hidden_size=32,
            num_layers=1,
            embedding_dim=16
        )
        
        # Test inference without training
        model.eval()
        
        with torch.no_grad():
            # Test encoding
            test_smiles = "CCO"
            tokens = tokenizer.tokenize(test_smiles)
            x = torch.tensor([tokens])
            
            mu, logvar = model.encode(x)
            assert mu.shape == (1, model.latent_dim)
            assert logvar.shape == (1, model.latent_dim)
            
            # Test sampling
            z = torch.randn(2, model.latent_dim)
            generated = model.sample(z, tokenizer, max_len=30)
            
            assert len(generated) == 2
            assert all(isinstance(smiles, str) for smiles in generated)
            assert all(len(smiles) > 0 for smiles in generated)
    
    def test_model_save_load_integration(self):
        """Test model save/load integration with training."""
        from amprenta_rag.ml.generative.model_io import save_model, load_model
        
        # Create and train small model
        tokenizer = build_vocabulary(DEMO_SMILES[:5], min_frequency=1)
        model = MoleculeVAE(
            vocab_size=tokenizer.vocab_size,
            latent_dim=16,
            hidden_size=32,
            num_layers=1
        )
        
        # Train for one step
        dataset = SMILESDataset(DEMO_SMILES[:3], tokenizer, max_length=30)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=2)
        optimizer = torch.optim.Adam(model.parameters())
        
        train_epoch(model, dataloader, optimizer, epoch=0, total_epochs=1)
        
        # Save model
        with tempfile.TemporaryDirectory() as temp_dir:
            model_path = Path(temp_dir) / "test_model"
            save_model(model, tokenizer, model_path)
            
            # Load model
            loaded_model, loaded_tokenizer = load_model(model_path)
            
            # Test loaded model
            assert loaded_model.vocab_size == model.vocab_size
            assert loaded_model.latent_dim == model.latent_dim
            assert loaded_tokenizer.vocab_size == tokenizer.vocab_size
            
            # Test inference with loaded model
            loaded_model.eval()
            with torch.no_grad():
                z = torch.randn(1, loaded_model.latent_dim)
                generated = loaded_model.sample(z, loaded_tokenizer)
                assert len(generated) == 1
                assert isinstance(generated[0], str)


class TestDemoModelCreation:
    """Test demo model creation specifically."""
    
    def test_demo_smiles_validity(self):
        """Test that demo SMILES are valid."""
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not available")
        
        valid_count = 0
        for smiles in DEMO_SMILES:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                valid_count += 1
        
        # At least 90% should be valid
        validity_ratio = valid_count / len(DEMO_SMILES)
        assert validity_ratio >= 0.9, f"Only {validity_ratio:.2%} of demo SMILES are valid"
    
    def test_demo_smiles_diversity(self):
        """Test that demo SMILES have good diversity."""
        # Check length diversity
        lengths = [len(smiles) for smiles in DEMO_SMILES]
        assert min(lengths) >= 3, "Some SMILES too short"
        assert max(lengths) <= 50, "Some SMILES too long"
        
        # Check character diversity
        all_chars = set("".join(DEMO_SMILES))
        expected_chars = set("CNOSPFClBr()[]=#+-c123456789")
        
        # Should have most common chemistry characters
        common_overlap = len(all_chars & expected_chars) / len(expected_chars)
        assert common_overlap >= 0.4, f"Demo SMILES missing common chemistry chars: {common_overlap:.2%}"
    
    def test_demo_vocabulary_size(self):
        """Test that demo vocabulary is reasonable size."""
        tokenizer = build_vocabulary(DEMO_SMILES, min_frequency=1)
        
        # Should be reasonable size for demo
        assert 15 <= tokenizer.vocab_size <= 100, f"Demo vocab size {tokenizer.vocab_size} not reasonable"
        
        # Should contain essential tokens
        essential_tokens = ["C", "N", "O", "(", ")", "c", "1", "="]
        for token in essential_tokens:
            assert token in tokenizer.vocab, f"Essential token '{token}' missing from vocab"
