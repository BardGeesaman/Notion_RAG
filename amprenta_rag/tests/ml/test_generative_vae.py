"""Tests for generative VAE module."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from uuid import uuid4

import pytest
import torch

from amprenta_rag.ml.generative.tokenizer import SMILESTokenizer
from amprenta_rag.ml.generative.vae import MoleculeVAE
from amprenta_rag.ml.generative.sampling import LatentSampler
from amprenta_rag.ml.generative.model_io import save_model, load_model


class TestSMILESTokenizer:
    """Test SMILES tokenizer functionality."""
    
    def test_tokenizer_encode_decode(self):
        """Test round-trip tokenization preserves SMILES."""
        tokenizer = SMILESTokenizer()
        
        # Simple SMILES examples
        test_smiles = [
            "CCO",  # ethanol
            "C1=CC=CC=C1",  # benzene
            "CC(=O)O",  # acetic acid
            "C[C@H](N)C(=O)O",  # alanine with stereochemistry
        ]
        
        for smiles in test_smiles:
            # Encode then decode
            token_ids = tokenizer.tokenize(smiles, add_special_tokens=False)
            decoded = tokenizer.detokenize(token_ids, remove_special_tokens=False)
            
            assert decoded == smiles, f"Round-trip failed for {smiles}: got {decoded}"
    
    def test_tokenizer_special_chars(self):
        """Test handling of stereochemistry and multi-character tokens."""
        tokenizer = SMILESTokenizer()
        
        # Test multi-character tokens
        smiles_with_multi = "CCl"  # Contains "Cl"
        tokens = tokenizer._parse_tokens(smiles_with_multi)
        assert "Cl" in tokens, f"Multi-char token 'Cl' not recognized in {tokens}"
        
        # Test stereochemistry - use case where @@ appears standalone
        smiles_stereo = "C@@C"  # Contains standalone "@@"
        tokens = tokenizer._parse_tokens(smiles_stereo)
        assert "@@" in tokens, f"Stereochemistry token '@@' not recognized in {tokens}"
        
        # Test bracketed stereochemistry
        smiles_bracket = "C[C@@H](O)C"  # Contains "[C@@H]"
        tokens = tokenizer._parse_tokens(smiles_bracket)
        assert "[C@@H]" in tokens, f"Bracketed stereochemistry '[C@@H]' not recognized in {tokens}"
        
        # Test round-trip with complex SMILES
        complex_smiles = "C[C@@H](Cl)Br"
        token_ids = tokenizer.tokenize(complex_smiles, add_special_tokens=False)
        decoded = tokenizer.detokenize(token_ids, remove_special_tokens=False)
        assert decoded == complex_smiles
    
    def test_tokenizer_vocab_persistence(self):
        """Test saving and loading vocabulary."""
        tokenizer = SMILESTokenizer()
        
        with tempfile.TemporaryDirectory() as temp_dir:
            vocab_path = Path(temp_dir) / "vocab.json"
            
            # Save vocabulary
            tokenizer.save_vocab(vocab_path)
            assert vocab_path.exists()
            
            # Load into new tokenizer
            new_tokenizer = SMILESTokenizer()
            new_tokenizer.load_vocab(vocab_path)
            
            # Verify vocabularies match
            assert new_tokenizer.vocab == tokenizer.vocab
            assert new_tokenizer.vocab_size == tokenizer.vocab_size
            
            # Test tokenization consistency
            test_smiles = "CCO"
            tokens1 = tokenizer.tokenize(test_smiles)
            tokens2 = new_tokenizer.tokenize(test_smiles)
            assert tokens1 == tokens2


class TestMoleculeVAE:
    """Test VAE model functionality."""
    
    @pytest.fixture
    def model_setup(self):
        """Create model and tokenizer for testing."""
        tokenizer = SMILESTokenizer()
        model = MoleculeVAE(
            vocab_size=tokenizer.vocab_size,
            latent_dim=64,  # Smaller for testing
            hidden_size=128,
            num_layers=1,
            embedding_dim=32
        )
        return model, tokenizer
    
    def test_vae_encode_shape(self, model_setup):
        """Test encoder output shape."""
        model, tokenizer = model_setup
        
        # Create dummy input
        batch_size = 2
        seq_len = 10
        x = torch.randint(0, tokenizer.vocab_size, (batch_size, seq_len))
        
        # Encode
        mu, logvar = model.encode(x)
        
        # Check shapes
        assert mu.shape == (batch_size, model.latent_dim)
        assert logvar.shape == (batch_size, model.latent_dim)
    
    def test_vae_decode_shape(self, model_setup):
        """Test decoder output shape."""
        model, tokenizer = model_setup
        
        # Create dummy latent vectors
        batch_size = 2
        max_len = 15
        z = torch.randn(batch_size, model.latent_dim)
        
        # Decode
        logits = model.decode(z, max_len)
        
        # Check shape
        assert logits.shape == (batch_size, max_len, model.vocab_size)
    
    def test_vae_forward_loss(self, model_setup):
        """Test forward pass and loss computation."""
        model, tokenizer = model_setup
        
        # Create dummy input
        batch_size = 2
        seq_len = 10
        x = torch.randint(1, tokenizer.vocab_size, (batch_size, seq_len))  # Avoid pad token
        
        # Forward pass
        recon_logits, mu, logvar = model(x)
        
        # Check shapes
        assert recon_logits.shape == (batch_size, seq_len - 1, model.vocab_size)
        assert mu.shape == (batch_size, model.latent_dim)
        assert logvar.shape == (batch_size, model.latent_dim)
        
        # Compute losses (basic check that they're computable)
        target = x[:, 1:]  # Shift for reconstruction loss
        recon_loss = torch.nn.functional.cross_entropy(
            recon_logits.reshape(-1, model.vocab_size),
            target.reshape(-1),
            reduction='mean'
        )
        
        kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp()) / batch_size
        
        total_loss = recon_loss + kl_loss
        
        # Check losses are finite
        assert torch.isfinite(recon_loss)
        assert torch.isfinite(kl_loss) 
        assert torch.isfinite(total_loss)
    
    def test_vae_sample(self, model_setup):
        """Test SMILES generation from latent vectors."""
        model, tokenizer = model_setup
        
        # Create latent vectors
        batch_size = 2
        z = torch.randn(batch_size, model.latent_dim)
        
        # Generate SMILES
        smiles_list = model.sample(z, tokenizer, max_len=20)
        
        # Check output
        assert len(smiles_list) == batch_size
        assert all(isinstance(smiles, str) for smiles in smiles_list)


class TestLatentSampler:
    """Test latent space sampling strategies."""
    
    @pytest.fixture
    def sampler(self):
        """Create sampler for testing."""
        return LatentSampler(latent_dim=64, device="cpu")
    
    def test_sampler_random(self, sampler):
        """Test random sampling."""
        n_samples = 5
        samples = sampler.random_sample(n_samples)
        
        # Check output
        assert len(samples) == n_samples
        assert all(isinstance(z, torch.Tensor) for z in samples)
        assert all(z.shape == (sampler.latent_dim,) for z in samples)
    
    def test_sampler_interpolate(self, sampler):
        """Test linear interpolation."""
        z1 = torch.randn(sampler.latent_dim)
        z2 = torch.randn(sampler.latent_dim)
        steps = 5
        
        interpolated = sampler.interpolate(z1, z2, steps)
        
        # Check output
        assert len(interpolated) == steps
        assert all(z.shape == (sampler.latent_dim,) for z in interpolated)
        
        # Check endpoints
        assert torch.allclose(interpolated[0], z1, atol=1e-6)
        assert torch.allclose(interpolated[-1], z2, atol=1e-6)
    
    def test_sampler_perturb(self, sampler):
        """Test noise perturbation."""
        z_base = torch.randn(sampler.latent_dim)
        n_perturbations = 3
        noise_scale = 0.1
        
        perturbed = sampler.perturb(z_base, noise_scale, n_perturbations)
        
        # Check output
        assert len(perturbed) == n_perturbations
        assert all(z.shape == (sampler.latent_dim,) for z in perturbed)
        
        # Check perturbations are different from original
        for z_pert in perturbed:
            assert not torch.allclose(z_pert, z_base, atol=1e-6)


class TestModelIO:
    """Test model save/load functionality."""
    
    @pytest.fixture
    def model_setup(self):
        """Create model and tokenizer for testing."""
        tokenizer = SMILESTokenizer()
        model = MoleculeVAE(
            vocab_size=tokenizer.vocab_size,
            latent_dim=64,
            hidden_size=128,
            num_layers=1,
            embedding_dim=32
        )
        return model, tokenizer
    
    def test_model_io_roundtrip(self, model_setup):
        """Test model save and load preserves weights."""
        original_model, tokenizer = model_setup
        
        with tempfile.TemporaryDirectory() as temp_dir:
            model_path = Path(temp_dir) / "test_model"
            
            # Save model
            save_model(original_model, tokenizer, model_path)
            
            # Check files exist
            assert (model_path / "model.pt").exists()
            assert (model_path / "config.json").exists()
            assert (model_path / "vocab.json").exists()
            
            # Load model
            loaded_model, loaded_tokenizer = load_model(model_path)
            
            # Check model architecture matches
            assert loaded_model.vocab_size == original_model.vocab_size
            assert loaded_model.latent_dim == original_model.latent_dim
            assert loaded_model.hidden_size == original_model.hidden_size
            
            # Check tokenizer matches
            assert loaded_tokenizer.vocab == tokenizer.vocab
            
            # Check model weights match (test a few parameters)
            original_params = dict(original_model.named_parameters())
            loaded_params = dict(loaded_model.named_parameters())
            
            for name in ["embedding.weight", "fc_mu.weight"]:
                assert torch.allclose(
                    original_params[name], 
                    loaded_params[name], 
                    atol=1e-6
                ), f"Parameter {name} differs after save/load"
            
            # Test inference consistency
            test_input = torch.randint(0, tokenizer.vocab_size, (1, 10))
            
            original_model.eval()
            loaded_model.eval()
            
            with torch.no_grad():
                mu1, logvar1 = original_model.encode(test_input)
                mu2, logvar2 = loaded_model.encode(test_input)
                
                assert torch.allclose(mu1, mu2, atol=1e-6)
                assert torch.allclose(logvar1, logvar2, atol=1e-6)
