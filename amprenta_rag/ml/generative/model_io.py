"""PyTorch model save/load utilities for generative models."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Tuple

import torch

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.generative.vae import MoleculeVAE
from amprenta_rag.ml.generative.tokenizer import SMILESTokenizer

logger = get_logger(__name__)


def save_model(
    model: MoleculeVAE, 
    tokenizer: SMILESTokenizer, 
    path: str | Path
) -> None:
    """Save VAE model and tokenizer to disk.
    
    Args:
        model: Trained MoleculeVAE model
        tokenizer: SMILES tokenizer
        path: Directory path to save model files
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    
    # Save model state dict
    model_path = path / "model.pt"
    torch.save(model.state_dict(), model_path)
    
    # Save model config
    config = {
        "vocab_size": model.vocab_size,
        "latent_dim": model.latent_dim,
        "hidden_size": model.hidden_size,
        "num_layers": model.num_layers,
        "embedding_dim": model.embedding_dim,
    }
    config_path = path / "config.json"
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    
    # Save tokenizer vocabulary
    vocab_path = path / "vocab.json"
    tokenizer.save_vocab(vocab_path)
    
    logger.info(f"Saved model to {path}")


def load_model(path: str | Path) -> Tuple[MoleculeVAE, SMILESTokenizer]:
    """Load VAE model and tokenizer from disk.
    
    Args:
        path: Directory path containing model files
        
    Returns:
        Tuple of (model, tokenizer)
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"Model directory not found: {path}")
    
    # Load model config
    config_path = path / "config.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Model config not found: {config_path}")
    
    with open(config_path, "r") as f:
        config = json.load(f)
    
    # Load tokenizer
    vocab_path = path / "vocab.json"
    if not vocab_path.exists():
        raise FileNotFoundError(f"Vocabulary file not found: {vocab_path}")
    
    tokenizer = SMILESTokenizer()
    tokenizer.load_vocab(vocab_path)
    
    # Initialize model with config
    model = MoleculeVAE(**config)
    
    # Load model weights
    model_path = path / "model.pt"
    if not model_path.exists():
        raise FileNotFoundError(f"Model weights not found: {model_path}")
    
    state_dict = torch.load(model_path, map_location="cpu")
    model.load_state_dict(state_dict)
    
    logger.info(f"Loaded model from {path}")
    return model, tokenizer


def save_checkpoint(
    model: MoleculeVAE,
    tokenizer: SMILESTokenizer,
    optimizer: torch.optim.Optimizer,
    epoch: int,
    loss: float,
    path: str | Path
) -> None:
    """Save training checkpoint.
    
    Args:
        model: MoleculeVAE model
        tokenizer: SMILES tokenizer
        optimizer: Training optimizer
        epoch: Current epoch number
        loss: Current loss value
        path: Path to save checkpoint
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    checkpoint = {
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "epoch": epoch,
        "loss": loss,
        "model_config": {
            "vocab_size": model.vocab_size,
            "latent_dim": model.latent_dim,
            "hidden_size": model.hidden_size,
            "num_layers": model.num_layers,
            "embedding_dim": model.embedding_dim,
        },
        "vocab": tokenizer.vocab,
    }
    
    torch.save(checkpoint, path)
    logger.info(f"Saved checkpoint to {path} (epoch {epoch}, loss {loss:.4f})")


def load_checkpoint(path: str | Path) -> Tuple[MoleculeVAE, SMILESTokenizer, torch.optim.Optimizer, int, float]:
    """Load training checkpoint.
    
    Args:
        path: Path to checkpoint file
        
    Returns:
        Tuple of (model, tokenizer, optimizer, epoch, loss)
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {path}")
    
    checkpoint = torch.load(path, map_location="cpu")
    
    # Reconstruct tokenizer
    tokenizer = SMILESTokenizer(vocab=checkpoint["vocab"])
    
    # Reconstruct model
    model = MoleculeVAE(**checkpoint["model_config"])
    model.load_state_dict(checkpoint["model_state_dict"])
    
    # Reconstruct optimizer (need to create new instance first)
    optimizer = torch.optim.Adam(model.parameters())
    optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
    
    epoch = checkpoint["epoch"]
    loss = checkpoint["loss"]
    
    logger.info(f"Loaded checkpoint from {path} (epoch {epoch}, loss {loss:.4f})")
    return model, tokenizer, optimizer, epoch, loss
