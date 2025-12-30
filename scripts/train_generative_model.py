#!/usr/bin/env python3
"""Training script for molecular VAE generative model."""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
from pathlib import Path
from typing import List, Optional

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from amprenta_rag.ml.generative.vae import MoleculeVAE
from amprenta_rag.ml.generative.tokenizer import SMILESTokenizer
from amprenta_rag.ml.generative.model_io import save_model, save_checkpoint

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class SMILESDataset(Dataset):
    """Dataset class for SMILES strings."""
    
    def __init__(self, smiles_list: List[str], tokenizer: SMILESTokenizer, max_length: int = 100):
        """Initialize dataset.
        
        Args:
            smiles_list: List of SMILES strings
            tokenizer: SMILES tokenizer
            max_length: Maximum sequence length
        """
        self.smiles_list = smiles_list
        self.tokenizer = tokenizer
        self.max_length = max_length
        
        # Pre-tokenize all SMILES
        self.tokenized_smiles = []
        valid_count = 0
        
        for smiles in smiles_list:
            try:
                tokens = tokenizer.tokenize(smiles, add_special_tokens=True)
                if len(tokens) <= max_length:
                    # Pad to max_length
                    padded = tokens + [tokenizer.pad_id] * (max_length - len(tokens))
                    self.tokenized_smiles.append(padded[:max_length])
                    valid_count += 1
                else:
                    logger.debug(f"Skipping SMILES too long: {smiles} ({len(tokens)} tokens)")
            except Exception as e:
                logger.debug(f"Failed to tokenize SMILES: {smiles}, error: {e}")
        
        logger.info(f"Dataset created: {valid_count}/{len(smiles_list)} valid SMILES")
    
    def __len__(self) -> int:
        """Get dataset size."""
        return len(self.tokenized_smiles)
    
    def __getitem__(self, idx: int) -> torch.Tensor:
        """Get tokenized SMILES at index."""
        return torch.tensor(self.tokenized_smiles[idx], dtype=torch.long)


def load_smiles_data(input_path: str) -> List[str]:
    """Load SMILES from file.
    
    Args:
        input_path: Path to input file (txt, csv, or json)
        
    Returns:
        List of SMILES strings
    """
    input_path = Path(input_path)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    smiles_list = []
    
    if input_path.suffix.lower() == '.txt':
        with open(input_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    smiles_list.append(line)
    
    elif input_path.suffix.lower() == '.csv':
        import pandas as pd
        df = pd.read_csv(input_path)
        
        # Look for common SMILES column names
        smiles_col = None
        for col in ['smiles', 'SMILES', 'canonical_smiles', 'smiles_string']:
            if col in df.columns:
                smiles_col = col
                break
        
        if smiles_col is None:
            # Use first column
            smiles_col = df.columns[0]
            logger.warning(f"No SMILES column found, using first column: {smiles_col}")
        
        smiles_list = df[smiles_col].dropna().astype(str).tolist()
    
    elif input_path.suffix.lower() == '.json':
        with open(input_path, 'r') as f:
            data = json.load(f)
        
        if isinstance(data, list):
            smiles_list = data
        elif isinstance(data, dict) and 'smiles' in data:
            smiles_list = data['smiles']
        else:
            raise ValueError("JSON file must contain list of SMILES or dict with 'smiles' key")
    
    else:
        raise ValueError(f"Unsupported file format: {input_path.suffix}")
    
    logger.info(f"Loaded {len(smiles_list)} SMILES from {input_path}")
    return smiles_list


def build_vocabulary(smiles_list: List[str], min_frequency: int = 1) -> SMILESTokenizer:
    """Build vocabulary from SMILES list.
    
    Args:
        smiles_list: List of SMILES strings
        min_frequency: Minimum frequency for tokens
        
    Returns:
        SMILES tokenizer with built vocabulary
    """
    logger.info("Building vocabulary from SMILES data...")
    
    # Start with default vocabulary
    tokenizer = SMILESTokenizer()
    
    # Count token frequencies
    token_counts = {}
    valid_smiles = 0
    
    for smiles in smiles_list:
        try:
            tokens = tokenizer._parse_tokens(smiles)
            for token in tokens:
                token_counts[token] = token_counts.get(token, 0) + 1
            valid_smiles += 1
        except Exception:
            continue
    
    # Filter by minimum frequency
    frequent_tokens = {token: count for token, count in token_counts.items() 
                      if count >= min_frequency}
    
    # Build new vocabulary
    vocab = {
        "<PAD>": 0,
        "<SOS>": 1,
        "<EOS>": 2,
        "<UNK>": 3,
    }
    
    # Add frequent tokens
    for token in sorted(frequent_tokens.keys()):
        if token not in vocab:
            vocab[token] = len(vocab)
    
    # Create new tokenizer with built vocabulary
    new_tokenizer = SMILESTokenizer(vocab=vocab)
    
    logger.info(f"Vocabulary built: {len(vocab)} tokens from {valid_smiles} valid SMILES")
    return new_tokenizer


def compute_loss(recon_logits: torch.Tensor, target: torch.Tensor, 
                mu: torch.Tensor, logvar: torch.Tensor, beta: float = 1.0) -> dict:
    """Compute VAE loss (reconstruction + KL divergence).
    
    Args:
        recon_logits: Reconstruction logits (batch_size, seq_len, vocab_size)
        target: Target sequences (batch_size, seq_len)
        mu: Latent means (batch_size, latent_dim)
        logvar: Latent log variances (batch_size, latent_dim)
        beta: KL annealing factor
        
    Returns:
        Dictionary with loss components
    """
    batch_size = target.size(0)
    
    # Reconstruction loss
    # Shift target for teacher forcing (remove SOS, add EOS)
    target_shifted = target[:, 1:]  # Remove SOS token
    recon_flat = recon_logits.reshape(-1, recon_logits.size(-1))
    target_flat = target_shifted.reshape(-1)
    
    recon_loss = nn.functional.cross_entropy(
        recon_flat, target_flat, 
        ignore_index=0,  # Ignore padding tokens
        reduction='mean'
    )
    
    # KL divergence loss
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp()) / batch_size
    
    # Total loss
    total_loss = recon_loss + beta * kl_loss
    
    return {
        'total_loss': total_loss,
        'recon_loss': recon_loss,
        'kl_loss': kl_loss,
        'beta': beta
    }


def train_epoch(model: MoleculeVAE, dataloader: DataLoader, optimizer: optim.Optimizer,
                epoch: int, total_epochs: int, device: str = "cpu") -> dict:
    """Train model for one epoch.
    
    Args:
        model: VAE model
        dataloader: Training data loader
        optimizer: Optimizer
        epoch: Current epoch (0-indexed)
        total_epochs: Total number of epochs
        device: Training device
        
    Returns:
        Dictionary with epoch metrics
    """
    model.train()
    
    total_loss = 0.0
    total_recon_loss = 0.0
    total_kl_loss = 0.0
    num_batches = len(dataloader)
    
    # KL annealing - start at 0 and increase to 1 over first half of training
    beta = min(1.0, (epoch + 1) / (total_epochs * 0.5))
    
    for batch_idx, batch in enumerate(dataloader):
        batch = batch.to(device)
        
        # Forward pass
        recon_logits, mu, logvar = model(batch)
        
        # Compute loss
        loss_dict = compute_loss(recon_logits, batch, mu, logvar, beta)
        
        # Backward pass
        optimizer.zero_grad()
        loss_dict['total_loss'].backward()
        
        # Gradient clipping
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        
        optimizer.step()
        
        # Accumulate losses
        total_loss += loss_dict['total_loss'].item()
        total_recon_loss += loss_dict['recon_loss'].item()
        total_kl_loss += loss_dict['kl_loss'].item()
        
        # Log progress
        if batch_idx % 10 == 0:
            logger.debug(f"Epoch {epoch+1}/{total_epochs}, Batch {batch_idx}/{num_batches}, "
                        f"Loss: {loss_dict['total_loss'].item():.4f}")
    
    return {
        'epoch': epoch + 1,
        'total_loss': total_loss / num_batches,
        'recon_loss': total_recon_loss / num_batches,
        'kl_loss': total_kl_loss / num_batches,
        'beta': beta
    }


def main():
    """Main training function."""
    parser = argparse.ArgumentParser(description="Train molecular VAE model")
    parser.add_argument("--input", "-i", required=True, help="Input SMILES file (txt/csv/json)")
    parser.add_argument("--output", "-o", required=True, help="Output model directory")
    parser.add_argument("--epochs", type=int, default=50, help="Number of training epochs")
    parser.add_argument("--batch-size", type=int, default=32, help="Batch size")
    parser.add_argument("--latent-dim", type=int, default=256, help="Latent dimension")
    parser.add_argument("--hidden-size", type=int, default=512, help="Hidden size")
    parser.add_argument("--num-layers", type=int, default=2, help="Number of GRU layers")
    parser.add_argument("--learning-rate", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--max-length", type=int, default=100, help="Maximum SMILES length")
    parser.add_argument("--checkpoint-freq", type=int, default=10, help="Checkpoint frequency")
    parser.add_argument("--device", default="auto", help="Training device (cpu/cuda/auto)")
    
    args = parser.parse_args()
    
    # Set device
    if args.device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device = args.device
    
    logger.info(f"Training on device: {device}")
    
    # Load data
    logger.info("Loading SMILES data...")
    smiles_list = load_smiles_data(args.input)
    
    # Build vocabulary
    tokenizer = build_vocabulary(smiles_list, min_frequency=1)
    
    # Create dataset
    dataset = SMILESDataset(smiles_list, tokenizer, max_length=args.max_length)
    dataloader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True)
    
    logger.info(f"Dataset: {len(dataset)} samples, {len(dataloader)} batches")
    
    # Create model
    model = MoleculeVAE(
        vocab_size=tokenizer.vocab_size,
        latent_dim=args.latent_dim,
        hidden_size=args.hidden_size,
        num_layers=args.num_layers,
    ).to(device)
    
    logger.info(f"Model: {sum(p.numel() for p in model.parameters()):,} parameters")
    
    # Create optimizer
    optimizer = optim.Adam(model.parameters(), lr=args.learning_rate)
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Training loop
    logger.info("Starting training...")
    start_time = time.time()
    
    for epoch in range(args.epochs):
        epoch_metrics = train_epoch(model, dataloader, optimizer, epoch, args.epochs, device)
        
        logger.info(f"Epoch {epoch_metrics['epoch']:3d}/{args.epochs}: "
                   f"Loss={epoch_metrics['total_loss']:.4f} "
                   f"(Recon={epoch_metrics['recon_loss']:.4f}, "
                   f"KL={epoch_metrics['kl_loss']:.4f}, β={epoch_metrics['beta']:.3f})")
        
        # Save checkpoint
        if (epoch + 1) % args.checkpoint_freq == 0:
            checkpoint_path = output_dir / f"checkpoint_epoch_{epoch+1}.pt"
            save_checkpoint(model, tokenizer, optimizer, epoch, 
                          epoch_metrics['total_loss'], checkpoint_path)
            logger.info(f"Checkpoint saved: {checkpoint_path}")
    
    # Save final model
    model.cpu()  # Move to CPU for saving
    save_model(model, tokenizer, output_dir)
    
    total_time = time.time() - start_time
    logger.info(f"Training complete! Total time: {total_time:.1f}s")
    logger.info(f"Model saved to: {output_dir}")
    
    # Save training metadata
    metadata = {
        "training_args": vars(args),
        "vocab_size": tokenizer.vocab_size,
        "dataset_size": len(dataset),
        "total_epochs": args.epochs,
        "total_time": total_time,
        "final_loss": epoch_metrics['total_loss'],
    }
    
    with open(output_dir / "training_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)
    
    logger.info("Training metadata saved")


if __name__ == "__main__":
    main()
