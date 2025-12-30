"""VAE model for molecular generation."""

from __future__ import annotations

from typing import Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class MoleculeVAE(nn.Module):
    """Variational Autoencoder for molecular SMILES generation."""
    
    def __init__(
        self,
        vocab_size: int,
        latent_dim: int = 256,
        hidden_size: int = 512,
        num_layers: int = 2,
        embedding_dim: int = 128,
        dropout: float = 0.1
    ):
        """Initialize VAE model.
        
        Args:
            vocab_size: Size of vocabulary
            latent_dim: Dimensionality of latent space
            hidden_size: Hidden size for GRU layers
            num_layers: Number of GRU layers
            embedding_dim: Embedding dimension
            dropout: Dropout rate
        """
        super().__init__()
        
        self.vocab_size = vocab_size
        self.latent_dim = latent_dim
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.embedding_dim = embedding_dim
        
        # Embedding layer
        self.embedding = nn.Embedding(vocab_size, embedding_dim)
        
        # Encoder
        self.encoder_gru = nn.GRU(
            embedding_dim, 
            hidden_size, 
            num_layers,
            batch_first=True,
            dropout=dropout if num_layers > 1 else 0
        )
        
        # Latent space projections
        self.fc_mu = nn.Linear(hidden_size, latent_dim)
        self.fc_logvar = nn.Linear(hidden_size, latent_dim)
        
        # Decoder
        self.decoder_input = nn.Linear(latent_dim, hidden_size)
        self.decoder_gru = nn.GRU(
            embedding_dim,
            hidden_size,
            num_layers,
            batch_first=True,
            dropout=dropout if num_layers > 1 else 0
        )
        self.decoder_output = nn.Linear(hidden_size, vocab_size)
        
        logger.info(f"Initialized MoleculeVAE: vocab_size={vocab_size}, latent_dim={latent_dim}")
    
    def encode(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """Encode input sequences to latent parameters.
        
        Args:
            x: Input token sequences (batch_size, seq_len)
            
        Returns:
            Tuple of (mu, logvar) tensors (batch_size, latent_dim)
        """
        # Embed tokens
        embedded = self.embedding(x)  # (batch_size, seq_len, embedding_dim)
        
        # Encode with GRU
        _, hidden = self.encoder_gru(embedded)  # hidden: (num_layers, batch_size, hidden_size)
        
        # Use last layer's hidden state
        h = hidden[-1]  # (batch_size, hidden_size)
        
        # Project to latent parameters
        mu = self.fc_mu(h)  # (batch_size, latent_dim)
        logvar = self.fc_logvar(h)  # (batch_size, latent_dim)
        
        return mu, logvar
    
    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        """Reparameterization trick for sampling from latent distribution.
        
        Args:
            mu: Mean parameters (batch_size, latent_dim)
            logvar: Log variance parameters (batch_size, latent_dim)
            
        Returns:
            Sampled latent vectors (batch_size, latent_dim)
        """
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def decode(self, z: torch.Tensor, max_len: int = 100) -> torch.Tensor:
        """Decode latent vectors to token sequences.
        
        Args:
            z: Latent vectors (batch_size, latent_dim)
            max_len: Maximum sequence length to generate
            
        Returns:
            Token logits (batch_size, max_len, vocab_size)
        """
        batch_size = z.size(0)
        device = z.device
        
        # Initialize decoder hidden state
        h = self.decoder_input(z)  # (batch_size, hidden_size)
        h = h.unsqueeze(0).repeat(self.num_layers, 1, 1)  # (num_layers, batch_size, hidden_size)
        
        # Initialize with SOS token (assuming ID=1)
        input_token = torch.ones(batch_size, 1, dtype=torch.long, device=device)
        
        outputs = []
        
        for _ in range(max_len):
            # Embed current token
            embedded = self.embedding(input_token)  # (batch_size, 1, embedding_dim)
            
            # Decode with GRU
            output, h = self.decoder_gru(embedded, h)
            
            # Project to vocabulary
            logits = self.decoder_output(output)  # (batch_size, 1, vocab_size)
            outputs.append(logits)
            
            # Use argmax for next input (teacher forcing disabled in generation)
            input_token = torch.argmax(logits, dim=-1)
        
        return torch.cat(outputs, dim=1)  # (batch_size, max_len, vocab_size)
    
    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass for training.
        
        Args:
            x: Input token sequences (batch_size, seq_len)
            
        Returns:
            Tuple of (reconstruction_logits, mu, logvar)
        """
        # Encode
        mu, logvar = self.encode(x)
        
        # Sample from latent distribution
        z = self.reparameterize(mu, logvar)
        
        # Decode (using teacher forcing)
        seq_len = x.size(1)
        recon_logits = self._decode_with_teacher_forcing(z, x, seq_len)
        
        return recon_logits, mu, logvar
    
    def _decode_with_teacher_forcing(
        self, 
        z: torch.Tensor, 
        target: torch.Tensor, 
        seq_len: int
    ) -> torch.Tensor:
        """Decode with teacher forcing for training.
        
        Args:
            z: Latent vectors (batch_size, latent_dim)
            target: Target sequences (batch_size, seq_len)
            seq_len: Sequence length
            
        Returns:
            Reconstruction logits (batch_size, seq_len, vocab_size)
        """
        batch_size = z.size(0)
        device = z.device
        
        # Initialize decoder hidden state
        h = self.decoder_input(z)  # (batch_size, hidden_size)
        h = h.unsqueeze(0).repeat(self.num_layers, 1, 1)  # (num_layers, batch_size, hidden_size)
        
        # Prepare decoder input (shift target by 1)
        decoder_input = target[:, :-1]  # Remove last token
        embedded = self.embedding(decoder_input)  # (batch_size, seq_len-1, embedding_dim)
        
        # Decode all at once
        output, _ = self.decoder_gru(embedded, h)  # (batch_size, seq_len-1, hidden_size)
        logits = self.decoder_output(output)  # (batch_size, seq_len-1, vocab_size)
        
        return logits
    
    def sample(self, z: torch.Tensor, tokenizer, max_len: int = 100) -> list[str]:
        """Generate SMILES from latent vectors.
        
        Args:
            z: Latent vectors (batch_size, latent_dim)
            tokenizer: SMILES tokenizer for detokenization
            max_len: Maximum sequence length
            
        Returns:
            List of generated SMILES strings
        """
        self.eval()
        
        with torch.no_grad():
            # Decode to token logits
            logits = self.decode(z, max_len)  # (batch_size, max_len, vocab_size)
            
            # Convert to token IDs
            token_ids = torch.argmax(logits, dim=-1)  # (batch_size, max_len)
            
            # Convert to SMILES strings
            smiles_list = []
            for i in range(token_ids.size(0)):
                tokens = token_ids[i].cpu().numpy().tolist()
                smiles = tokenizer.detokenize(tokens)
                smiles_list.append(smiles)
        
        return smiles_list
