"""SMILES tokenizer for molecular generative models."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Multi-character tokens that should be treated as single units
# Order matters: longer tokens first to avoid partial matches
MULTI_CHAR_TOKENS = [
    "[C@@H]", "[C@H]", "[nH]", "[NH]", "[OH]", "[SH]", 
    "[N+]", "[O-]", "[S+]", "[n+]", "[o-]", "[s+]",
    "Cl", "Br", "Si", "Se", "@@"
]

# Common atoms and symbols for default vocabulary
DEFAULT_VOCAB_CHARS = [
    # Special tokens (must be first)
    "<PAD>", "<SOS>", "<EOS>", "<UNK>",
    # Common atoms
    "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H",
    # Ring numbers
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "0",
    # Bonds and structure
    "(", ")", "[", "]", "=", "#", "-", "+", 
    # Stereochemistry
    "@", "@@", "/", "\\",
    # Charges and other
    "+", "-", ".", "%"
] + MULTI_CHAR_TOKENS


class SMILESTokenizer:
    """Character-level SMILES tokenizer with vocabulary persistence."""
    
    def __init__(self, vocab: Optional[Dict[str, int]] = None):
        """Initialize tokenizer with vocabulary.
        
        Args:
            vocab: Token to ID mapping. If None, uses default vocabulary.
        """
        if vocab is None:
            vocab = {token: idx for idx, token in enumerate(DEFAULT_VOCAB_CHARS)}
        
        self.vocab = vocab
        self.reverse_vocab = {idx: token for token, idx in vocab.items()}
        
        # Special token IDs
        self.pad_id = vocab.get("<PAD>", 0)
        self.sos_id = vocab.get("<SOS>", 1) 
        self.eos_id = vocab.get("<EOS>", 2)
        self.unk_id = vocab.get("<UNK>", 3)
        
        logger.info(f"Initialized SMILES tokenizer with vocabulary size: {len(self.vocab)}")
    
    def _parse_tokens(self, smiles: str) -> List[str]:
        """Parse SMILES string into tokens, handling multi-character tokens.
        
        Args:
            smiles: SMILES string
            
        Returns:
            List of token strings
        """
        tokens = []
        i = 0
        
        while i < len(smiles):
            # Check for multi-character tokens first
            found_multi = False
            for multi_token in MULTI_CHAR_TOKENS:
                if smiles[i:i+len(multi_token)] == multi_token:
                    tokens.append(multi_token)
                    i += len(multi_token)
                    found_multi = True
                    break
            
            if not found_multi:
                # Single character token
                tokens.append(smiles[i])
                i += 1
        
        return tokens
    
    def tokenize(self, smiles: str, add_special_tokens: bool = True) -> List[int]:
        """Convert SMILES string to token IDs.
        
        Args:
            smiles: SMILES string to tokenize
            add_special_tokens: Whether to add SOS/EOS tokens
            
        Returns:
            List of token IDs
        """
        if not smiles:
            return []
        
        # Parse into tokens
        tokens = self._parse_tokens(smiles)
        
        # Convert to IDs
        token_ids = []
        if add_special_tokens:
            token_ids.append(self.sos_id)
        
        for token in tokens:
            token_id = self.vocab.get(token, self.unk_id)
            token_ids.append(token_id)
        
        if add_special_tokens:
            token_ids.append(self.eos_id)
        
        return token_ids
    
    def detokenize(self, token_ids: List[int], remove_special_tokens: bool = True) -> str:
        """Convert token IDs back to SMILES string.
        
        Args:
            token_ids: List of token IDs
            remove_special_tokens: Whether to remove SOS/EOS/PAD tokens
            
        Returns:
            SMILES string
        """
        tokens = []
        
        for token_id in token_ids:
            if remove_special_tokens and token_id in [self.pad_id, self.sos_id, self.eos_id]:
                continue
            
            token = self.reverse_vocab.get(token_id, "<UNK>")
            if token == "<UNK>":
                logger.warning(f"Unknown token ID: {token_id}")
                continue
            
            tokens.append(token)
        
        return "".join(tokens)
    
    def save_vocab(self, path: str | Path) -> None:
        """Save vocabulary to JSON file.
        
        Args:
            path: Path to save vocabulary
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(path, "w") as f:
            json.dump(self.vocab, f, indent=2)
        
        logger.info(f"Saved vocabulary to {path}")
    
    def load_vocab(self, path: str | Path) -> None:
        """Load vocabulary from JSON file.
        
        Args:
            path: Path to vocabulary file
        """
        path = Path(path)
        
        with open(path, "r") as f:
            vocab = json.load(f)
        
        self.vocab = vocab
        self.reverse_vocab = {idx: token for token, idx in vocab.items()}
        
        # Update special token IDs
        self.pad_id = vocab.get("<PAD>", 0)
        self.sos_id = vocab.get("<SOS>", 1)
        self.eos_id = vocab.get("<EOS>", 2)
        self.unk_id = vocab.get("<UNK>", 3)
        
        logger.info(f"Loaded vocabulary from {path}, size: {len(self.vocab)}")
    
    @property
    def vocab_size(self) -> int:
        """Get vocabulary size."""
        return len(self.vocab)
