#!/usr/bin/env python3
"""Quick demo model seeder for generative chemistry."""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import torch
import torch.optim as optim
from torch.utils.data import DataLoader

from amprenta_rag.ml.generative.vae import MoleculeVAE
from amprenta_rag.ml.generative.tokenizer import SMILESTokenizer
from amprenta_rag.ml.generative.model_io import save_model
from scripts.train_generative_model import SMILESDataset, train_epoch

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Demo drug-like molecules (hand-picked for diversity and validity)
DEMO_SMILES = [
    # Simple molecules
    "CCO",  # ethanol
    "CC(C)O",  # isopropanol
    "CCC",  # propane
    "CCCC",  # butane
    "C1CCCCC1",  # cyclohexane
    
    # Aromatics
    "c1ccccc1",  # benzene
    "Cc1ccccc1",  # toluene
    "CCc1ccccc1",  # ethylbenzene
    "c1ccc(O)cc1",  # phenol
    "Cc1ccc(O)cc1",  # p-cresol
    
    # Drugs and drug-like
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",  # ibuprofen
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # caffeine
    "CCN(CC)CCNC(=O)c1ccc(N)cc1",  # procaine
    "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",  # salbutamol
    
    # Heterocycles
    "c1ccncc1",  # pyridine
    "c1cccnc1",  # pyrimidine
    "c1cnccn1",  # pyrazine
    "c1csc(n1)",  # thiazole
    "c1coc(n1)",  # oxazole
    
    # Functional groups
    "CC(=O)C",  # acetone
    "CCC(=O)O",  # propanoic acid
    "CCOC(=O)C",  # ethyl acetate
    "CCN(CC)CC",  # triethylamine
    "CC(C)(C)O",  # tert-butanol
    
    # Pharmaceuticals
    "CN(C)CCc1c[nH]c2ccc(C[C@H](C(=O)O)N)cc12",  # tryptophan derivative
    "COc1cc2c(cc1OC)C(=O)C(CC1CCN(C)CC1)C2",  # codeine-like
    "CC1=CC(=O)C=C(C)C1=O",  # quinone derivative
    "CCc1ccc(OCC(O)CNC(C)C)cc1",  # beta-blocker like
    "Cc1ccc(S(=O)(=O)Nc2ncccn2)cc1",  # sulfonamide
    
    # Natural products
    "CC(C)=CCO",  # prenol
    "CC(C)=CCC=C(C)C",  # geraniol backbone
    "CC1=CC(=O)C(C(C)C)=CC1=O",  # terpene quinone
    "COc1cc(C=CCO)ccc1O",  # eugenol
    "CC(=O)c1ccc(O)cc1",  # p-hydroxyacetophenone
    
    # Miscellaneous drug-like
    "CCOc1ccc(CC(C)N)cc1",  # amphetamine derivative
    "CN1CCC[C@H]1c2cccnc2",  # nicotine
    "CC(C)NCC(O)c1ccc(O)c(O)c1",  # catecholamine
    "Cc1ccc(C(=O)Nc2ccccc2)cc1",  # benzamide
    "CCc1ccc(C(=O)O)cc1",  # phenylpropionic acid
    
    # Additional diversity
    "C1CCC(CC1)N",  # cyclohexylamine
    "CCc1c[nH]c2ccccc12",  # tryptamine
    "CC(C)c1ccc(C)cc1",  # cymene
    "COc1ccccc1C=O",  # anisaldehyde
    "Cc1cccc(C)c1O",  # xylenol
]


def create_demo_model() -> None:
    """Create and train a small demo model."""
    logger.info("Creating demo generative chemistry model...")
    
    # Model parameters (small for fast training)
    latent_dim = 64
    hidden_size = 128
    num_layers = 1
    embedding_dim = 64
    max_length = 80
    batch_size = 16
    epochs = 10
    learning_rate = 1e-3
    
    logger.info(f"Demo SMILES: {len(DEMO_SMILES)} molecules")
    
    # Build vocabulary from demo data
    tokenizer = SMILESTokenizer()
    
    # Create dataset
    dataset = SMILESDataset(DEMO_SMILES, tokenizer, max_length=max_length)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    logger.info(f"Dataset: {len(dataset)} samples, {len(dataloader)} batches")
    
    # Create small model
    model = MoleculeVAE(
        vocab_size=tokenizer.vocab_size,
        latent_dim=latent_dim,
        hidden_size=hidden_size,
        num_layers=num_layers,
        embedding_dim=embedding_dim,
    )
    
    param_count = sum(p.numel() for p in model.parameters())
    logger.info(f"Model: {param_count:,} parameters")
    
    # Create optimizer
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    
    # Training loop
    logger.info("Training demo model...")
    start_time = time.time()
    
    for epoch in range(epochs):
        epoch_metrics = train_epoch(model, dataloader, optimizer, epoch, epochs, device="cpu")
        
        if epoch % 2 == 0:  # Log every 2 epochs
            logger.info(f"Epoch {epoch_metrics['epoch']:2d}/{epochs}: "
                       f"Loss={epoch_metrics['total_loss']:.4f}")
    
    training_time = time.time() - start_time
    logger.info(f"Training complete in {training_time:.1f}s")
    
    # Save model
    output_dir = Path("models/generative/demo_vae")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    save_model(model, tokenizer, output_dir)
    
    # Save demo metadata
    metadata = {
        "description": "Demo VAE model for generative chemistry",
        "training_smiles": DEMO_SMILES,
        "model_params": {
            "vocab_size": tokenizer.vocab_size,
            "latent_dim": latent_dim,
            "hidden_size": hidden_size,
            "num_layers": num_layers,
            "embedding_dim": embedding_dim,
        },
        "training_params": {
            "epochs": epochs,
            "batch_size": batch_size,
            "learning_rate": learning_rate,
            "max_length": max_length,
        },
        "training_time": training_time,
        "final_loss": epoch_metrics['total_loss'],
        "parameter_count": param_count,
    }
    
    with open(output_dir / "demo_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)
    
    # Test generation
    logger.info("Testing model generation...")
    model.eval()
    
    with torch.no_grad():
        # Generate a few molecules
        z_samples = torch.randn(3, latent_dim)
        generated = model.sample(z_samples, tokenizer, max_len=50)
        
        logger.info("Generated molecules:")
        for i, smiles in enumerate(generated):
            logger.info(f"  {i+1}: {smiles}")
    
    # Calculate model size
    model_files = list(output_dir.glob("*.pt"))
    model_size = sum(f.stat().st_size for f in model_files) / (1024 * 1024)
    logger.info(f"Model saved to: {output_dir}")
    logger.info(f"Model size: {model_size:.1f} MB")
    
    # Verify constraints
    if training_time > 30:
        logger.warning(f"Training took {training_time:.1f}s (target: <30s)")
    if model_size > 10:
        logger.warning(f"Model size {model_size:.1f}MB (target: <10MB)")
    
    logger.info("Demo model creation complete! 🎉")


def verify_demo_model() -> bool:
    """Verify the demo model works correctly."""
    from amprenta_rag.ml.generative.model_io import load_model
    
    model_path = Path("models/generative/demo_vae")
    
    if not model_path.exists():
        logger.error("Demo model not found")
        return False
    
    try:
        logger.info("Verifying demo model...")
        
        # Load model
        model, tokenizer = load_model(model_path)
        
        # Test basic functionality
        model.eval()
        with torch.no_grad():
            # Test encoding
            test_smiles = "CCO"
            tokens = tokenizer.tokenize(test_smiles)
            x = torch.tensor([tokens])
            mu, logvar = model.encode(x)
            
            # Test sampling
            z = torch.randn(1, model.latent_dim)
            generated = model.sample(z, tokenizer)
            
            logger.info(f"Model verification successful!")
            logger.info(f"Test input: {test_smiles}")
            logger.info(f"Generated: {generated[0]}")
            
            return True
    
    except Exception as e:
        logger.error(f"Model verification failed: {e}")
        return False


def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Seed demo generative model")
    parser.add_argument("--verify-only", action="store_true", 
                       help="Only verify existing model")
    
    args = parser.parse_args()
    
    if args.verify_only:
        success = verify_demo_model()
        exit(0 if success else 1)
    
    # Create demo model
    create_demo_model()
    
    # Verify it works
    success = verify_demo_model()
    if not success:
        logger.error("Demo model verification failed!")
        exit(1)
    
    logger.info("Demo model ready for use! 🚀")


if __name__ == "__main__":
    main()
