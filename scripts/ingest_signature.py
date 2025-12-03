#!/usr/bin/env python3
"""
Ingest external lipid signature from file into Notion databases.

This script loads a signature definition from a TSV/CSV file and creates/updates
the corresponding pages in Notion:
- Lipid Signatures database
- Lipid Signature Components database
- Lipid Species database (canonical ontology)

Usage:
    python scripts/ingest_signature.py --signature-file path/to/signature.tsv
    python scripts/ingest_signature.py --signature-file signature.tsv --signature-type "Consortium" --version "1.0"
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.signature_ingestion import ingest_signature_from_file
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Ingest external lipid signature into Notion databases."
    )
    parser.add_argument(
        "--signature-file",
        required=True,
        type=Path,
        help="Path to signature TSV/CSV file",
    )
    parser.add_argument(
        "--signature-type",
        default="Literature-derived",
        choices=["Consortium", "Literature-derived", "Open Dataset", "Other"],
        help="Signature type (default: Literature-derived)",
    )
    parser.add_argument(
        "--data-ownership",
        default="Public",
        help="Data ownership (default: Public)",
    )
    parser.add_argument(
        "--version",
        default=None,
        help="Version string (optional)",
    )
    parser.add_argument(
        "--description",
        default=None,
        help="Description text (optional)",
    )
    parser.add_argument(
        "--disease-context",
        nargs="*",
        default=None,
        help="Disease context strings (optional, applied to all components)",
    )
    parser.add_argument(
        "--matrix",
        nargs="*",
        default=None,
        help="Matrix strings (optional, applied to all components)",
    )
    
    args = parser.parse_args()
    
    signature_file = args.signature_file
    
    if not signature_file.exists():
        print(f"Error: Signature file not found: {signature_file}", file=sys.stderr)
        return 1
    
    if not signature_file.is_file():
        print(f"Error: Path is not a file: {signature_file}", file=sys.stderr)
        return 1
    
    try:
        result = ingest_signature_from_file(
            signature_path=signature_file,
            signature_type=args.signature_type,
            data_ownership=args.data_ownership,
            version=args.version,
            description=args.description,
            disease_context=args.disease_context,
            matrix=args.matrix,
        )
        
        print("\n" + "=" * 80)
        print("✅ Signature Ingestion Complete")
        print("=" * 80)
        print(f"\nSignature Page ID: {result['signature_page_id']}")
        print(f"Components Created: {result['component_count']}")
        print(f"Lipid Species Created/Linked: {result['species_count']}")
        
        if result['warnings']:
            print(f"\n⚠️  Warnings ({len(result['warnings'])}):")
            for warning in result['warnings']:
                print(f"  - {warning}")
        
        print("\n" + "=" * 80)
        
        return 0 if not result['warnings'] else 1
        
    except Exception as e:
        print(f"\n❌ Error ingesting signature: {e}", file=sys.stderr)
        logger.exception("Error ingesting signature")
        return 1


if __name__ == "__main__":
    sys.exit(main())

