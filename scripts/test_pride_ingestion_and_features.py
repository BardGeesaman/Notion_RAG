#!/usr/bin/env python3
"""
Test PRIDE ingestion and feature extraction end-to-end.

Tests:
1. Harvest a PRIDE study (creates dataset in Postgres)
2. Extract protein features from the dataset
3. Verify features are linked
"""

import sys
from pathlib import Path
from uuid import UUID

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset
from amprenta_rag.ingestion.repository_feature_extraction import (
    extract_features_from_repository_dataset,
)
from amprenta_rag.ingestion.repositories.discovery import fetch_study_metadata
from amprenta_rag.ingestion.postgres_integration import create_or_update_dataset_in_postgres
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)


def test_pride_ingestion(project_id: str):
    """Test harvesting a PRIDE study and creating a dataset."""
    print("\n" + "="*70)
    print("TEST 1: PRIDE Study Ingestion")
    print("="*70 + "\n")

    print(f"Harvesting PRIDE project: {project_id}")

    try:
        # Fetch study metadata
        metadata = fetch_study_metadata(
            repository="PRIDE",
            study_id=project_id,
        )

        if not metadata:
            print(f"\n❌ Failed to fetch metadata for {project_id}")
            return None

        print(f"\n✅ Metadata fetched successfully!")
        print(f"   Title: {metadata.title}")
        print(f"   Organism: {', '.join(metadata.organism) if metadata.organism else 'N/A'}")

        # Create dataset in Postgres
        with db_session() as db:
            dataset = create_or_update_dataset_in_postgres(
                name=metadata.title,
                omics_type=OmicsType.PROTEOMICS,
                external_ids={"study_id": project_id, "accession": project_id, "repository": "PRIDE"},
                description=metadata.summary,
                organism=metadata.organism,
                disease=metadata.disease,
                db=db,
            )

            db.commit()
            db.refresh(dataset)

            print(f"\n✅ Dataset created in Postgres!")
            print(f"   Dataset ID: {dataset.id}")
            print(f"   Dataset Name: {dataset.name}")
            print(f"   Omics Type: {dataset.omics_type}")

            return dataset.id

    except Exception as e:
        print(f"\n❌ Ingestion failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_pride_feature_extraction(dataset_id: UUID, project_id: str):
    """Test extracting protein features from a PRIDE dataset."""
    print("\n" + "="*70)
    print("TEST 2: PRIDE Feature Extraction")
    print("="*70 + "\n")

    print(f"Extracting protein features for dataset: {dataset_id}")
    print(f"PRIDE project: {project_id}")

    try:
        # Extract features
        linked_count = extract_features_from_repository_dataset(
            dataset_id=dataset_id,
            repository="PRIDE",
            study_id=project_id,
            download_dir=Path("/tmp"),
        )

        if linked_count > 0:
            print(f"\n✅ Feature extraction successful!")
            print(f"   Linked {linked_count} proteins to dataset")
        else:
            print(f"\n⚠️  No features extracted (this may be OK if project has no suitable files)")
            return

        # Verify features in database
        with db_session() as db:
            dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()

            if dataset:
                # Access linked features through relationship
                features = dataset.features or []

                print(f"\n✅ Verification:")
                print(f"   Total features linked: {len(features)}")

                # Filter for proteins
                proteins = [f for f in features if f.feature_type == "protein"]

                if proteins:
                    print(f"   Proteins: {len(proteins)}")
                    print(f"\n   Sample proteins (first 10):")
                    for i, feature in enumerate(proteins[:10], 1):
                        print(f"      {i}. {feature.name}")
                else:
                    print(f"   No proteins found (features: {len(features)})")
            else:
                print(f"\n⚠️  Dataset not found in database")

    except Exception as e:
        print(f"\n❌ Feature extraction failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("PRIDE INGESTION & FEATURE EXTRACTION TEST")
    print("="*70)

    # Use a PRIDE project that we know has TSV files
    project_id = "PXD071156"  # From our earlier test

    print(f"\nUsing PRIDE project: {project_id}")
    print("(This project has TSV files suitable for protein extraction)")

    # Test 1: Ingestion
    dataset_id = test_pride_ingestion(project_id)

    if not dataset_id:
        print("\n❌ Cannot proceed with feature extraction - ingestion failed")
        return

    # Test 2: Feature Extraction
    test_pride_feature_extraction(dataset_id, project_id)

    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()

