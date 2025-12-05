#!/usr/bin/env python3
"""
Comprehensive testing script for all major RAG system features.

Tests all implemented features with real Notion data to ensure
production readiness.
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.analysis.dataset_comparison import compare_datasets
from amprenta_rag.analysis.program_signature_maps import generate_program_signature_map
from amprenta_rag.analysis.enrichment import enrich_dataset_pathways
from amprenta_rag.chemistry.notion_integration import create_compound_feature_page
from amprenta_rag.chemistry.schema import Compound
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.reporting.evidence_report import generate_dataset_evidence_report
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type

logger = get_logger(__name__)


class TestResult:
    """Represents the result of a single test."""
    
    def __init__(self, name: str):
        self.name = name
        self.passed = False
        self.error: Optional[str] = None
        self.details: Dict = {}
    
    def success(self, details: Optional[Dict] = None):
        self.passed = True
        if details:
            self.details.update(details)
    
    def failure(self, error: str, details: Optional[Dict] = None):
        self.passed = False
        self.error = error
        if details:
            self.details.update(details)


def get_test_datasets() -> List[str]:
    """Get a list of dataset IDs with features for testing."""
    cfg = get_config()
    from amprenta_rag.clients.notion_client import notion_headers
    import requests
    
    url = f"{cfg.notion.base_url}/databases/{cfg.notion.exp_data_db_id}/query"
    payload = {"page_size": 20}
    
    try:
        resp = requests.post(url, headers=notion_headers(), json=payload, timeout=30)
        resp.raise_for_status()
        datasets = resp.json().get("results", [])
        
        dataset_ids = []
        for dataset in datasets:
            page_id = dataset.get("id", "")
            features = extract_dataset_features_by_type(page_id, use_cache=True)
            total = sum(len(f) for f in features.values())
            if total > 0:
                dataset_ids.append(page_id)
        
        return dataset_ids
    except Exception as e:
        logger.warning("Error fetching datasets: %r", e)
        return []


def test_feature_extraction() -> TestResult:
    """Test feature extraction from datasets."""
    result = TestResult("Feature Extraction")
    
    try:
        dataset_ids = get_test_datasets()
        if not dataset_ids:
            result.failure("No datasets with features found")
            return result
        
        test_dataset_id = dataset_ids[0]
        features_by_type = extract_dataset_features_by_type(
            test_dataset_id,
            use_cache=True,
        )
        
        total_features = sum(len(f) for f in features_by_type.values())
        feature_types = [ft for ft, fs in features_by_type.items() if fs]
        
        result.success({
            "dataset_id": test_dataset_id,
            "total_features": total_features,
            "feature_types": feature_types,
            "features_by_type": {ft: len(fs) for ft, fs in features_by_type.items() if fs},
        })
        
    except Exception as e:
        result.failure(str(e))
        logger.exception("Feature extraction test failed")
    
    return result


def test_pathway_enrichment() -> TestResult:
    """Test pathway enrichment with a real dataset."""
    result = TestResult("Pathway Enrichment")
    
    try:
        dataset_ids = get_test_datasets()
        if not dataset_ids:
            result.failure("No datasets with features found")
            return result
        
        test_dataset_id = dataset_ids[0]
        
        # Extract features to verify we have data
        features_by_type = extract_dataset_features_by_type(test_dataset_id, use_cache=True)
        total = sum(len(f) for f in features_by_type.values())
        
        if total == 0:
            result.failure("Dataset has no extractable features")
            return result
        
        # Test pathway enrichment (with lenient threshold for testing)
        enrichment_results = enrich_dataset_pathways(
            dataset_page_id=test_dataset_id,
            p_value_threshold=0.1,  # More lenient for testing
            pathway_sources=["KEGG"],  # Just KEGG for faster testing
        )
        
        result.success({
            "dataset_id": test_dataset_id,
            "features_tested": total,
            "pathways_found": len(enrichment_results),
        })
        
    except Exception as e:
        result.failure(str(e))
        logger.exception("Pathway enrichment test failed")
    
    return result


def test_program_signature_maps() -> TestResult:
    """Test program signature map generation."""
    result = TestResult("Program Signature Maps")
    
    try:
        cfg = get_config()
        if not hasattr(cfg.notion, "programs_db_id") or not cfg.notion.programs_db_id:
            result.failure("Programs DB not configured")
            return result
        
        from amprenta_rag.clients.notion_client import notion_headers
        import requests
        
        url = f"{cfg.notion.base_url}/databases/{cfg.notion.programs_db_id}/query"
        payload = {"page_size": 5}
        resp = requests.post(url, headers=notion_headers(), json=payload, timeout=30)
        resp.raise_for_status()
        programs = resp.json().get("results", [])
        
        if not programs:
            result.failure("No programs found for testing")
            return result
        
        test_program_id = programs[0].get("id", "")
        
        # Generate program signature map
        # Check function signature
        import inspect
        sig = inspect.signature(generate_program_signature_map)
        params = list(sig.parameters.keys())
        
        # Call with correct parameters
        if 'top_n_signatures' in params:
            program_map = generate_program_signature_map(
                program_page_id=test_program_id,
                top_n_signatures=5,
                use_cache=True,
            )
        else:
            program_map = generate_program_signature_map(
                program_page_id=test_program_id,
            )
        
        result.success({
            "program_id": test_program_id,
            "signature_scores": len(program_map.signature_scores) if hasattr(program_map, 'signature_scores') else 0,
            "top_signatures": len(program_map.top_signatures) if hasattr(program_map, 'top_signatures') else 0,
            "omics_coverage": program_map.omics_coverage.total_datasets if program_map.omics_coverage else 0,
        })
        
    except Exception as e:
        result.failure(str(e))
        logger.exception("Program signature maps test failed")
    
    return result


def test_dataset_comparison() -> TestResult:
    """Test dataset comparison functionality."""
    result = TestResult("Dataset Comparison")
    
    try:
        dataset_ids = get_test_datasets()
        
        if len(dataset_ids) < 2:
            result.failure("Need at least 2 datasets with features for comparison")
            return result
        
        # Compare datasets
        comparison = compare_datasets(
            dataset1_id=dataset_ids[0],
            dataset2_id=dataset_ids[1],
        )
        
        # Get similarity from the comparison object
        similarity = getattr(comparison, 'overall_jaccard_similarity', None) or getattr(comparison, 'jaccard_similarity', 0.0)
        shared = getattr(comparison, 'shared_features_count', None) or getattr(comparison, 'shared_features', 0)
        
        result.success({
            "dataset1_id": dataset_ids[0],
            "dataset2_id": dataset_ids[1],
            "jaccard_similarity": similarity,
            "shared_features": shared,
        })
        
    except Exception as e:
        result.failure(str(e))
        logger.exception("Dataset comparison test failed")
    
    return result


def test_evidence_reports() -> TestResult:
    """Test evidence report generation."""
    result = TestResult("Evidence Reports")
    
    try:
        dataset_ids = get_test_datasets()
        if not dataset_ids:
            result.failure("No datasets found for testing")
            return result
        
        test_dataset_id = dataset_ids[0]
        
        # Generate evidence report
        report = generate_dataset_evidence_report(
            dataset_page_id=test_dataset_id,
            top_k_chunks=20,
        )
        
        result.success({
            "dataset_id": test_dataset_id,
            "report_length": len(report.summary) if hasattr(report, 'summary') else 0,
            "has_content": hasattr(report, 'summary') and len(report.summary) > 0,
        })
        
    except Exception as e:
        result.failure(str(e))
        logger.exception("Evidence reports test failed")
    
    return result


def test_chemistry_integration() -> TestResult:
    """Test chemistry compound promotion with RDKit."""
    result = TestResult("Chemistry Integration (RDKit)")
    
    try:
        # Test compound page creation with RDKit
        test_compound = Compound(
            compound_id="TEST-RDKIT-001",
            smiles="CCO",  # Ethanol
            inchi_key=None,  # Will be generated by RDKit
            canonical_smiles=None,  # Will be generated by RDKit
            molecular_formula=None,  # Will be generated by RDKit
            molecular_weight=None,  # Will be generated by RDKit
            logp=None,  # Will be generated by RDKit
            hbd_count=None,  # Will be generated by RDKit
            hba_count=None,  # Will be generated by RDKit
            rotatable_bonds=None,  # Will be generated by RDKit
        )
        
        # Normalize using RDKit
        from amprenta_rag.chemistry.normalization import (
            normalize_smiles,
            compute_molecular_descriptors,
        )
        
        canonical, inchi_key, formula = normalize_smiles(test_compound.smiles)
        descriptors = compute_molecular_descriptors(test_compound.smiles)
        
        # Update compound with RDKit-generated values
        test_compound.canonical_smiles = canonical
        test_compound.inchi_key = inchi_key
        test_compound.molecular_formula = formula
        test_compound.molecular_weight = descriptors.get("molecular_weight")
        test_compound.logp = descriptors.get("logp")
        test_compound.hbd_count = descriptors.get("hbd_count")
        test_compound.hba_count = descriptors.get("hba_count")
        test_compound.rotatable_bonds = descriptors.get("rotatable_bonds")
        
        # Verify RDKit worked
        if not inchi_key or not canonical or not formula:
            result.failure("RDKit normalization failed - missing values")
            return result
        
        # Create Notion page
        page_id = create_compound_feature_page(test_compound)
        
        if page_id:
            result.success({
                "compound_id": test_compound.compound_id,
                "notion_page_id": page_id,
                "inchi_key": inchi_key,
                "canonical_smiles": canonical,
                "molecular_formula": formula,
                "molecular_weight": test_compound.molecular_weight,
                "logp": test_compound.logp,
            })
        else:
            result.failure("Could not create compound page (may need DB configuration)")
        
    except Exception as e:
        result.failure(str(e))
        logger.exception("Chemistry integration test failed")
    
    return result


def test_rdkit_functionality() -> TestResult:
    """Test RDKit functionality directly."""
    result = TestResult("RDKit Functionality")
    
    try:
        from amprenta_rag.chemistry.normalization import (
            normalize_smiles,
            compute_molecular_descriptors,
        )
        
        test_smiles = "CCO"  # Ethanol
        canonical, inchi_key, formula = normalize_smiles(test_smiles)
        descriptors = compute_molecular_descriptors(test_smiles)
        
        # Verify all RDKit features work
        # Note: CCO is already canonical, so it may not change
        checks = {
            "canonical_smiles": canonical is not None and len(canonical) > 0,
            "inchi_key": inchi_key is not None and len(inchi_key) == 27,
            "molecular_formula": formula is not None,
            "molecular_weight": descriptors.get("molecular_weight") is not None,
            "logp": descriptors.get("logp") is not None,
            "hbd_count": descriptors.get("hbd_count") is not None,
            "hba_count": descriptors.get("hba_count") is not None,
        }
        
        all_passed = all(checks.values())
        
        if all_passed:
            result.success({
                "inchi_key": inchi_key,
                "canonical_smiles": canonical,
                "molecular_formula": formula,
                "descriptors": {k: v for k, v in descriptors.items() if v is not None},
            })
        else:
            failed = [k for k, v in checks.items() if not v]
            result.failure(f"RDKit features failed: {', '.join(failed)}")
        
    except Exception as e:
        result.failure(str(e))
        logger.exception("RDKit functionality test failed")
    
    return result


def main():
    """Run all tests."""
    print("\n" + "=" * 80)
    print("COMPREHENSIVE FEATURE TESTING")
    print("=" * 80 + "\n")
    
    tests = [
        test_rdkit_functionality,
        test_feature_extraction,
        test_pathway_enrichment,
        test_program_signature_maps,
        test_dataset_comparison,
        test_evidence_reports,
        test_chemistry_integration,
    ]
    
    results: List[TestResult] = []
    
    for test_func in tests:
        print(f"Running: {test_func.__name__}...")
        try:
            result = test_func()
            results.append(result)
            
            if result.passed:
                print(f"  ✅ {result.name}: PASSED")
                if result.details:
                    for key, value in result.details.items():
                        if isinstance(value, (list, dict)):
                            print(f"     {key}: {len(value) if hasattr(value, '__len__') else 'N/A'}")
                        else:
                            print(f"     {key}: {value}")
            else:
                print(f"  ❌ {result.name}: FAILED")
                if result.error:
                    print(f"     Error: {result.error}")
        except Exception as e:
            print(f"  ❌ {test_func.__name__}: EXCEPTION - {e}")
            result = TestResult(test_func.__name__)
            result.failure(str(e))
            results.append(result)
        
        print()
    
    # Summary
    print("=" * 80)
    print("TEST SUMMARY")
    print("=" * 80 + "\n")
    
    passed = sum(1 for r in results if r.passed)
    total = len(results)
    
    print(f"Total Tests: {total}")
    print(f"Passed: {passed}")
    print(f"Failed: {total - passed}")
    print(f"Success Rate: {passed/total*100:.1f}%\n")
    
    if passed == total:
        print("✅ ALL TESTS PASSED!")
    else:
        print("⚠️  Some tests failed. Review errors above.")
        print("\nFailed tests:")
        for r in results:
            if not r.passed:
                print(f"  - {r.name}: {r.error}")
    
    print("\n" + "=" * 80 + "\n")


if __name__ == "__main__":
    main()

