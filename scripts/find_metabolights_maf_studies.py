#!/usr/bin/env python3
"""
Find MetaboLights studies with MAF (Metabolite Assignment File) files.

Uses the "List-then-Check" approach:
1. Get list of study IDs
2. Check each study's file list for m_*.tsv files
3. Return studies that have MAF files
"""

import sys
import time
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import requests
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

BASE_URL = "https://www.ebi.ac.uk/metabolights/ws"


def find_studies_with_maf(limit=20):
    """
    Scans the most recent public studies to find ones with 
    Metabolite Assignment Files (m_*.tsv).
    
    Args:
        limit: Maximum number of studies to check
        
    Returns:
        List of dictionaries with study_id and maf_files
    """
    print("\n" + "="*70)
    print("SEARCHING FOR METABOLIGHTS STUDIES WITH MAF FILES")
    print("="*70 + "\n")
    
    # 1. Get list of all public studies
    print("Step 1: Fetching study list...")
    try:
        resp = requests.get(f"{BASE_URL}/studies", timeout=30)
        resp.raise_for_status()
        
        data = resp.json()
        
        # The API might return different structures
        # Try common patterns
        if isinstance(data, list):
            all_studies = data
        elif isinstance(data, dict):
            all_studies = data.get('content', [])
            if not all_studies:
                all_studies = data.get('studies', [])
        else:
            print(f"⚠️  Unexpected response structure: {type(data)}")
            return []
            
        print(f"✅ Found {len(all_studies)} total studies")
        
        # Extract study IDs if they're in dictionary format
        study_ids = []
        for item in all_studies:
            if isinstance(item, str):
                study_ids.append(item)
            elif isinstance(item, dict):
                study_id = item.get('studyIdentifier') or item.get('id') or item.get('accession')
                if study_id:
                    study_ids.append(study_id)
        
        if not study_ids:
            print("⚠️  Could not extract study IDs from response")
            return []
            
        # Sort and take most recent (highest number)
        def extract_number(study_id):
            try:
                return int(study_id.replace('MTBLS', ''))
            except Exception as e:
                logger.warning("Failed to parse study number from %s: %r", study_id, e)
                return 0
        
        study_ids.sort(key=extract_number, reverse=True)
        study_ids_to_check = study_ids[:limit]
        
        print(f"✅ Checking top {len(study_ids_to_check)} most recent studies\n")
        
    except Exception as e:
        print(f"❌ Failed to fetch study list: {e}")
        return []
    
    found_studies = []
    failed_studies = []
    skipped_studies = []
    
    # 2. Iterate through studies and check for MAF files
    print("Step 2: Checking each study for MAF files...\n")
    
    for i, study_id in enumerate(study_ids_to_check, 1):
        time.sleep(1.0)  # Rate limiting - be polite!
        
        print(f"[{i}/{len(study_ids_to_check)}] Checking {study_id}...", end=" ")
        
        try:
            # Get study details to find HTTP URL (since files endpoint is broken)
            study_url = f"{BASE_URL}/studies/{study_id}"
            study_resp = requests.get(study_url, timeout=30)
            
            # CASE 1: Success
            if study_resp.status_code == 200:
                study_data = study_resp.json()
                mtbls_study = study_data.get('mtblsStudy', {})
                http_url = mtbls_study.get('studyHttpUrl', '')
                
                if not http_url:
                    print("No HTTP URL")
                    skipped_studies.append(study_id)
                    continue
                
                # Convert to HTTPS
                https_url = http_url.replace('http://ftp.', 'https://ftp.')
                
                # Check investigation file which lists all files
                import re
                inv_url = f"{https_url}/i_Investigation.txt"
                
                try:
                    inv_resp = requests.get(inv_url, timeout=30)
                    
                    if inv_resp.status_code == 200:
                        inv_content = inv_resp.text
                        
                        # 3. Check for m_*.tsv files (MAF pattern: ^m_.*\.tsv$)
                        # Look for files matching the pattern
                        maf_files = re.findall(r'(m_[^\s\t\n]+\.tsv)', inv_content, re.IGNORECASE)
                        
                        if maf_files:
                            # Filter to only .tsv files (not .txt)
                            maf_files = [f for f in maf_files if f.endswith('.tsv')]
                            
                            if maf_files:
                                print(f"✅ FOUND {len(maf_files)} MAF file(s)!")
                                print(f"   Files: {list(set(maf_files))[:3]}")
                                found_studies.append({
                                    'id': study_id,
                                    'maf_files': list(set(maf_files)),
                                    'url': https_url
                                })
                            else:
                                print("No MAF files")
                        else:
                            print("No MAF files")
                    elif inv_resp.status_code == 404:
                        print("Investigation file not found")
                    elif inv_resp.status_code >= 500:
                        print("Server error on investigation file")
                        failed_studies.append(f"{study_id} (investigation)")
                    else:
                        print(f"Investigation status {inv_resp.status_code}")
                except Exception as inv_e:
                    print(f"Error accessing investigation: {inv_e}")
                    continue
            
            # CASE 2: Known "Private/Missing" codes
            elif study_resp.status_code in [403, 404]:
                print(f"Private or missing (403/404)")
                skipped_studies.append(study_id)
            
            # CASE 3: The dreaded 500 Server Error - skip and continue
            elif study_resp.status_code >= 500:
                print(f"Server error (500) - skipping")
                failed_studies.append(study_id)
            
            else:
                print(f"Status {study_resp.status_code}")
                skipped_studies.append(study_id)
                
        except requests.exceptions.RequestException as e:
            print(f"Connection error - skipping")
            failed_studies.append(f"{study_id} (connection)")
        except Exception as e:
            print(f"Error: {e}")
            failed_studies.append(f"{study_id} (exception)")
            continue
    
    print(f"\n" + "="*70)
    print(f"SCAN COMPLETE")
    print("="*70)
    print(f"\n✅ Found: {len(found_studies)} studies with MAF files")
    print(f"⚠️  Skipped: {len(skipped_studies)} studies (private/missing)")
    print(f"❌ Failed: {len(failed_studies)} studies (server errors)")
    print("="*70 + "\n")
    
    return found_studies


def main():
    """Run the search."""
    import sys
    
    # Get limit from command line or use default
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    
    found_studies = find_studies_with_maf(limit=limit)
    
    if found_studies:
        print("Studies with MAF files:")
        for study in found_studies:
            print(f"\n  {study['id']}:")
            for maf_file in study['maf_files']:
                print(f"    - {maf_file}")
        
        print(f"\n✅ Use one of these study IDs for testing:")
        print(f"   python scripts/test_metabolights_feature_extraction.py {found_studies[0]['id']}")
    else:
        print("⚠️  No studies with MAF files found in the checked range")
        print("   Try increasing the limit or checking different study IDs")


if __name__ == "__main__":
    main()

