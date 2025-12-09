# GEO and PRIDE Repository Fixes

## Issues Identified

### GEO Repository
- Problem: Using "gds" database doesn't return results with current search terms
- Solution: Need to use direct GEO Series access or different database

### PRIDE Repository  
- Problem: API endpoints returning HTML instead of JSON
- Solution: Need to use correct API endpoints or alternative access method

## Recommendations

For now, these repositories may require:
1. Different API endpoints
2. API authentication keys
3. Alternative access methods (direct website parsing, FTP)

The current implementations work for discovery but may need study IDs to be provided directly rather than discovered automatically.

## Workaround

Users can still import studies by providing known study IDs:
- GEO: Provide GSE IDs directly (e.g., GSE12345)
- PRIDE: Provide PXD IDs directly (e.g., PXD012345)

