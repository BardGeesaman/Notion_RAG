# Multi-Omics Signature Ingestion + Scoring - Implementation Plan

## Overview

Extend the signature system from lipid-only to full multi-omics support (genes, proteins, metabolites, lipids).

## Current Status

✅ Feature type inference module created  
✅ Signature loader extended to support multi-omics  
⏳ Component creation needs extension  
⏳ Scoring engine needs extension  
⏳ Signature pages need Modalities field  

## Implementation Steps

### Phase 1: Core Infrastructure (In Progress)
1. ✅ Feature type inference (`feature_type_inference.py`)
2. ✅ Extended signature loader (`signature_loader.py`)
3. ⏳ Fix duplicate `__post_init__` in SignatureComponent

### Phase 2: Component Creation (Next)
1. Extend `find_or_create_component_page()` to:
   - Support all feature types
   - Add "Feature Type" property
   - Link to appropriate feature database

2. Create multi-omics component linking:
   - Link gene components → Gene Features DB
   - Link protein components → Protein Features DB
   - Link metabolite components → Metabolite Features DB
   - Link lipid components → Lipid Species DB (existing)

### Phase 3: Signature Page Updates
1. Add "Modalities" multi-select field to Signature pages
2. Auto-populate from component feature types

### Phase 4: Scoring Engine Extension
1. Extend scoring to work across all omics types
2. Match features based on omics type of dataset

### Phase 5: Embedding Updates
1. Update signature embedding to include modalities
2. Include feature types in text representation

## Files to Modify/Create

1. ✅ `amprenta_rag/signatures/feature_type_inference.py` - Created
2. ✅ `amprenta_rag/signatures/signature_loader.py` - Extended
3. ⏳ `amprenta_rag/ingestion/signature_notion_crud.py` - Needs extension
4. ⏳ `amprenta_rag/ingestion/signature_linking.py` - Needs extension
5. ⏳ `amprenta_rag/ingestion/signature_matching.py` - Needs extension
6. ⏳ `amprenta_rag/ingestion/signature_embedding.py` - Needs update

