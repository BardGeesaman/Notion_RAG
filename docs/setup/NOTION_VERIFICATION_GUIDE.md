# Notion Verification Guide - Repository Ingestion

**Date**: 2025-12-04  
**Test Page ID**: `2bfadf61-42ab-81d8-9046-ca3d38d424a4`

---

## Where to Look in Notion

### 1. Open the Database

Go to your **"Experimental Data Assets"** database in Notion.

This is the database where all experimental datasets are stored.

---

### 2. Find the Test Page

**Option A: Search**
- Search for: `ST004217` (the MW Study ID)
- Or search for: `Lipid Alterations` (part of the title)

**Option B: Direct Link**
- Page ID: `2bfadf61-42ab-81d8-9046-ca3d38d424a4`
- Direct URL: `https://www.notion.so/2bfadf6181d89046ca3d38d424a4`
  (Replace with your workspace URL if needed)

**Option C: Look for Title**
- Title starts with: "Lipid Alterations in ASAH1-Deficient Cells: Insights into Ceramide Accumulation and Lysosomal Dysfunction"

---

## Properties to Verify

### ✅ Experiment Name
- **Expected**: "Lipid Alterations in ASAH1-Deficient Cells: Insights into Ceramide Accumulation and Lysosomal Dysfunction"
- **Status**: Should be populated

### ✅ Summary
- **Expected**: Should contain "MW Study ID: ST004217"
- **Status**: Should have the study ID and summary text

### ✅ Omics Type
- **Expected**: "Lipidomics"
- **Status**: Should be set to Lipidomics

### ✅ Signature Match Score
- **Expected**: `0.65`
- **Status**: Should show a numeric value of 0.65
- **This confirms**: Signature matching worked!

### ✅ Related Signature(s)
- **Expected**: Should show 2 signatures:
  1. ALS-CSF-Core-6Ceramides
  2. test_signature_verification
- **Status**: Should be a relation field with 2 linked pages
- **This confirms**: Signature matching and linking worked!

### ✅ Embedding IDs
- **Expected**: Should have content (297 characters)
- **Status**: Should contain embedding IDs
- **This confirms**: RAG embeddings were created and written!

---

## Page Content to Verify

### ✅ mwTab Data Section

Scroll down in the page content. You should see:

1. **Heading**: "mwTab Data" (heading level 2)
2. **Code Blocks**: Multiple code blocks containing the mwTab data
3. **This confirms**: The harvest script successfully added mwTab data to the page

---

## Feature Links to Verify

### ✅ Lipid Species Links

If you check any **Lipid Species** pages in your Lipid Species database, you should see:
- The dataset page linked in the "Experimental Data Assets" or "Datasets" relation
- **This confirms**: Feature linking worked (749 features were linked)

### ✅ Metabolite Features Links

If any metabolite features were created, they should also link back to this dataset.

---

## What Success Looks Like

If you see all of the above:

✅ **Page exists** with correct title  
✅ **Summary** contains study ID  
✅ **Omics Type** is set correctly  
✅ **Signature Match Score** is 0.65  
✅ **2 signatures** are linked  
✅ **Embedding IDs** are present  
✅ **mwTab data** is in page content  
✅ **Features** are linked (check a few Lipid Species pages)  

**Then the repository ingestion is working perfectly!** 🎉

---

## Troubleshooting

### If Signature Match Score is Missing
- Check if the property exists in your database schema
- The property name should be exactly "Signature Match Score"

### If Related Signature(s) is Empty
- Check if signatures exist in your Lipid Signatures database
- Verify the signature matching threshold (default: 0.3)

### If Embedding IDs is Empty
- The ingestion may have been interrupted
- Re-run ingestion: `python scripts/ingest_dataset.py --dataset-page-id 2bfadf61-42ab-81d8-9046-ca3d38d424a4 --force`

### If mwTab Data is Missing
- Check if the MW API returned data for ST004217
- The harvest script should have added it automatically

---

## Quick Verification Checklist

- [ ] Page exists in Experimental Data Assets database
- [ ] Title matches expected value
- [ ] Summary contains "MW Study ID: ST004217"
- [ ] Omics Type = "Lipidomics"
- [ ] Signature Match Score = 0.65
- [ ] Related Signature(s) shows 2 signatures
- [ ] Embedding IDs has content
- [ ] Page content has "mwTab Data" section
- [ ] At least one Lipid Species page links to this dataset

**If all checked**: ✅ Repository ingestion is working perfectly!

