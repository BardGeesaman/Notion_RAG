# Cross-Omics RAG Reasoning - Testing Guide

## Quick Test Options

### 1. Test Feature Summary (Easiest - No page IDs needed!)

Test with a known feature that's been ingested:

```bash
# Test with a lipid feature
python test_cross_omics.py --feature "lipid:Cer(d18:1/16:0)"

# Test with a metabolite feature
python test_cross_omics.py --feature "metabolite:Glutamate"

# Test with a gene feature
python test_cross_omics.py --feature "gene:TP53"

# Test with a protein feature
python test_cross_omics.py --feature "protein:P04637"
```

### 2. Test Using CLI Script

```bash
# Feature summary
python scripts/rag_query.py --cross-omics-feature "lipid:Cer(d18:1/16:0)" --top-k 10

# Program summary (need program page ID)
python scripts/rag_query.py --cross-omics-program <program_page_id> --top-k 10

# Signature summary (need signature page ID)
python scripts/rag_query.py --cross-omics-signature <signature_page_id> --top-k 10

# Dataset summary (need dataset page ID)
python scripts/rag_query.py --cross-omics-dataset <dataset_page_id> --top-k 10
```

## Finding Test Data

### Find Signature Page IDs

You can query Notion or use existing signature names:
- Look for "ALS-CSF-Core-6Ceramides" or other signatures in your Notion
- Or check the signature ingestion logs

### Find Dataset Page IDs

From previous ingestion reports:
- Check `LIPIDOMICS_INGESTION_IMPLEMENTATION.md`
- Check `ST004396_INGESTION_REPORT.md`
- Or list datasets: `python scripts/list_lipidomics_datasets.py`

### Find Program Page IDs

- Query Notion "Programs" database
- Check program ingestion logs

### Find Feature Page IDs

From `FEATURE_LINKING_COMPLETE_SUCCESS.md`:
- Metabolites: Glutamate, Glutamine, Serine, Lactate
- These have been linked and should have chunks

## Expected Behavior

### Successful Test Output Should Show:

1. **Feature Summary**:
   - Feature name and type
   - Number of linked datasets
   - Chunk counts by omics type
   - LLM-generated summary with:
     - Per-omics findings
     - Cross-omics convergence
     - Disease/context

2. **Program Summary**:
   - Program name
   - Linked experiments/datasets
   - Cross-omics findings across all modalities

3. **Signature Summary**:
   - Signature name and modalities
   - Matched datasets
   - Feature patterns

4. **Dataset Summary**:
   - Dataset name and omics type
   - Linked signatures
   - Cross-omics context

### Error Messages (Normal):

- "No sufficient multi-omics context found" - Expected if:
  - Feature/dataset hasn't been ingested
  - No chunks available in Pinecone
  - No linked datasets/experiments

## Troubleshooting

### ModuleNotFoundError

Activate your Python environment first:
```bash
# Try one of these:
source venv/bin/activate
source .venv/bin/activate
source myenv/bin/activate
# Or check what virtual environments exist:
ls -la | grep -E "venv|env"
```

### No Chunks Found

Ensure:
1. Datasets/features have been ingested into Pinecone
2. Feature linking has been completed
3. Chunks were successfully upserted

### OpenAI API Errors

Check:
1. `OPENAI_API_KEY` is set in `.env`
2. API key is valid
3. You have credits/quota available

## Test Checklist

- [ ] Feature summary works for at least one feature type
- [ ] CLI script shows help text with cross-omics options
- [ ] Error handling works (invalid IDs, missing data)
- [ ] Logging appears correctly
- [ ] LLM generates coherent summaries

## Next Steps After Testing

1. **If tests pass**: System is ready for production use!
2. **If errors**: Check logs, fix configuration, verify data ingestion
3. **If slow**: Consider optimizing chunk retrieval or caching

---

**Ready to test? Start with a feature summary - it's the easiest!**

```bash
python test_cross_omics.py --feature "lipid:Cer(d18:1/16:0)"
```

