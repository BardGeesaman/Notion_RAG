# Cross-Omics RAG Reasoning

Comprehensive guide to the cross-omics reasoning capabilities in the Amprenta RAG System.

## Overview

The cross-omics reasoning system generates intelligent, multi-omics summaries by:
- Retrieving evidence from multiple omics types
- Analyzing cross-omics patterns and convergence
- Synthesizing insights using LLM-powered reasoning
- Providing disease, model system, and matrix context

## Features

### ✅ Implemented Capabilities

1. **Program-Level Summaries**
   - Aggregates evidence from all linked experiments and datasets
   - Analyzes multi-omics patterns across the program
   - Provides disease and model system context

2. **Signature-Level Summaries**
   - Analyzes how signatures appear across datasets
   - Identifies cross-omics support for signature components
   - Contextualizes signatures within disease/model systems

3. **Feature-Level Summaries**
   - Summarizes evidence for specific features (genes, proteins, metabolites, lipids)
   - Identifies cross-omics patterns involving the feature
   - Provides disease and matrix context

4. **Dataset-Level Summaries**
   - Generates comprehensive summaries of datasets
   - Highlights multi-omics evidence
   - Contextualizes within disease and model systems

### Context Awareness

The system automatically extracts and incorporates:

- **Disease Context**: Disease names, conditions being studied
- **Model System Context**: In vitro, in vivo, patient-derived samples
- **Matrix Context**: CSF, plasma, serum, tissue, etc.
- **Comparative Context**: Disease vs control, treatment responses

## Usage

### Program Summary

```bash
python scripts/rag_query.py \
  --cross-omics-program PROGRAM_PAGE_ID \
  --top-k-per-omics 20
```

```python
from amprenta_rag.query.cross_omics_reasoning import cross_omics_program_summary

summary = cross_omics_program_summary(
    program_page_id="your-program-id",
    top_k_per_omics=20,
)
print(summary)
```

### Signature Summary

```bash
python scripts/rag_query.py \
  --cross-omics-signature SIGNATURE_PAGE_ID \
  --top-k-datasets 20 \
  --top-k-chunks 100
```

```python
from amprenta_rag.query.cross_omics_reasoning import cross_omics_signature_summary

summary = cross_omics_signature_summary(
    signature_page_id="your-signature-id",
    top_k_datasets=20,
    top_k_chunks=100,
)
```

### Feature Summary

```bash
python scripts/rag_query.py \
  --cross-omics-feature "Cer(d18:1/16:0)" lipid \
  --top-k-datasets 20
```

```python
from amprenta_rag.query.cross_omics_reasoning import cross_omics_feature_summary

summary = cross_omics_feature_summary(
    feature_name="Cer(d18:1/16:0)",
    feature_type="lipid",
    top_k_datasets=20,
)
```

### Dataset Summary

```bash
python scripts/rag_query.py \
  --cross-omics-dataset DATASET_PAGE_ID \
  --top-k-chunks 100
```

```python
from amprenta_rag.query.cross_omics_reasoning import cross_omics_dataset_summary

summary = cross_omics_dataset_summary(
    dataset_page_id="your-dataset-id",
    top_k_chunks=100,
)
```

## Architecture

### Modular Structure

```
amprenta_rag/query/cross_omics/
├── __init__.py              # Public API exports
├── program_summary.py       # Program-level summaries
├── signature_summary.py     # Signature-level summaries
├── feature_summary.py       # Feature-level summaries
├── dataset_summary.py       # Dataset-level summaries
├── context_extraction.py    # Context extraction utilities
├── helpers.py              # Shared helper functions
├── prompt_templates.py     # LLM prompt templates
└── synthesis.py            # LLM synthesis engine
```

### Data Flow

```
1. Identify Objects (Programs/Signatures/Features/Datasets)
   ↓
2. Extract Related Objects (Experiments, Datasets, Features)
   ↓
3. Extract Context (Disease, Model System, Matrix)
   ↓
4. Retrieve Chunks from Pinecone (with metadata filters)
   ↓
5. Group by Omics Type
   ↓
6. Build Enhanced Prompt (with context)
   ↓
7. Synthesize Summary via LLM
   ↓
8. Return Markdown Summary
```

## Context Extraction

### Disease Context

The system extracts disease information from:
- Program properties
- Experiment properties
- Dataset metadata
- Literature metadata

### Model System Context

Identifies model systems from:
- Experiment "Model Systems" property
- Dataset metadata
- Sample type information

### Matrix Context

Extracts matrix/tissue information from:
- Dataset "Matrix" property
- Sample type metadata
- Experimental context

### Comparative Context

Identifies comparative analyses from:
- Experiment designs (control vs treatment)
- Dataset comparisons
- Temporal patterns

## Prompt Templates

The system uses sophisticated prompt templates that:

1. **Provide Structure**: Guide LLM to generate well-organized summaries
2. **Include Context**: Automatically incorporate disease, model, matrix context
3. **Request Convergence Analysis**: Ask for cross-omics pattern identification
4. **Request Divergence Analysis**: Identify conflicting or modality-specific signals

### Example Prompt Structure

```
High-level context:
- Disease context: [extracted from metadata]
- Model system: [extracted from metadata]
- Matrix: [extracted from metadata]

Per-omics findings:
- Lipidomics: [key changes]
- Metabolomics: [key changes]
- Proteomics: [key changes]
- Transcriptomics: [key changes]

Cross-omics convergence:
- Features consistently changing across omics
- Convergent pathway signals

Cross-omics divergence:
- Conflicting signals
- Modality-specific patterns
```

## Performance

### Optimization Features

- **Chunk Filtering**: Uses metadata filters to retrieve relevant chunks only
- **Top-K Limiting**: Configurable limits per omics type
- **Context Caching**: Caches context extraction results
- **Parallel Retrieval**: Retrieves chunks in parallel when possible

### Typical Performance

- Program Summary: 10-30 seconds
- Signature Summary: 5-15 seconds
- Feature Summary: 5-10 seconds
- Dataset Summary: 3-8 seconds

## Best Practices

1. **Use Appropriate Top-K Values**
   - Start with default values (20 per omics type)
   - Increase for comprehensive summaries
   - Decrease for faster summaries

2. **Verify Object Relationships**
   - Ensure Programs have linked Experiments/Datasets
   - Verify Signatures have linked components
   - Check Features are linked to datasets

3. **Review Context Extraction**
   - Verify disease/model/matrix context is correctly extracted
   - Check Notion database properties are properly configured

4. **Iterate on Prompts**
   - Review generated summaries
   - Adjust prompt templates if needed
   - Fine-tune context extraction

## Limitations

1. **LLM Dependency**: Requires OpenAI API access
2. **Context Quality**: Depends on metadata quality in Notion
3. **Chunk Retrieval**: Limited by Pinecone index completeness
4. **Token Limits**: Very large contexts may be truncated

## Future Enhancements

Potential enhancements (TIER 1.3):

1. **Deeper Disease Context**
   - Disease progression stage awareness
   - Disease-specific pathway knowledge
   - Disease-model associations

2. **Enhanced Comparative Analysis**
   - Temporal pattern recognition
   - Treatment response analysis
   - Cross-disease comparisons

3. **Better Context Extraction**
   - More sophisticated metadata parsing
   - Relationship graph traversal
   - Context validation

4. **Customizable Templates**
   - User-defined prompt templates
   - Custom summary formats
   - Domain-specific templates

## Examples

See [Usage Examples](USAGE_EXAMPLES.md) for detailed examples of cross-omics reasoning.

## Troubleshooting

### No Context Extracted

**Issue**: Summary lacks disease/model/matrix context

**Solution**:
1. Verify Notion database properties are populated
2. Check property names match expected schema
3. Review context extraction logs

### Empty Summary

**Issue**: Summary is empty or minimal

**Solution**:
1. Verify objects have linked datasets/experiments
2. Check Pinecone index has chunks
3. Increase `top_k` parameters
4. Verify metadata filters aren't too restrictive

### Slow Performance

**Issue**: Summary generation takes too long

**Solution**:
1. Reduce `top_k` parameters
2. Check network connectivity
3. Verify API rate limits aren't being hit
4. Consider caching results

See [Troubleshooting Guide](TROUBLESHOOTING.md) for more details.

## API Reference

See [API Reference](API_REFERENCE.md) for detailed function signatures and parameters.

