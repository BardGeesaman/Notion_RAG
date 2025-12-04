# Documentation

**Last Updated**: 2025-12-03

Welcome to the Amprenta RAG multi-omics platform documentation.

---

## Documentation Index

### Getting Started

- **[Architecture Overview](ARCHITECTURE.md)**: High-level system architecture and component overview
- **[Configuration Guide](CONFIGURATION.md)**: Environment setup and database configuration
- **[Usage Examples](USAGE_EXAMPLES.md)**: Practical code examples for common tasks

### Reference

- **[API Reference](API_REFERENCE.md)**: Comprehensive API documentation for all modules

---

## Quick Start

1. **Configure Environment**
   ```bash
   cp .env.example .env
   # Edit .env with your API keys and database IDs
   ```

2. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Ingest Data**
   ```python
   from amprenta_rag.ingestion.lipidomics_ingestion import ingest_lipidomics_file
   result = ingest_lipidomics_file("data/lipidomics.csv", "ST004396")
   ```

4. **Query RAG**
   ```python
   from amprenta_rag.query.rag_engine import query_rag
   results = query_rag("ceramide signature in ALS", top_k=10)
   ```

---

## Documentation Structure

```
docs/
├── README.md              # This file
├── ARCHITECTURE.md        # System architecture
├── API_REFERENCE.md        # API documentation
├── USAGE_EXAMPLES.md      # Code examples
└── CONFIGURATION.md       # Setup guide
```

---

## Key Concepts

### Multi-Omics Support

The platform supports four omics types:
- **Lipidomics**: Lipid species (Cer, SM, HexCer, etc.)
- **Metabolomics**: Metabolites (amino acids, nucleotides, etc.)
- **Proteomics**: Proteins (UniProt IDs, gene names)
- **Transcriptomics**: Genes (HGNC symbols, Ensembl IDs)

### Signature System

Signatures are collections of features with directions (↑/↓) and weights:
- Multi-omics signatures can include any combination of feature types
- Automatic scoring against datasets
- Direction-aware matching

### Feature Linking

Features are automatically linked to datasets:
- Creates/finds feature pages in Notion
- Maintains bidirectional relations
- Idempotent and non-destructive

### RAG Integration

All ingested data is embedded and indexed:
- Text representations built from metadata
- Chunked and embedded via OpenAI
- Stored in Pinecone for semantic search

---

## Support

For questions or issues:
1. Check the [Configuration Guide](CONFIGURATION.md) for setup issues
2. Review [Usage Examples](USAGE_EXAMPLES.md) for code patterns
3. Consult [API Reference](API_REFERENCE.md) for function details

---

## See Also

- Project README: [../README.md](../README.md)
- Test Suite: [../amprenta_rag/tests/](../amprenta_rag/tests/)

