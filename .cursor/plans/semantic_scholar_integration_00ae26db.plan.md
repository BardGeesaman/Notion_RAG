---
name: Semantic Scholar Integration
overview: Integrate Semantic Scholar and OpenAlex APIs to enable citation graph analysis, AI-powered paper relevance scoring, and author/institution metadata enrichment for the scientific paper management system.
todos:
  - id: ss-batch-1
    content: "Batch 1: Repository Clients (semantic_scholar.py, openalex.py)"
    status: in_progress
  - id: ss-batch-2
    content: "Batch 2: Database Schema (extend Paper, add Citation table)"
    status: pending
    dependencies:
      - ss-batch-1
  - id: ss-batch-3
    content: "Batch 3: API Endpoints (citations, references, enrich)"
    status: pending
    dependencies:
      - ss-batch-2
---

# Semantic Scholar /

OpenAlex Integration

## Overview

Add academic literature intelligence to the existing paper ingestion system by integrating:

- **Semantic Scholar API** - Citation graphs, AI summaries (TLDR), influential citations
- **OpenAlex API** - Open metadata, author affiliations, institutions, concepts

## Current State

Existing paper system (from Scientific Paper Ingestion feature):

- [papers.py](amprenta_rag/api/routers/papers.py) - CRUD endpoints
- [paper_search.py](scripts/dashboard/pages/paper_search.py) - UI page
- PubMed and bioRxiv ingestion working
- Author/affiliation linking implemented

## Integration Architecture

```javascript
┌─────────────────────────────────────────────────────────┐
│                    External APIs                         │
├─────────────────────────────────────────────────────────┤
│  Semantic Scholar API  │  OpenAlex API  │  Existing     │
│  - Paper search        │  - Works       │  - PubMed     │
│  - Citations           │  - Authors     │  - bioRxiv    │
│  - References          │  - Institutions│               │
│  - TLDR summaries      │  - Concepts    │               │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│              Repository Layer (New)                      │
├─────────────────────────────────────────────────────────┤
│  SemanticScholarRepository  │  OpenAlexRepository       │
│  - search_papers()          │  - get_work()             │
│  - get_citations()          │  - get_author()           │
│  - get_references()         │  - get_institution()      │
│  - get_paper_details()      │  - search_works()         │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│              Database Models (Extend)                    │
├─────────────────────────────────────────────────────────┤
│  Paper (extend)        │  Citation (new)                │
│  + semantic_scholar_id │  - citing_paper_id             │
│  + openalex_id         │  - cited_paper_id              │
│  + tldr_summary        │  - is_influential              │
│  + citation_count      │  - context                     │
│  + influential_count   │                                │
└─────────────────────────────────────────────────────────┘
```



## Batch Organization

### Batch 1: Repository Clients

**Files to create:**

- `amprenta_rag/ingestion/papers/semantic_scholar.py`
- `amprenta_rag/ingestion/papers/openalex.py`

**Features:**

- Semantic Scholar: search, paper details, citations, references, TLDR
- OpenAlex: work lookup, author details, institution details
- Rate limiting and error handling
- Response caching

### Batch 2: Database Schema

**Files to modify/create:**

- `amprenta_rag/models/content.py` - Extend Paper model
- `alembic/versions/xxx_add_citation_tracking.py` - Migration

**Schema additions:**

- Paper: semantic_scholar_id, openalex_id, tldr_summary, citation_count, influential_citations
- Citation table: citing_paper_id, cited_paper_id, is_influential, context

### Batch 3: API Endpoints

**Files to modify/create:**

- `amprenta_rag/api/routers/papers.py` - Add citation endpoints

**Endpoints:**