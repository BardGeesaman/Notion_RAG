# Amprenta RAG Architecture (v0.1)

Last updated: [YYYY-MM-DD]
Owner: [Name]

## 1. Goals

- Provide Amprenta teammates with a unified assistant that can answer questions using:
  - Notion content (pages & databases)
  - External documents ([e.g. Google Drive, PDFs])
  - Structured data ([e.g. Postgres, analytics])
- Priorities:
  - High factual accuracy
  - Clear citations back to sources (URLs, Notion links, etc.)
  - Easy to add new data sources

## 2. High-level architecture (text diagram)

1. **Data Sources**
   - Notion:
     - [List key DBs and page trees, e.g. Projects, Docs, RAG System]
   - Files:
     - [Drive folders / S3 buckets / etc.]
   - Structured data:
     - [Database tables, APIs]

2. **Ingestion layer**
   - Connectors:
     - Notion API connector: pulls pages/DB entries via [your method/tool].
     - File connector: scans configured folders and extracts text.
     - Structured data connector: queries [DB / APIs].
   - Sync:
     - Full sync: [describe when/how]
     - Incremental sync: [describe change detection strategy]

3. **Preprocessing & chunking**
   - Steps:
     - Clean text (remove boilerplate, navigation, etc.).
     - Normalize formatting (headings, bullet lists, tables).
     - Chunk strategy:
       - For Notion: [by heading / block / combination]
       - For docs: [fixed token size / semantic paragraphs]
     - Attach metadata:
       - `source_system` (notion / drive / db)
       - `source_id` (Notion page ID, file path, row ID)
       - `source_url` or `notion_url`
       - `doc_type`, `owner`, `sensitivity`, timestamps

4. **Embeddings & storage**
   - Embedding model:
     - [e.g. `text-embedding-3-large` (OpenAI)]
   - Vector store:
     - [Name & tech, e.g. Postgres + pgvector, Pinecone, Weaviate, etc.]
   - Index layout:
     - One unified index vs per-source indexes:
       - [Describe your choice]
     - Additional dense/keyword index if used (e.g. BM25, hybrid search).

5. **Retrieval pipeline**
   - Inputs:
     - User query
     - Conversation context (recent turns)
     - Optional filters (source, date range, sensitivity)
   - Steps:
     1) Rewrite/expand query if needed (query understanding step).
     2) Vector retrieval: top-k [N] documents from vector store.
     3) Optional re-ranking:
        - [Model or heuristic used]
     4) Apply filters:
        - Permission filters
        - Source filters
        - Date range
   - Outputs:
     - Ranked list of chunks + metadata, with source URLs.

6. **Answer generation**
   - LLM: [which model(s) you use]
   - System prompt outline:
     - Role: “You are an assistant for Amprenta…”
     - Rules:
       - Only use retrieved context as ground truth.
       - Show citations with source titles/URLs.
       - Say “I don’t know” when context is missing.
   - Templates:
     - Q&A
     - Summarization of a doc or set of docs
     - “Where is X stored?” / navigation help

7. **Evaluation & monitoring**
   - Logging:
     - Store anonymized queries, retrieved docs, and responses.
   - Evaluation:
     - Golden question set in `RAG Evaluations` Notion DB:
       - `Question`
       - `Expected answer (short)`
       - `Must-use sources`
       - `Eval owner`
   - Metrics:
     - Retrieval quality (did we hit the right docs?)
     - Answer correctness (manual or LLM-graded)
     - Latency & cost

## 3. MVP scope

- Data sources in MVP:
  - [List the 1–3 Notion DBs + any other critical sources]
- Out of scope for MVP:
  - [e.g. complex row-level access control, real-time sync, etc.]

## 4. Open questions

- [Which exact vector DB/service will we use first?]
- [Do we need per-client isolation in the index?]
- [How will we deploy this (CLI tool, API, internal app, etc.)?]