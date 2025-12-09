# AI Agent Tools - UI Integration Plan

**Status**: ⚠️ **CRITICAL - MUST HAVE**  
**Last Updated**: 2025-12-04

## Overview

The web UI must include comprehensive AI agent tools to enable users to interact with the RAG system, cross-omics reasoning, and other AI-powered features directly through the browser interface.

## Current AI Capabilities (CLI-Only)

### ✅ Existing Features

1. **RAG Query System** (`scripts/rag_query.py`)
   - Natural language semantic search over Pinecone
   - Filters: source-type, disease, target, lipid, signature
   - AI-synthesized answers from retrieved chunks
   - Context display with citations

2. **Cross-Omics Reasoning** (`amprenta_rag/query/cross_omics_reasoning.py`)
   - AI-powered program summaries
   - AI-powered signature summaries
   - AI-powered feature summaries (genes, proteins, metabolites, lipids)
   - AI-powered dataset summaries

3. **Signature Similarity** (`rag_query.py --signature-score`)
   - Find matching signatures for datasets
   - Score and rank signature matches

4. **Metadata Classification** (`classify_literature_metadata.py`)
   - AI-powered metadata extraction
   - Disease, target, modality classification

## UI Integration Plan

### Phase 1: Streamlit Dashboard Integration (IMMEDIATE)

#### 1.1 RAG Query Interface Page

**New Dashboard Page: "AI Search" or "Query AI"**

Features:
- Large text input for natural language queries
- Filter sidebar:
  - Source type (Literature, Email, Experiment, Dataset, Signature)
  - Disease filter
  - Target filter
  - Lipid filter
  - Signature filter
- Results display:
  - Top matches with relevance scores
  - Expandable context chunks
  - AI-synthesized answer in prominent box
  - Source citations with links
- Export options (copy answer, download results)

**Implementation:**
```python
# In run_dashboard.py
if page == "AI Search":
    st.header("🤖 AI-Powered Search")
    
    query = st.text_input("Ask a question about your data...", "")
    
    col1, col2 = st.columns([3, 1])
    with col2:
        source_type = st.multiselect("Source Type", ["Literature", "Email", ...])
        disease = st.text_input("Disease filter", "")
        target = st.text_input("Target filter", "")
    
    if st.button("Search") and query:
        result = query_rag(
            query=query,
            source_types=source_type,
            disease=disease,
            target=target,
            top_k=10,
        )
        
        # Display AI answer
        st.info(result.answer)
        
        # Display matches
        for match in result.matches:
            with st.expander(f"[{match.score:.3f}] {match.title}"):
                st.write(match.snippet)
                st.caption(f"Source: {match.source}")
```

#### 1.2 Cross-Omics Reasoning Widget

**Add to existing entity detail pages:**

- **On Program Pages:**
  - "Ask AI about this Program" button
  - Calls `cross_omics_program_summary(program_id)`
  - Displays AI-generated summary

- **On Signature Pages:**
  - "Ask AI about this Signature" button
  - Calls `cross_omics_signature_summary(signature_id)`
  - Displays cross-omics analysis

- **On Feature Pages:**
  - "Ask AI about this Feature" button
  - Calls `cross_omics_feature_summary(feature_type, feature_name)`
  - Shows feature across omics types

- **On Dataset Pages:**
  - "Ask AI about this Dataset" button
  - Calls `cross_omics_dataset_summary(dataset_id)`
  - Shows dataset summary

**Implementation:**
```python
# In dataset detail view
if dataset:
    col1, col2 = st.columns([3, 1])
    with col2:
        if st.button("🤖 Ask AI about Dataset"):
            with st.spinner("Generating AI summary..."):
                summary = cross_omics_dataset_summary(
                    dataset_page_id=dataset.notion_page_id
                )
                st.info(summary)
```

#### 1.3 Signature Matching Interface

**On Dataset Pages:**

- "Find Matching Signatures" button
- Calls `signature_similarity_query(dataset_id)`
- Displays:
  - Ranked list of matching signatures
  - Similarity scores
  - Signature details
  - Visual similarity indicator

### Phase 2: Enhanced AI Features (SHORT-TERM)

#### 2.1 Interactive Chat Interface

**Chat-style UI:**
- Conversation history
- Follow-up questions
- Context-aware responses
- Export conversation

#### 2.2 Batch AI Operations

- Bulk classification interface
- Batch cross-omics summaries
- Bulk signature matching

#### 2.3 AI Insights Dashboard

- Automated insights generation
- Trend analysis
- Anomaly detection
- Recommendations

### Phase 3: Next.js Frontend Integration (LONG-TERM)

#### 3.1 Modern AI Chat Interface

- Real-time streaming responses
- Markdown rendering
- Code syntax highlighting
- Interactive charts in responses

#### 3.2 Advanced AI Features

- Multi-modal queries (text + visual)
- Voice input
- AI-powered data entry
- Automated report generation

## Technical Architecture

### Backend API Endpoints Needed

```python
# FastAPI endpoints for AI tools
POST /api/v1/ai/query
  - Query: string
  - Filters: object
  - Returns: RAGQueryResult

POST /api/v1/ai/cross-omics/program/{id}
  - Returns: AI summary

POST /api/v1/ai/cross-omics/signature/{id}
  - Returns: AI summary

POST /api/v1/ai/cross-omics/feature
  - Feature type: string
  - Feature name: string
  - Returns: AI summary

POST /api/v1/ai/cross-omics/dataset/{id}
  - Returns: AI summary

GET /api/v1/ai/signatures/match/{dataset_id}
  - Returns: List of matching signatures with scores
```

### Frontend Components

**Streamlit (Phase 1):**
- `ai_query_page.py` - RAG query interface
- `ai_widgets.py` - Reusable AI components
- Integration into existing dashboard pages

**Next.js (Phase 3):**
- `AIChatInterface.tsx` - Chat component
- `AISummaryWidget.tsx` - Summary display
- `SignatureMatcher.tsx` - Matching interface
- API client for AI endpoints

## User Experience Flow

### Example: User queries about a program

1. User navigates to "Programs" page
2. Clicks on a program (e.g., "ALS Research Program")
3. Sees program details
4. Clicks "🤖 Ask AI about this Program"
5. AI generates cross-omics summary:
   - Which omics types are represented
   - Key findings across omics
   - Related datasets and features
   - Signature associations
6. User can ask follow-up questions
7. User exports or saves summary

### Example: Natural language search

1. User navigates to "AI Search" page
2. Types: "What are the key lipid signatures associated with ALS?"
3. Applies filters (optional):
   - Source type: Literature, Dataset
   - Disease: ALS
4. Clicks "Search"
5. System:
   - Retrieves relevant chunks from Pinecone
   - Ranks and filters results
   - Generates AI-synthesized answer
6. User sees:
   - Clear answer at top
   - Ranked source matches below
   - Can expand context for each match
7. User can refine query or ask follow-up

## Implementation Priority

### Must Have (Phase 1 - Week 1-2):

1. ✅ RAG Query Interface page in Streamlit
2. ✅ Cross-omics reasoning widgets on entity pages
3. ✅ Signature matching interface on dataset pages

### Should Have (Phase 2 - Week 3-4):

4. Interactive chat interface
5. Batch AI operations
6. Improved visualizations for AI results

### Nice to Have (Phase 3 - Month 2+):

7. Next.js frontend with advanced AI features
8. Multi-modal queries
9. Automated insights

## Success Metrics

- **Adoption**: % of users using AI features
- **Query Volume**: Number of AI queries per day
- **User Satisfaction**: Feedback on AI answer quality
- **Time Saved**: Reduction in manual data exploration time

## Notes

- **API Rate Limits**: Consider OpenAI API costs and rate limits
- **Caching**: Cache AI responses where appropriate
- **Error Handling**: Graceful degradation if AI services unavailable
- **Privacy**: Ensure sensitive data is handled appropriately

---

**CRITICAL**: AI agent tools are essential for user adoption. The CLI-only approach limits accessibility. Prioritize UI integration immediately.

