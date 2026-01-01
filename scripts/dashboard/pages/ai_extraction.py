"""AI Extraction Tools Dashboard.

Provides UI for:
- OCR text extraction
- Web page scraping
- Entity normalization
"""

import streamlit as st

from scripts.dashboard.core.api import _api_post
from scripts.dashboard.core.state import require_auth


def render_ai_extraction_page():
    """Main entry point for AI Extraction Tools page."""
    require_auth()
    
    st.title("🤖 AI Extraction Tools")
    st.caption("Extract text from documents, scrape web pages, and normalize entity names")
    
    # Tab layout
    tab_ocr, tab_scraper, tab_normalizer = st.tabs([
        "📄 OCR",
        "🌐 Web Scraper",
        "🔬 Entity Normalizer"
    ])
    
    with tab_ocr:
        _render_ocr_tab()
    
    with tab_scraper:
        _render_scraper_tab()
    
    with tab_normalizer:
        _render_normalizer_tab()


def _render_ocr_tab():
    """OCR text extraction from images/PDFs."""
    st.subheader("📄 OCR Text Extraction")
    st.caption("Extract text from scanned documents and images using Tesseract OCR")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        uploaded_file = st.file_uploader(
            "Upload file (PDF, PNG, JPEG, TIFF)",
            type=["pdf", "png", "jpg", "jpeg", "tiff", "gif"],
            key="ocr_file"
        )
    
    with col2:
        language = st.selectbox(
            "Language",
            ["eng", "fra", "deu", "spa", "ita", "por", "nld", "rus", "chi_sim", "jpn"],
            format_func=lambda x: {
                "eng": "English",
                "fra": "French",
                "deu": "German",
                "spa": "Spanish",
                "ita": "Italian",
                "por": "Portuguese",
                "nld": "Dutch",
                "rus": "Russian",
                "chi_sim": "Chinese (Simplified)",
                "jpn": "Japanese",
            }.get(x, x),
            key="ocr_language"
        )
    
    if uploaded_file and st.button("🔍 Extract Text", key="ocr_extract", type="primary"):
        with st.spinner("Extracting text..."):
            # Send file to API
            import httpx
            import os
            
            API_BASE = os.environ.get("API_URL", "http://localhost:8000")
            
            try:
                with httpx.Client(timeout=60) as client:
                    files = {"file": (uploaded_file.name, uploaded_file.getvalue(), uploaded_file.type)}
                    response = client.post(
                        f"{API_BASE}/api/v1/extraction/ocr?language={language}",
                        files=files,
                    )
                    result = response.json()
                
                if result.get("success"):
                    st.success(f"✅ Extracted {result.get('word_count', 0)} words")
                    st.text_area(
                        "Extracted Text",
                        value=result.get("text", ""),
                        height=300,
                        key="ocr_result"
                    )
                else:
                    st.error(f"❌ {result.get('error', 'OCR failed')}")
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")


def _render_scraper_tab():
    """Web page content scraping."""
    st.subheader("🌐 Web Scraper")
    st.caption("Extract content from web pages for RAG ingestion")
    
    url = st.text_input(
        "URL",
        placeholder="https://example.com/article",
        key="scraper_url"
    )
    
    if url and st.button("🔍 Scrape", key="scraper_extract", type="primary"):
        with st.spinner("Scraping content..."):
            result = _api_post("/extraction/scrape", {"url": url})
            
            if result and result.get("success"):
                st.success(f"✅ Scraped {result.get('word_count', 0)} words")
                
                col1, col2 = st.columns(2)
                with col1:
                    st.write("**Title:**", result.get("title", "N/A"))
                with col2:
                    st.write("**Author:**", result.get("author", "N/A"))
                
                with st.expander("📖 Full Content", expanded=True):
                    st.write(result.get("content", ""))
            else:
                st.error(f"❌ {result.get('error', 'Scraping failed') if result else 'API error'}")


def _render_normalizer_tab():
    """Entity name normalization."""
    st.subheader("🔬 Entity Normalizer")
    st.caption("Normalize gene, compound, and disease names to standard identifiers")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        entity_type = st.selectbox(
            "Entity Type",
            ["gene", "compound", "disease"],
            key="normalize_type"
        )
    
    with col2:
        entity_name = st.text_input(
            "Entity Name",
            placeholder="e.g., BRCA1, aspirin, diabetes",
            key="normalize_name"
        )
    
    if entity_name and st.button("🔍 Lookup", key="normalize_lookup", type="primary"):
        with st.spinner("Looking up..."):
            result = _api_post("/extraction/normalize", {
                "entity_type": entity_type,
                "name": entity_name
            })
            
            if result and result.get("success"):
                st.success(f"✅ Found match for '{entity_name}'")
                
                data = result.get("result", {})
                
                if entity_type == "gene":
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.metric("Symbol", data.get("symbol", "N/A"))
                        st.metric("UniProt ID", data.get("uniprot_id", "N/A"))
                    with col_b:
                        st.metric("Entrez ID", data.get("entrez_id", "N/A"))
                        st.metric("Organism", data.get("organism", "N/A"))
                    st.write("**Full Name:**", data.get("name", "N/A"))
                
                elif entity_type == "compound":
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.metric("PubChem CID", data.get("cid", "N/A"))
                        st.metric("MW", f"{data.get('molecular_weight', 'N/A')}")
                    with col_b:
                        st.metric("Formula", data.get("molecular_formula", "N/A"))
                    st.write("**SMILES:**", data.get("smiles", "N/A"))
                    st.write("**InChI Key:**", data.get("inchi_key", "N/A"))
                
                elif entity_type == "disease":
                    st.metric("Name", data.get("name", "N/A"))
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.metric("MeSH ID", data.get("mesh_id", "N/A"))
                    with col_b:
                        st.metric("DOID", data.get("doid", "N/A"))
            else:
                st.warning(f"⚠️ {result.get('error', 'Not found') if result else 'API error'}")