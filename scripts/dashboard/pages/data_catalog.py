"""Data Catalog Dashboard - Browse entities, search columns, view lineage, manage glossary."""

from __future__ import annotations

import streamlit as st
from uuid import UUID
import requests
import pandas as pd

st.set_page_config(page_title="Data Catalog", page_icon="📚", layout="wide")
st.title("📚 Data Catalog")

API_BASE = "http://localhost:8000/api/v1/catalog"

# Initialize session state
if "selected_entity" not in st.session_state:
    st.session_state.selected_entity = None

# ============================================================================
# TAB LAYOUT
# ============================================================================

tab1, tab2, tab3, tab4 = st.tabs([
    "📊 Browse Entities",
    "🔍 Column Search", 
    "🌳 Data Lineage",
    "📖 Business Glossary"
])


# ============================================================================
# TAB 1: BROWSE ENTITIES
# ============================================================================

with tab1:
    st.subheader("Entity Types")
    
    # Filters
    col1, col2, col3 = st.columns([2, 2, 1])
    with col1:
        category_filter = st.selectbox(
            "Category",
            ["All", "Core", "Chemistry", "Omics", "Admin", "Other"],
            index=0
        )
    with col2:
        search_filter = st.text_input("Search entities", placeholder="e.g., Compound")
    with col3:
        if st.button("🔄 Refresh Catalog", help="Discover new models"):
            # Call refresh endpoint
            try:
                resp = requests.post(f"{API_BASE}/refresh")
                if resp.ok:
                    result = resp.json()
                    st.success(f"Refreshed! {result.get('entries_updated', 0)} entries updated, {result.get('lineage_edges_created', 0)} edges created")
                    st.rerun()
                else:
                    st.error(f"Refresh failed: {resp.status_code}")
            except Exception as e:
                st.error(f"Refresh failed: {e}")
    
    # Fetch entries
    params = {}
    if category_filter != "All":
        params["category"] = category_filter
    if search_filter:
        params["search"] = search_filter
    
    try:
        resp = requests.get(f"{API_BASE}/entries", params=params, timeout=10)
        if resp.ok:
            entries = resp.json()
        else:
            st.error(f"API error: {resp.status_code} - {resp.text[:200]}")
            entries = []
    except requests.exceptions.ConnectionError:
        st.warning("⚠️ Cannot connect to API server. Make sure FastAPI is running on port 8000.")
        entries = []
    except Exception as e:
        st.error(f"Failed to fetch entries: {e}")
        entries = []
    
    # Display as cards grouped by category
    if entries:
        # Group by category
        by_category = {}
        for e in entries:
            cat = e.get("category", "Other")
            by_category.setdefault(cat, []).append(e)
        
        for cat, items in sorted(by_category.items()):
            st.markdown(f"### {cat}")
            cols = st.columns(3)
            for i, entry in enumerate(items):
                with cols[i % 3]:
                    with st.container(border=True):
                        st.markdown(f"**{entry['display_name']}**")
                        st.caption(f"`{entry['table_name']}`")
                        if entry.get('description'):
                            desc = entry.get('description', '')
                            st.write(desc[:100] + "..." if len(desc) > 100 else desc)
                        
                        col_a, col_b = st.columns(2)
                        row_count = entry.get('row_count')
                        col_a.metric("Rows", f"{row_count:,}" if row_count else "N/A")
                        # Calculate column count from API response
                        col_b.metric("Columns", "?")  # Will be updated when viewing details
                        
                        if st.button("View Columns", key=f"view_{entry['entity_type']}"):
                            st.session_state.selected_entity = entry['entity_type']
                            st.rerun()
        
        # Show columns for selected entity
        if st.session_state.get('selected_entity'):
            st.markdown("---")
            entity_type = st.session_state.selected_entity
            st.subheader(f"Columns: {entity_type}")
            
            try:
                resp = requests.get(f"{API_BASE}/entries/{entity_type}")
                if resp.ok:
                    detail = resp.json()
                    columns = detail.get('columns', [])
                    
                    # Display as table
                    if columns:
                        df = pd.DataFrame([{
                            "Column": c['column_name'],
                            "Type": c['data_type'],
                            "Nullable": "✓" if c['is_nullable'] else "✗",
                            "PK": "🔑" if c['is_primary_key'] else "",
                            "FK": f"→ {c['foreign_key_target']}" if c['is_foreign_key'] else "",
                            "Description": c.get('description', '')[:50] or "-"
                        } for c in columns])
                        
                        st.dataframe(df, use_container_width=True, hide_index=True)
                        st.caption(f"Total columns: {len(columns)}")
                    else:
                        st.info("No columns found")
                    
                    if st.button("Close"):
                        st.session_state.selected_entity = None
                        st.rerun()
                else:
                    st.error(f"Failed to load columns: {resp.status_code}")
            except Exception as e:
                st.error(f"Failed to load columns: {e}")
    else:
        st.info("No entities found. Click 'Refresh Catalog' to discover models.")


# ============================================================================
# TAB 2: COLUMN SEARCH
# ============================================================================

with tab2:
    st.subheader("Search Columns Across All Entities")
    
    search_query = st.text_input(
        "Search",
        placeholder="e.g., patient_id, ic50, created_at",
        help="Search by column name or description"
    )
    
    if search_query and len(search_query) >= 2:
        try:
            resp = requests.get(f"{API_BASE}/columns/search", params={"q": search_query, "limit": 100})
            results = resp.json() if resp.ok else []
        except Exception as e:
            st.error(f"Search failed: {e}")
            results = []
        
        if results:
            st.success(f"Found {len(results)} columns matching '{search_query}'")
            
            df = pd.DataFrame([{
                "Entity": r['entity_type'],
                "Column": r['column_name'],
                "Type": r['data_type'],
                "Description": r.get('description', '-')[:80]
            } for r in results])
            
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            # Quick stats
            st.markdown("### Usage Summary")
            entity_counts = df['Entity'].value_counts()
            st.bar_chart(entity_counts)
        else:
            st.warning(f"No columns found matching '{search_query}'")
    elif search_query:
        st.info("Enter at least 2 characters to search")
    else:
        st.info("Enter a search term to find columns across all entities")


# ============================================================================
# TAB 3: DATA LINEAGE
# ============================================================================

with tab3:
    st.subheader("Data Lineage Graph")
    
    col1, col2, col3 = st.columns([2, 2, 1])
    with col1:
        entity_type = st.text_input("Entity Type", placeholder="e.g., Dataset")
    with col2:
        entity_id = st.text_input("Entity ID", placeholder="UUID")
    with col3:
        depth = st.slider("Depth", 1, 5, 3)
    
    direction = st.radio("Direction", ["both", "upstream", "downstream"], horizontal=True)
    
    if st.button("Generate Lineage Graph") and entity_type and entity_id:
        try:
            resp = requests.get(
                f"{API_BASE}/lineage/{entity_type}/{entity_id}",
                params={"depth": depth, "direction": direction}
            )
            if resp.ok:
                graph = resp.json()
                nodes = graph.get('nodes', [])
                edges = graph.get('edges', [])
                
                if nodes:
                    st.success(f"Found {len(nodes)} nodes, {len(edges)} edges")
                    
                    # Simple visualization using graphviz
                    try:
                        import graphviz
                        dot = graphviz.Digraph()
                        dot.attr(rankdir='LR')
                        
                        # Color coding by category
                        colors = {
                            'Core': '#4CAF50',
                            'Chemistry': '#2196F3', 
                            'Omics': '#9C27B0',
                            'Admin': '#FF9800',
                            'Other': '#607D8B'
                        }
                        
                        for node in nodes:
                            node_data = node.get('data', {})
                            node_id = node_data.get('id', '')
                            label = node_data.get('label', node_id)
                            is_center = node_data.get('is_center', False)
                            
                            # Highlight center node
                            if is_center:
                                dot.node(node_id, label, style='filled', fillcolor='#FFD700', fontcolor='black')
                            else:
                                dot.node(node_id, label, style='filled', fillcolor='#E3F2FD', fontcolor='black')
                        
                        for edge in edges:
                            edge_data = edge.get('data', {})
                            source = edge_data.get('source', '')
                            target = edge_data.get('target', '')
                            rel = edge_data.get('relationship', '')
                            dot.edge(source, target, label=rel)
                        
                        st.graphviz_chart(dot)
                        
                        # Node list
                        with st.expander("Node Details"):
                            for node in nodes:
                                data = node.get('data', {})
                                st.write(f"- **{data.get('label')}** (Type: {data.get('type', 'N/A')})")
                    
                    except ImportError:
                        st.warning("Graphviz not available. Showing node/edge data:")
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            st.markdown("**Nodes:**")
                            for node in nodes:
                                data = node.get('data', {})
                                st.write(f"- {data.get('label')} ({data.get('type')})")
                        
                        with col2:
                            st.markdown("**Edges:**")
                            for edge in edges:
                                data = edge.get('data', {})
                                st.write(f"- {data.get('source')} → {data.get('target')} ({data.get('relationship')})")
                else:
                    st.info("No lineage data found for this entity")
            else:
                st.error(f"API error: {resp.status_code} - {resp.text}")
        except Exception as e:
            st.error(f"Failed to fetch lineage: {e}")
    else:
        st.info("Enter entity type and ID to view lineage graph")
        
        # Help text
        with st.expander("ℹ️ How to use"):
            st.markdown("""
            1. Enter the **entity type** (e.g., `Dataset`, `Experiment`, `Compound`)
            2. Enter the **entity ID** (UUID from the database)
            3. Adjust **depth** (how many hops to traverse)
            4. Select **direction**:
               - `upstream`: Where did this data come from?
               - `downstream`: Where does this data flow to?
               - `both`: Full lineage in both directions
            """)


# ============================================================================
# TAB 4: BUSINESS GLOSSARY
# ============================================================================

with tab4:
    st.subheader("Business Glossary")
    
    # Search and filter
    col1, col2 = st.columns([3, 1])
    with col1:
        glossary_search = st.text_input("Search terms", placeholder="e.g., IC50")
    with col2:
        glossary_category = st.selectbox(
            "Category",
            ["All", "Chemistry", "Biology", "Statistics", "Other"]
        )
    
    # Fetch terms
    params = {}
    if glossary_search:
        params["search"] = glossary_search
    if glossary_category != "All":
        params["category"] = glossary_category
    
    try:
        resp = requests.get(f"{API_BASE}/glossary", params=params)
        terms = resp.json() if resp.ok else []
    except Exception as e:
        st.error(f"Failed to fetch glossary: {e}")
        terms = []
    
    # Add new term form
    with st.expander("➕ Add New Term"):
        with st.form("add_term"):
            new_term = st.text_input("Term *")
            new_definition = st.text_area("Definition *")
            new_category = st.selectbox("Category", ["Chemistry", "Biology", "Statistics", "Other"])
            new_synonyms = st.text_input("Synonyms (comma-separated)")
            
            if st.form_submit_button("Create Term"):
                if new_term and new_definition:
                    payload = {
                        "term": new_term,
                        "definition": new_definition,
                        "category": new_category,
                        "synonyms": [s.strip() for s in new_synonyms.split(",") if s.strip()] if new_synonyms else []
                    }
                    try:
                        resp = requests.post(f"{API_BASE}/glossary", json=payload)
                        if resp.ok:
                            st.success(f"Created term: {new_term}")
                            st.rerun()
                        else:
                            error_detail = resp.json().get('detail', 'Unknown error') if resp.headers.get('content-type', '').startswith('application/json') else resp.text
                            st.error(f"Failed: {error_detail}")
                    except Exception as e:
                        st.error(f"Error: {e}")
                else:
                    st.warning("Term and definition are required")
    
    # Display terms
    if terms:
        st.markdown(f"### Found {len(terms)} terms")
        
        for term in terms:
            with st.container(border=True):
                col1, col2 = st.columns([4, 1])
                with col1:
                    st.markdown(f"### {term['term']}")
                    st.write(term['definition'])
                    
                    if term.get('synonyms'):
                        st.caption(f"**Synonyms:** {', '.join(term['synonyms'])}")
                    if term.get('category'):
                        st.caption(f"**Category:** {term['category']}")
                    if term.get('created_at'):
                        st.caption(f"**Created:** {term['created_at'][:10]}")
                
                with col2:
                    if st.button("🗑️", key=f"del_{term['id']}", help="Delete term"):
                        try:
                            resp = requests.delete(f"{API_BASE}/glossary/{term['id']}")
                            if resp.ok:
                                st.success("Deleted")
                                st.rerun()
                            else:
                                st.error(f"Delete failed: {resp.status_code}")
                        except Exception as e:
                            st.error(f"Delete failed: {e}")
    else:
        if glossary_search or glossary_category != "All":
            st.info("No terms found matching your filters")
        else:
            st.info("No glossary terms found. Add your first term above!")


# ============================================================================
# REGISTER IN PAGE_REGISTRY
# ============================================================================

def main():
    """Entry point for page registry."""
    return  # Page renders on import
