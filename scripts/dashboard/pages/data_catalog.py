"""Data Catalog Dashboard - Browse entities, search columns, view lineage, manage glossary."""

from __future__ import annotations

import streamlit as st
from uuid import UUID
import requests
import pandas as pd

API_BASE = "http://localhost:8000/api/v1/catalog"


def main():
    """Entry point for page registry."""
    st.set_page_config(page_title="Data Catalog", page_icon="📚", layout="wide")
    st.title("📚 Data Catalog")

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
        st.header("Browse Data Entities")
        
        col1, col2 = st.columns([1, 2])
        
        with col1:
            # Filter controls
            st.subheader("Filters")
            
            entity_type = st.selectbox(
                "Entity Type",
                options=["All", "Table", "View", "Dataset", "Model"],
                index=0
            )
            
            program_filter = st.text_input("Program Filter", placeholder="e.g., compound-discovery")
            
            search_term = st.text_input("Search", placeholder="Search entity names...")
            
            if st.button("Refresh Catalog"):
                try:
                    response = requests.post(f"{API_BASE}/refresh")
                    if response.status_code == 200:
                        st.success("Catalog refreshed successfully!")
                        st.rerun()
                    else:
                        st.error(f"Failed to refresh: {response.status_code}")
                except Exception as e:
                    st.error(f"Error: {e}")
        
        with col2:
            # Entity list
            st.subheader("Entities")
            
            try:
                # Build query parameters
                params = {}
                if entity_type != "All":
                    params["entity_type"] = entity_type.lower()
                if program_filter:
                    params["program"] = program_filter
                if search_term:
                    params["search"] = search_term
                
                response = requests.get(f"{API_BASE}/entries", params=params)
                
                if response.status_code == 200:
                    entities = response.json()
                    
                    if entities:
                        for entity in entities[:20]:  # Limit to 20 for performance
                            with st.expander(f"{entity['display_name']} ({entity['entity_type']})"):
                                st.write(f"**Schema:** {entity.get('schema_name', 'N/A')}")
                                st.write(f"**Description:** {entity.get('description', 'No description')}")
                                st.write(f"**Last Updated:** {entity.get('last_updated', 'Unknown')}")
                                
                                if entity.get('tags'):
                                    st.write(f"**Tags:** {', '.join(entity['tags'])}")
                                
                                if st.button(f"View Details", key=f"details_{entity['id']}"):
                                    st.session_state.selected_entity = entity
                                    st.rerun()
                    else:
                        st.info("No entities found matching your criteria")
                else:
                    st.error(f"Failed to load entities: {response.status_code}")
                    
            except Exception as e:
                st.error(f"Error loading entities: {e}")
        
        # Selected entity details
        if st.session_state.selected_entity:
            st.subheader(f"Entity Details: {st.session_state.selected_entity['display_name']}")
            
            entity = st.session_state.selected_entity
            
            col1, col2 = st.columns(2)
            with col1:
                st.write(f"**Type:** {entity['entity_type']}")
                st.write(f"**Schema:** {entity.get('schema_name', 'N/A')}")
                st.write(f"**Created:** {entity.get('created_at', 'Unknown')}")
            
            with col2:
                st.write(f"**Program:** {entity.get('program', 'N/A')}")
                st.write(f"**Owner:** {entity.get('owner', 'Unknown')}")
                st.write(f"**Last Updated:** {entity.get('last_updated', 'Unknown')}")
            
            if entity.get('description'):
                st.write(f"**Description:** {entity['description']}")
            
            # Show columns
            try:
                response = requests.get(f"{API_BASE}/entries/{entity['id']}/columns")
                if response.status_code == 200:
                    columns = response.json()
                    if columns:
                        st.subheader("Columns")
                        df = pd.DataFrame(columns)
                        st.dataframe(df, use_container_width=True)
                    else:
                        st.info("No column metadata available")
            except Exception as e:
                st.error(f"Error loading columns: {e}")


    # ============================================================================
    # TAB 2: COLUMN SEARCH
    # ============================================================================

    with tab2:
        st.header("Column Search")
        
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.subheader("Search Filters")
            
            column_name = st.text_input("Column Name", placeholder="e.g., compound_id")
            data_type = st.selectbox(
                "Data Type",
                options=["All", "text", "integer", "float", "boolean", "date", "json"],
                index=0
            )
            has_description = st.checkbox("Has Description")
            
        with col2:
            st.subheader("Search Results")
            
            # Perform search when filters change
            try:
                params = {}
                if column_name:
                    params["name"] = column_name
                if data_type != "All":
                    params["data_type"] = data_type
                if has_description:
                    params["has_description"] = "true"
                
                if params:  # Only search if there are filters
                    response = requests.get(f"{API_BASE}/columns/search", params=params)
                    
                    if response.status_code == 200:
                        results = response.json()
                        
                        if results:
                            for result in results[:50]:  # Limit results
                                with st.expander(f"{result['entity_name']}.{result['column_name']}"):
                                    st.write(f"**Type:** {result['data_type']}")
                                    st.write(f"**Entity:** {result['entity_name']} ({result['entity_type']})")
                                    
                                    if result.get('description'):
                                        st.write(f"**Description:** {result['description']}")
                                    
                                    if result.get('sample_values'):
                                        st.write(f"**Sample Values:** {', '.join(map(str, result['sample_values'][:5]))}")
                                    
                                    # Update column metadata
                                    new_desc = st.text_area(
                                        "Update Description",
                                        value=result.get('description', ''),
                                        key=f"desc_{result['id']}"
                                    )
                                    
                                    if st.button(f"Update", key=f"update_{result['id']}"):
                                        try:
                                            update_data = {"description": new_desc}
                                            update_response = requests.put(
                                                f"{API_BASE}/columns/{result['id']}", 
                                                json=update_data
                                            )
                                            if update_response.status_code == 200:
                                                st.success("Description updated!")
                                                st.rerun()
                                            else:
                                                st.error("Failed to update description")
                                        except Exception as e:
                                            st.error(f"Error updating: {e}")
                        else:
                            st.info("No columns found matching your criteria")
                    else:
                        st.error(f"Search failed: {response.status_code}")
                else:
                    st.info("Enter search criteria to find columns")
                    
            except Exception as e:
                st.error(f"Error performing search: {e}")


    # ============================================================================
    # TAB 3: DATA LINEAGE
    # ============================================================================

    with tab3:
        st.header("Data Lineage")
        
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.subheader("Lineage Options")
            
            # Entity selector for lineage
            entity_name = st.text_input("Entity Name", placeholder="e.g., experiments")
            
            if st.button("Detect Lineage"):
                try:
                    response = requests.post(f"{API_BASE}/lineage/detect")
                    if response.status_code == 200:
                        st.success("Lineage detection completed!")
                        st.rerun()
                    else:
                        st.error("Failed to detect lineage")
                except Exception as e:
                    st.error(f"Error: {e}")
        
        with col2:
            st.subheader("Lineage Graph")
            
            if entity_name:
                try:
                    params = {"entity_name": entity_name}
                    response = requests.get(f"{API_BASE}/lineage", params=params)
                    
                    if response.status_code == 200:
                        lineage_data = response.json()
                        
                        if lineage_data.get("nodes"):
                            # Display lineage information
                            st.write("**Entities in Lineage:**")
                            for node in lineage_data["nodes"]:
                                st.write(f"- {node['name']} ({node['type']})")
                            
                            st.write("**Relationships:**")
                            for edge in lineage_data.get("edges", []):
                                st.write(f"- {edge['source']} → {edge['target']} ({edge['relationship_type']})")
                            
                            # Simple visualization using text
                            st.subheader("Lineage Visualization")
                            st.text("Upstream ← [Current Entity] → Downstream")
                            
                            # You could integrate graphviz here for better visualization
                            # For now, showing a simple text representation
                            
                        else:
                            st.info("No lineage data found for this entity")
                    else:
                        st.error(f"Failed to get lineage: {response.status_code}")
                        
                except Exception as e:
                    st.error(f"Error loading lineage: {e}")
            else:
                st.info("Enter an entity name to view its lineage")


    # ============================================================================
    # TAB 4: BUSINESS GLOSSARY
    # ============================================================================

    with tab4:
        st.header("Business Glossary")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("Add New Term")
            
            with st.form("add_term"):
                term_name = st.text_input("Term")
                definition = st.text_area("Definition")
                category = st.selectbox(
                    "Category",
                    options=["General", "Chemistry", "Biology", "Data Science", "Business"]
                )
                synonyms = st.text_input("Synonyms (comma-separated)")
                
                if st.form_submit_button("Add Term"):
                    try:
                        term_data = {
                            "term": term_name,
                            "definition": definition,
                            "category": category,
                            "synonyms": [s.strip() for s in synonyms.split(",") if s.strip()]
                        }
                        
                        response = requests.post(f"{API_BASE}/glossary", json=term_data)
                        
                        if response.status_code == 201:
                            st.success(f"Term '{term_name}' added successfully!")
                            st.rerun()
                        else:
                            st.error(f"Failed to add term: {response.status_code}")
                            
                    except Exception as e:
                        st.error(f"Error adding term: {e}")
        
        with col2:
            st.subheader("Search Terms")
            
            search_query = st.text_input("Search glossary", placeholder="Enter term or definition...")
            
            try:
                params = {"search": search_query} if search_query else {}
                response = requests.get(f"{API_BASE}/glossary", params=params)
                
                if response.status_code == 200:
                    terms = response.json()
                    
                    if terms:
                        for term in terms:
                            with st.expander(f"{term['term']} ({term['category']})"):
                                st.write(f"**Definition:** {term['definition']}")
                                
                                if term.get('synonyms'):
                                    st.write(f"**Synonyms:** {', '.join(term['synonyms'])}")
                                
                                st.write(f"**Created:** {term.get('created_at', 'Unknown')}")
                                
                                # Update/Delete options
                                col1, col2 = st.columns(2)
                                
                                with col1:
                                    if st.button(f"Edit", key=f"edit_{term['id']}"):
                                        # Could implement inline editing
                                        st.info("Edit functionality could be implemented here")
                                
                                with col2:
                                    if st.button(f"Delete", key=f"delete_{term['id']}"):
                                        try:
                                            delete_response = requests.delete(f"{API_BASE}/glossary/{term['id']}")
                                            if delete_response.status_code == 204:
                                                st.success("Term deleted!")
                                                st.rerun()
                                            else:
                                                st.error("Failed to delete term")
                                        except Exception as e:
                                            st.error(f"Error deleting: {e}")
                    else:
                        st.info("No terms found" if search_query else "No terms in glossary yet")
                else:
                    st.error(f"Failed to load glossary: {response.status_code}")
                    
            except Exception as e:
                st.error(f"Error loading glossary: {e}")
