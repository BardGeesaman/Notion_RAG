"""AI Expert Training Console - Admin interface for expert management."""

from __future__ import annotations

import streamlit as st
import requests
import pandas as pd
from datetime import datetime
from uuid import UUID

API_BASE = "http://localhost:8000/api/v1/experts"


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_experts():
    """Fetch available experts."""
    try:
        resp = requests.get(f"{API_BASE}/")
        return resp.json() if resp.ok else []
    except:
        return []


def get_training_examples(expert_id=None):
    """Fetch training examples."""
    try:
        params = {"expert_id": expert_id} if expert_id else {}
        resp = requests.get(f"{API_BASE}/training/examples", params=params)
        return resp.json() if resp.ok else []
    except:
        return []


def create_training_example(data):
    """Create new training example."""
    try:
        resp = requests.post(f"{API_BASE}/training/examples", json=data)
        return resp.json() if resp.ok else None
    except:
        return None


def get_knowledge_docs(expert_id=None):
    """Fetch knowledge documents."""
    try:
        params = {"expert_id": expert_id} if expert_id else {}
        resp = requests.get(f"{API_BASE}/knowledge", params=params)
        return resp.json() if resp.ok else []
    except:
        return []


def upload_knowledge_doc(data):
    """Upload knowledge document."""
    try:
        resp = requests.post(f"{API_BASE}/knowledge", json=data)
        return resp.json() if resp.ok else None
    except:
        return None


def get_feedback_summary(expert_id=None):
    """Get feedback summary."""
    try:
        params = {"expert_id": expert_id} if expert_id else {}
        resp = requests.get(f"{API_BASE}/feedback/summary", params=params)
        return resp.json() if resp.ok else {}
    except:
        return {}


def update_expert_config(expert_id, config):
    """Update expert configuration."""
    try:
        resp = requests.put(f"{API_BASE}/{expert_id}/config", json=config)
        return resp.ok
    except:
        return False


def export_training_data(expert_id=None, format_type="jsonl"):
    """Export training data."""
    try:
        params = {
            "expert_id": expert_id,
            "format": format_type
        }
        resp = requests.get(f"{API_BASE}/training/export", params=params)
        return resp.text if resp.ok else None
    except:
        return None


def main():
    """Entry point for page registry."""
    st.set_page_config(page_title="AI Expert Training", page_icon="🎓", layout="wide")

    # Initialize session state
    if "selected_expert_id" not in st.session_state:
        st.session_state.selected_expert_id = None
    if "selected_expert" not in st.session_state:
        st.session_state.selected_expert = None

    st.title("🎓 AI Expert Training Console")
    st.markdown("Admin interface for managing expert training data, knowledge base, and performance.")

    # ============================================================================
    # EXPERT SELECTOR
    # ============================================================================

    experts = get_experts()

    if not experts:
        st.warning("⚠️ No experts available. Please check the API connection.")
        # Continue with page layout even without experts
        experts = []

    # Expert selection
    col1, col2 = st.columns([1, 3])

    with col1:
        if experts:
            expert_options = {f"{exp['name']} ({', '.join(exp['specializations'][:2])})": exp for exp in experts}
            options = ["All Experts"] + list(expert_options.keys())
        else:
            expert_options = {}
            options = ["No Experts Available"]
            
        selected_expert_name = st.selectbox(
            "Select Expert",
            options=options,
            index=0
        )
        
        if selected_expert_name != "All Experts" and selected_expert_name != "No Experts Available":
            st.session_state.selected_expert = expert_options[selected_expert_name]
            st.session_state.selected_expert_id = st.session_state.selected_expert['id']
        else:
            st.session_state.selected_expert = None
            st.session_state.selected_expert_id = None

    with col2:
        if st.session_state.selected_expert:
            expert = st.session_state.selected_expert
            st.info(f"**{expert['name']}** - {', '.join(expert['specializations'])}")
            st.write(expert.get('description', 'No description available'))
        else:
            st.info("View training data across all experts")

    # ============================================================================
    # TAB LAYOUT
    # ============================================================================

    tab1, tab2, tab3, tab4 = st.tabs([
        "📚 Training Examples",
        "📖 Knowledge Base",
        "📊 Feedback Review",
        "⚙️ Expert Configuration"
    ])

    # ============================================================================
    # TAB 1: TRAINING EXAMPLES
    # ============================================================================

    with tab1:
        st.header("Training Examples Management")
        
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.subheader("Add Training Example")
            
            with st.form("add_example"):
                input_text = st.text_area("User Input", placeholder="Enter user question or prompt...")
                output_text = st.text_area("Expected Output", placeholder="Enter expert response...")
                
                if experts:
                    expert_id_for_example = st.selectbox(
                        "Target Expert",
                        options=[exp['id'] for exp in experts],
                        format_func=lambda x: next(exp['name'] for exp in experts if exp['id'] == x),
                        index=0
                    )
                else:
                    st.info("No experts available to create training examples for.")
                    expert_id_for_example = None
                
                category = st.selectbox(
                    "Category",
                    options=["General", "Chemistry", "Biology", "Data Analysis", "Safety", "Regulatory"]
                )
                
                difficulty = st.selectbox("Difficulty", options=["Easy", "Medium", "Hard"])
                
                tags = st.text_input("Tags (comma-separated)", placeholder="molecular-design, toxicity, etc.")
                
                if st.form_submit_button("Add Example"):
                    if input_text and output_text:
                        example_data = {
                            "expert_id": expert_id_for_example,
                            "input_text": input_text,
                            "output_text": output_text,
                            "category": category,
                            "difficulty": difficulty.lower(),
                            "tags": [tag.strip() for tag in tags.split(",") if tag.strip()],
                            "is_approved": True  # Auto-approve admin-created examples
                        }
                        
                        result = create_training_example(example_data)
                        if result:
                            st.success("Training example added successfully!")
                            st.rerun()
                        else:
                            st.error("Failed to add training example")
                    else:
                        st.error("Please provide both input and output text")
        
        with col2:
            st.subheader("Existing Training Examples")
            
            examples = get_training_examples(st.session_state.selected_expert_id)
            
            if examples:
                # Filter and search
                search_term = st.text_input("Search examples", placeholder="Search by input text...")
                category_filter = st.selectbox(
                    "Filter by Category",
                    options=["All"] + list(set(ex.get('category', 'General') for ex in examples))
                )
                
                filtered_examples = examples
                if search_term:
                    filtered_examples = [ex for ex in filtered_examples 
                                       if search_term.lower() in ex.get('input_text', '').lower()]
                if category_filter != "All":
                    filtered_examples = [ex for ex in filtered_examples 
                                       if ex.get('category') == category_filter]
                
                st.write(f"Showing {len(filtered_examples)} of {len(examples)} examples")
                
                for i, example in enumerate(filtered_examples[:20]):  # Limit display
                    with st.expander(f"Example {i+1}: {example.get('category', 'General')}"):
                        st.write("**Input:**")
                        st.write(example.get('input_text', ''))
                        
                        st.write("**Expected Output:**")
                        st.write(example.get('output_text', ''))
                        
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.write(f"**Difficulty:** {example.get('difficulty', 'Unknown')}")
                        with col2:
                            approval_status = "✅ Approved" if example.get('is_approved') else "⏳ Pending"
                            st.write(f"**Status:** {approval_status}")
                        with col3:
                            if example.get('tags'):
                                st.write(f"**Tags:** {', '.join(example['tags'])}")
                        
                        st.write(f"**Created:** {example.get('created_at', 'Unknown')}")
            else:
                st.info("No training examples found")

    # ============================================================================
    # TAB 2: KNOWLEDGE BASE
    # ============================================================================

    with tab2:
        st.header("Knowledge Base Management")
        
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.subheader("Upload Knowledge Document")
            
            with st.form("upload_doc"):
                doc_title = st.text_input("Document Title")
                doc_content = st.text_area("Content", height=200, placeholder="Enter document content or paste text...")
                
                doc_expert_id = st.selectbox(
                    "Target Expert",
                    options=[exp['id'] for exp in experts],
                    format_func=lambda x: next(exp['name'] for exp in experts if exp['id'] == x),
                    index=0
                )
                
                doc_type = st.selectbox(
                    "Document Type",
                    options=["Reference", "Protocol", "Guidelines", "Research", "FAQ"]
                )
                
                doc_tags = st.text_input("Tags (comma-separated)")
                
                if st.form_submit_button("Upload Document"):
                    if doc_title and doc_content:
                        doc_data = {
                            "expert_id": doc_expert_id,
                            "title": doc_title,
                            "content": doc_content,
                            "doc_type": doc_type.lower(),
                            "tags": [tag.strip() for tag in doc_tags.split(",") if tag.strip()]
                        }
                        
                        result = upload_knowledge_doc(doc_data)
                        if result:
                            st.success("Document uploaded successfully!")
                            st.rerun()
                        else:
                            st.error("Failed to upload document")
                    else:
                        st.error("Please provide title and content")
        
        with col2:
            st.subheader("Knowledge Documents")
            
            docs = get_knowledge_docs(st.session_state.selected_expert_id)
            
            if docs:
                doc_search = st.text_input("Search documents", placeholder="Search by title or content...")
                doc_type_filter = st.selectbox(
                    "Filter by Type",
                    options=["All"] + list(set(doc.get('doc_type', 'reference') for doc in docs))
                )
                
                filtered_docs = docs
                if doc_search:
                    filtered_docs = [doc for doc in filtered_docs 
                                   if doc_search.lower() in doc.get('title', '').lower() or
                                      doc_search.lower() in doc.get('content', '').lower()]
                if doc_type_filter != "All":
                    filtered_docs = [doc for doc in filtered_docs 
                                   if doc.get('doc_type') == doc_type_filter]
                
                st.write(f"Showing {len(filtered_docs)} of {len(docs)} documents")
                
                for doc in filtered_docs[:15]:  # Limit display
                    with st.expander(f"{doc.get('title', 'Untitled')} ({doc.get('doc_type', 'reference')})"):
                        content_preview = doc.get('content', '')[:300]
                        if len(doc.get('content', '')) > 300:
                            content_preview += "..."
                        st.write(content_preview)
                        
                        if doc.get('tags'):
                            st.write(f"**Tags:** {', '.join(doc['tags'])}")
                        
                        st.write(f"**Uploaded:** {doc.get('created_at', 'Unknown')}")
            else:
                st.info("No knowledge documents found")

    # ============================================================================
    # TAB 3: FEEDBACK REVIEW
    # ============================================================================

    with tab3:
        st.header("Feedback Analysis")
        
        feedback_data = get_feedback_summary(st.session_state.selected_expert_id)
        
        # Summary metrics (always show with defaults)
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Feedback", feedback_data.get('total_feedback', 0) if feedback_data else 0)
        
        with col2:
            avg_rating = feedback_data.get('avg_rating', 0) if feedback_data else 0
            st.metric("Average Rating", f"{avg_rating:.1f}/5.0")
        
        with col3:
            recent_corrections = len(feedback_data.get('recent_corrections', [])) if feedback_data else 0
            st.metric("Recent Corrections", recent_corrections)
        
        with col4:
            improvement_trend = feedback_data.get('improvement_trend', 0) if feedback_data else 0
            st.metric("Trend", f"{improvement_trend:+.1f}%")
        
        if feedback_data:
            
            # Detailed feedback
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Top Feedback Tags")
                top_tags = feedback_data.get('top_tags', {})
                if top_tags:
                    for tag, count in list(top_tags.items())[:10]:
                        st.write(f"• {tag}: {count} mentions")
                else:
                    st.info("No feedback tags available")
            
            with col2:
                st.subheader("Recent Corrections")
                corrections = feedback_data.get('recent_corrections', [])
                if corrections:
                    for correction in corrections[:10]:
                        st.write(f"• {correction}")
                else:
                    st.info("No recent corrections")
            
            # Rating distribution
            if 'rating_distribution' in feedback_data:
                st.subheader("Rating Distribution")
                rating_dist = feedback_data['rating_distribution']
                
                # Simple bar chart using columns
                for rating in [5, 4, 3, 2, 1]:
                    count = rating_dist.get(str(rating), 0)
                    percentage = (count / feedback_data.get('total_feedback', 1)) * 100
                    st.write(f"⭐ {rating}: {count} ({percentage:.1f}%)")
        else:
            st.info("No feedback data available")

    # ============================================================================
    # TAB 4: EXPERT CONFIGURATION
    # ============================================================================

    with tab4:
        st.header("Expert Configuration")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Basic Configuration")
            
            if st.session_state.selected_expert:
                expert = st.session_state.selected_expert
                
                with st.form("expert_config"):
                    model_name = st.selectbox(
                        "Base Model",
                        options=["gpt-4", "gpt-3.5-turbo", "claude-3-sonnet", "claude-3-haiku"],
                        index=0
                    )
                    
                    temperature = st.slider("Temperature", 0.0, 2.0, 0.7, 0.1)
                    max_tokens = st.number_input("Max Tokens", 100, 4000, 1000, 100)
                    
                    system_prompt = st.text_area(
                        "System Prompt",
                        value=expert.get('system_prompt', ''),
                        height=150,
                        placeholder="Enter the system prompt for this expert..."
                    )
                    
                    is_active = st.checkbox("Expert Active", value=expert.get('is_active', True))
                    
                    if st.form_submit_button("Update Configuration"):
                        config_data = {
                            "model_name": model_name,
                            "temperature": temperature,
                            "max_tokens": max_tokens,
                            "system_prompt": system_prompt,
                            "is_active": is_active
                        }
                        
                        if update_expert_config(expert['id'], config_data):
                            st.success("Configuration updated successfully!")
                            st.rerun()
                        else:
                            st.error("Failed to update configuration")
            
            with col2:
                st.subheader("Training Data Export")
                
                export_format = st.selectbox(
                    "Export Format",
                    options=["jsonl", "csv", "json"],
                    index=0
                )
                
                if st.button("Export Training Data"):
                    with st.spinner("Exporting data..."):
                        export_data = export_training_data(expert['id'], export_format)
                        
                        if export_data:
                            st.download_button(
                                label=f"Download {export_format.upper()}",
                                data=export_data,
                                file_name=f"{expert['name']}_training_data.{export_format}",
                                mime=f"application/{export_format}"
                            )
                        else:
                            st.error("Failed to export training data")
                
                if st.session_state.selected_expert:
                    st.subheader("Expert Statistics")
                    
                    # Get expert-specific stats
                    examples = get_training_examples(expert['id'])
                    docs = get_knowledge_docs(expert['id'])
                    feedback = get_feedback_summary(expert['id'])
                    
                    st.metric("Training Examples", len(examples))
                    st.metric("Knowledge Documents", len(docs))
                    st.metric("Total Feedback", feedback.get('total_feedback', 0))
                    
                    if feedback.get('avg_rating'):
                        st.metric("Average Rating", f"{feedback['avg_rating']:.1f}/5.0")
                
                else:
                    st.info("Select a specific expert to configure settings")
            
            # Global statistics when "All Experts" is selected
            st.subheader("Global Statistics")
            
            total_examples = len(get_training_examples())
            total_docs = len(get_knowledge_docs())
            global_feedback = get_feedback_summary()
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("Total Training Examples", total_examples)
            
            with col2:
                st.metric("Total Knowledge Documents", total_docs)
            
            with col3:
                st.metric("Total Feedback", global_feedback.get('total_feedback', 0))
            
            # Expert summary table
            if experts:
                st.subheader("Expert Overview")
                
                expert_summary = []
                for expert in experts:
                    examples_count = len(get_training_examples(expert['id']))
                    docs_count = len(get_knowledge_docs(expert['id']))
                    feedback_summary = get_feedback_summary(expert['id'])
                    
                    expert_summary.append({
                        "Name": expert['name'],
                        "Specializations": ', '.join(expert['specializations']),
                        "Active": "✅" if expert.get('is_active', True) else "❌",
                        "Training Examples": examples_count,
                        "Knowledge Docs": docs_count,
                        "Avg Rating": f"{feedback_summary.get('avg_rating', 0):.1f}/5.0" if feedback_summary.get('avg_rating') else "N/A"
                    })
                
                df = pd.DataFrame(expert_summary)
                st.dataframe(df, use_container_width=True)
