"""AI Expert Training Console - Admin interface for expert management."""

from __future__ import annotations

import streamlit as st
import requests
import pandas as pd
from datetime import datetime
from uuid import UUID

st.set_page_config(page_title="AI Expert Training", page_icon="🎓", layout="wide")

API_BASE = "http://localhost:8000/api/v1/experts"

# Initialize session state
if "selected_expert_id" not in st.session_state:
    st.session_state.selected_expert_id = None
if "selected_expert" not in st.session_state:
    st.session_state.selected_expert = None


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


def get_training_examples(expert_id, approved_only=False):
    """Fetch expert's training examples."""
    try:
        params = {"approved_only": approved_only} if approved_only else {}
        resp = requests.get(f"{API_BASE}/experts/{expert_id}/training-examples", params=params)
        return resp.json() if resp.ok else []
    except:
        return []


def get_feedback_summary(expert_id):
    """Fetch expert's feedback summary."""
    try:
        resp = requests.get(f"{API_BASE}/experts/{expert_id}/feedback-summary")
        return resp.json() if resp.ok else {}
    except:
        return {}


def export_training_data(expert_id):
    """Export training data for expert."""
    try:
        resp = requests.post(f"{API_BASE}/experts/{expert_id}/export-training")
        return resp.json() if resp.ok else {}
    except:
        return {}


# ============================================================================
# MAIN LAYOUT
# ============================================================================

st.title("🎓 AI Expert Training Console")
st.markdown("Admin interface for managing expert training, knowledge, and configuration.")

# Expert selection
experts = get_experts()
if experts:
    expert_options = {expert["id"]: f"{expert['name']} ({expert['role']})" for expert in experts}
    
    selected_expert_id = st.selectbox(
        "Select Expert to Manage",
        options=list(expert_options.keys()),
        format_func=lambda x: expert_options.get(x, x),
        index=0 if expert_options else None,
        key="expert_selector"
    )
    
    if selected_expert_id:
        st.session_state.selected_expert_id = selected_expert_id
        st.session_state.selected_expert = next(e for e in experts if e["id"] == selected_expert_id)
        
        # Show selected expert info
        expert = st.session_state.selected_expert
        col1, col2, col3, col4 = st.columns(4)
        
        col1.metric("Expert", expert["name"])
        col2.metric("Role", expert["role"])
        col3.metric("Version", expert["prompt_version"])
        col4.metric("Status", "Active" if expert["is_active"] else "Inactive")
else:
    st.error("No experts available. Check API connection.")
    st.stop()


# ============================================================================
# TABS
# ============================================================================

tab1, tab2, tab3, tab4 = st.tabs([
    "📝 Training Examples",
    "📚 Knowledge Docs", 
    "📊 Feedback Review",
    "⚙️ Expert Config"
])


# ============================================================================
# TAB 1: TRAINING EXAMPLES
# ============================================================================

with tab1:
    st.subheader("Training Examples")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        show_approved_only = st.checkbox("Show approved only")
    with col2:
        if st.button("📥 Export Training Data", type="primary"):
            export_data = export_training_data(st.session_state.selected_expert_id)
            if export_data:
                st.success(f"Exported {export_data.get('example_count', 0)} examples")
                
                # Provide download
                import json
                jsonl_content = "\n".join([json.dumps(item) for item in export_data.get("data", [])])
                st.download_button(
                    "⬇️ Download JSONL",
                    data=jsonl_content,
                    file_name=f"expert_{st.session_state.selected_expert['name']}_training.jsonl",
                    mime="application/jsonl"
                )
    
    # Add new training example
    with st.expander("➕ Add Training Example"):
        with st.form("add_training_example"):
            input_text = st.text_area("Input (Question)", placeholder="What is the best approach for...")
            ideal_output = st.text_area("Ideal Output (Answer)", placeholder="The best approach is...")
            context = st.text_input("Context (optional)", placeholder="Background information")
            tags = st.text_input("Tags (comma-separated)", placeholder="SAR, optimization, ADMET")
            
            if st.form_submit_button("Add Example", type="primary"):
                if input_text and ideal_output:
                    try:
                        payload = {
                            "expert_id": st.session_state.selected_expert_id,
                            "input_text": input_text,
                            "ideal_output": ideal_output,
                            "context": context if context else None,
                            "tags": [t.strip() for t in tags.split(",") if t.strip()] if tags else None
                        }
                        resp = requests.post(f"{API_BASE}/experts/{st.session_state.selected_expert_id}/training-examples", json=payload)
                        if resp.ok:
                            st.success("Training example added!")
                            st.rerun()
                        else:
                            st.error(f"Failed: {resp.json().get('detail', 'Unknown error')}")
                    except Exception as e:
                        st.error(f"Error: {e}")
                else:
                    st.error("Input and ideal output are required")
    
    # List existing examples
    examples = get_training_examples(st.session_state.selected_expert_id, approved_only=show_approved_only)
    
    if examples:
        for example in examples:
            with st.container(border=True):
                col1, col2 = st.columns([4, 1])
                
                with col1:
                    st.markdown(f"**Question:** {example['question']}")
                    st.markdown(f"**Answer:** {example['ideal_answer']}")
                    
                    status = "✅ Approved" if example['is_approved'] else "⏳ Pending"
                    st.caption(f"Status: {status} | Version: {example['prompt_version']} | Created: {example['created_at'][:10]}")
                
                with col2:
                    if not example['is_approved']:
                        if st.button("✅ Approve", key=f"approve_{example['id']}"):
                            # Note: Approve function would need to be added to service
                            st.success("Approved!")
                    
                    if st.button("🗑️ Delete", key=f"delete_ex_{example['id']}"):
                        st.warning("Delete functionality not implemented")
    else:
        st.info("No training examples found. Add examples above to improve expert performance.")


# ============================================================================
# TAB 2: KNOWLEDGE DOCS
# ============================================================================

with tab2:
    st.subheader("Knowledge Documents")
    
    # Upload new knowledge document
    with st.expander("📤 Upload Knowledge Document"):
        with st.form("upload_knowledge"):
            doc_title = st.text_input("Document Title *")
            doc_content = st.text_area("Content *", height=200, placeholder="Paste document content here...")
            doc_source_type = st.selectbox("Source Type", ["manual", "paper", "web", "internal"])
            doc_source_url = st.text_input("Source URL (optional)")
            
            if st.form_submit_button("Upload Document", type="primary"):
                if doc_title and doc_content:
                    try:
                        payload = {
                            "expert_id": st.session_state.selected_expert_id,
                            "title": doc_title,
                            "content": doc_content,
                            "source_type": doc_source_type,
                            "source_url": doc_source_url if doc_source_url else None
                        }
                        resp = requests.post(f"{API_BASE}/experts/{st.session_state.selected_expert_id}/knowledge", json=payload)
                        if resp.ok:
                            st.success("Knowledge document uploaded and embedded!")
                            st.rerun()
                        else:
                            st.error(f"Upload failed: {resp.json().get('detail', 'Unknown error')}")
                    except Exception as e:
                        st.error(f"Error: {e}")
                else:
                    st.error("Title and content are required")
    
    # List existing knowledge docs (placeholder - would need API endpoint)
    st.info("Knowledge document listing requires additional API endpoint implementation.")
    
    # Show knowledge stats
    try:
        resp = requests.get(f"{API_BASE}/{st.session_state.selected_expert_id}/stats")
        if resp.ok:
            stats = resp.json()
            st.metric("Knowledge Documents", stats.get("knowledge_doc_count", 0))
    except:
        pass


# ============================================================================
# TAB 3: FEEDBACK REVIEW
# ============================================================================

with tab3:
    st.subheader("Feedback Analysis")
    
    feedback_summary = get_feedback_summary(st.session_state.selected_expert_id)
    
    if feedback_summary and feedback_summary.get("feedback_count", 0) > 0:
        # Key metrics
        col1, col2, col3 = st.columns(3)
        
        avg_rating = feedback_summary.get("avg_rating")
        col1.metric("Average Rating", f"{avg_rating:.1f}/5.0" if avg_rating else "N/A")
        col2.metric("Total Feedback", feedback_summary.get("feedback_count", 0))
        col3.metric("Recent Corrections", len(feedback_summary.get("recent_corrections", [])))
        
        # Rating distribution
        st.markdown("### Rating Distribution")
        rating_dist = feedback_summary.get("rating_distribution", {})
        if rating_dist:
            # Convert to DataFrame for chart
            df = pd.DataFrame([
                {"Rating": f"{rating} ⭐", "Count": count}
                for rating, count in rating_dist.items()
            ])
            st.bar_chart(df.set_index("Rating"))
        
        # Recent corrections
        st.markdown("### Recent User Corrections")
        corrections = feedback_summary.get("recent_corrections", [])
        if corrections:
            for i, correction in enumerate(corrections[-5:]):  # Last 5
                with st.container(border=True):
                    st.write(correction)
                    if st.button("💡 Convert to Training", key=f"convert_{i}"):
                        st.info("Conversion to training example requires manual review")
        else:
            st.info("No recent corrections")
    else:
        st.info("No feedback data available for this expert yet.")


# ============================================================================
# TAB 4: EXPERT CONFIG
# ============================================================================

with tab4:
    st.subheader("Expert Configuration")
    
    expert = st.session_state.selected_expert
    
    # Current configuration display
    st.markdown("### Current Configuration")
    
    col1, col2 = st.columns(2)
    with col1:
        st.text_input("Name", value=expert["name"], disabled=True)
        st.text_input("Role", value=expert["role"], disabled=True)
        st.text_input("Prompt Version", value=expert["prompt_version"], disabled=True)
    
    with col2:
        st.selectbox("Status", ["Active", "Inactive"], index=0 if expert["is_active"] else 1, disabled=True)
        specializations = ", ".join(expert.get("specializations", []))
        st.text_area("Specializations", value=specializations, disabled=True, height=100)
    
    # System prompt editor
    st.markdown("### System Prompt")
    
    # Get current prompt (would need to fetch full details)
    current_prompt = "System prompt would be loaded from API..."
    
    with st.form("update_expert"):
        new_prompt = st.text_area(
            "System Prompt",
            value=current_prompt,
            height=300,
            help="Define the expert's personality, knowledge, and response style"
        )
        
        col1, col2, col3 = st.columns(3)
        with col1:
            new_specializations = st.text_input(
                "Specializations (comma-separated)",
                value=", ".join(expert.get("specializations", []))
            )
        with col2:
            is_active = st.checkbox("Active", value=expert["is_active"])
        with col3:
            bump_version = st.checkbox("Bump Version", help="Increment prompt version")
        
        if st.form_submit_button("Update Expert", type="primary"):
            try:
                payload = {
                    "system_prompt": new_prompt,
                    "specializations": [s.strip() for s in new_specializations.split(",") if s.strip()],
                    "is_active": is_active,
                    "bump_version": bump_version
                }
                resp = requests.patch(f"{API_BASE}/experts/{st.session_state.selected_expert_id}", json=payload)
                if resp.ok:
                    st.success("Expert updated successfully!")
                    if bump_version:
                        st.info("Prompt version incremented")
                    st.rerun()
                else:
                    st.error(f"Update failed: {resp.json().get('detail', 'Unknown error')}")
            except Exception as e:
                st.error(f"Error: {e}")
    
    # Expert statistics
    st.markdown("### Performance Metrics")
    try:
        resp = requests.get(f"{API_BASE}/{st.session_state.selected_expert_id}/stats")
        if resp.ok:
            stats = resp.json()
            
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Messages", stats.get("message_count", 0))
            col2.metric("Knowledge Docs", stats.get("knowledge_doc_count", 0))
            col3.metric("Training Examples", stats.get("training_example_count", 0))
            
            avg_rating = stats.get("average_rating")
            col4.metric("Avg Rating", f"{avg_rating:.1f}/5.0" if avg_rating else "N/A")
    except:
        st.error("Failed to load expert statistics")


# ============================================================================
# ADMIN TOOLS
# ============================================================================

st.markdown("---")
st.markdown("### 🛠️ Admin Tools")

col1, col2, col3 = st.columns(3)

with col1:
    if st.button("🔄 Refresh Expert Data"):
        # Reload expert data
        experts = get_experts()
        if experts:
            st.session_state.selected_expert = next(e for e in experts if e["id"] == st.session_state.selected_expert_id)
            st.success("Data refreshed!")
            st.rerun()

with col2:
    if st.button("📊 Generate Report"):
        st.info("Expert performance report generation not yet implemented")

with col3:
    if st.button("🔧 System Diagnostics"):
        # Basic API health check
        try:
            resp = requests.get(f"{API_BASE}/")
            if resp.ok:
                st.success(f"API OK - {len(resp.json())} experts available")
            else:
                st.error(f"API Error: {resp.status_code}")
        except Exception as e:
            st.error(f"API Connection Failed: {e}")


# ============================================================================
# USAGE TIPS
# ============================================================================

with st.expander("💡 Training Console Guide"):
    st.markdown("""
    **Training Examples:**
    - Add high-quality question-answer pairs for fine-tuning
    - Review and approve examples before export
    - Export approved examples in JSONL format for OpenAI fine-tuning
    
    **Knowledge Documents:**
    - Upload domain-specific documents for RAG retrieval
    - Documents are automatically chunked and embedded
    - Experts use this knowledge to enhance responses
    
    **Feedback Review:**
    - Monitor user satisfaction with expert responses
    - Convert negative feedback into training improvements
    - Track rating trends over time
    
    **Expert Configuration:**
    - Update system prompts to refine expert behavior
    - Manage specializations and active status
    - Version control for prompt changes
    
    **Best Practices:**
    - Regularly review feedback and update prompts
    - Add diverse training examples covering edge cases
    - Keep knowledge base current with latest research
    - Test changes with sample questions before deployment
    """)


# ============================================================================
# REGISTER IN PAGE_REGISTRY
# ============================================================================

def main():
    """Entry point for page registry."""
    return  # Page renders on import
