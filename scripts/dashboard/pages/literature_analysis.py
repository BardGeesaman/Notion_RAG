"""Literature Critical Analysis page."""
from __future__ import annotations

import json
from uuid import UUID

import pandas as pd
import streamlit as st

from amprenta_rag.analysis.literature_critique import (
    detect_contradictions,
    extract_unanswered_questions,
    generate_critique,
)
from amprenta_rag.database.models import LiteratureCritique
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_literature_analysis_page() -> None:
    """Render the Literature Analysis page."""
    st.header("📚 Literature Critical Analysis")
    st.markdown("Analyze scientific literature for strengths, weaknesses, gaps, and contradictions.")
    
    tab1, tab2, tab3 = st.tabs(["Critique", "Unanswered Questions", "Contradiction Finder"])
    
    with tab1:
        render_critique_tab()
    with tab2:
        render_questions_tab()
    with tab3:
        render_contradiction_tab()


def render_critique_tab() -> None:
    """Render the Critique tab."""
    st.subheader("Generate Critical Analysis")
    st.markdown("Paste scientific text to analyze for strengths, weaknesses, limitations, and methodology quality.")
    
    text = st.text_area(
        "Scientific Text",
        height=300,
        placeholder="Paste paper abstract, methods section, or full text here...",
        key="critique_text"
    )
    
    col1, col2 = st.columns([1, 1])
    with col1:
        analyze_btn = st.button("🔍 Analyze", type="primary", use_container_width=True)
    with col2:
        if st.session_state.get("critique_result"):
            save_btn = st.button("💾 Save Critique", use_container_width=True)
        else:
            save_btn = False
    
    if analyze_btn and text:
        with st.spinner("Analyzing text..."):
            with db_session() as db:
                result = generate_critique(text)
                st.session_state["critique_result"] = result
                st.session_state["critique_text"] = text
                st.rerun()
    
    if save_btn and st.session_state.get("critique_result"):
        user = get_current_user()
        if not user:
            st.error("You must be logged in to save critiques.")
        else:
            with db_session() as db:
                critique = LiteratureCritique(
                    source_type="text",
                    source_id=None,
                    source_text=st.session_state.get("critique_text", ""),
                    critique_type="critique",
                    content=st.session_state["critique_result"],
                    created_by_id=UUID(user.get("id")) if user.get("id") and user.get("id") != "test" else None,
                )
                db.add(critique)
                db.commit()
                st.success("Critique saved!")
                st.rerun()
    
    if st.session_state.get("critique_result"):
        result = st.session_state["critique_result"]
        
        # Methodology Score
        score = result.get("methodology_score", 50)
        st.metric("Methodology Score", f"{score}/100", delta=None)
        
        # Strengths
        strengths = result.get("strengths", [])
        if strengths:
            st.markdown("### ✅ Strengths")
            for i, strength in enumerate(strengths, 1):
                st.markdown(f"{i}. {strength}")
        
        # Weaknesses
        weaknesses = result.get("weaknesses", [])
        if weaknesses:
            st.markdown("### ⚠️ Weaknesses")
            for i, weakness in enumerate(weaknesses, 1):
                st.markdown(f"{i}. {weakness}")
        
        # Limitations
        limitations = result.get("limitations", [])
        if limitations:
            st.markdown("### 🔍 Limitations")
            for i, limitation in enumerate(limitations, 1):
                st.markdown(f"{i}. {limitation}")


def render_questions_tab() -> None:
    """Render the Unanswered Questions tab."""
    st.subheader("Extract Unanswered Questions")
    st.markdown("Identify research gaps and unanswered questions from scientific text.")
    
    text = st.text_area(
        "Scientific Text",
        height=300,
        placeholder="Paste text to analyze for gaps and open questions...",
        key="questions_text"
    )
    
    if st.button("❓ Extract Questions", type="primary"):
        if not text:
            st.error("Please enter text to analyze.")
        else:
            with st.spinner("Extracting questions..."):
                questions = extract_unanswered_questions(text)
                st.session_state["extracted_questions"] = questions
                st.rerun()
    
    if st.session_state.get("extracted_questions"):
        questions = st.session_state["extracted_questions"]
        st.markdown(f"### Found {len(questions)} Question(s)")
        for i, question in enumerate(questions, 1):
            st.markdown(f"{i}. {question}")


def render_contradiction_tab() -> None:
    """Render the Contradiction Finder tab."""
    st.subheader("Detect Contradictions")
    st.markdown("Compare two scientific texts to identify conflicting claims.")
    
    col1, col2 = st.columns(2)
    with col1:
        text1 = st.text_area(
            "Paper A",
            height=250,
            placeholder="First paper text...",
            key="contradiction_text1"
        )
    with col2:
        text2 = st.text_area(
            "Paper B",
            height=250,
            placeholder="Second paper text...",
            key="contradiction_text2"
        )
    
    if st.button("🔍 Compare", type="primary"):
        if not text1 or not text2:
            st.error("Please enter text in both fields.")
        else:
            with st.spinner("Comparing texts for contradictions..."):
                contradictions = detect_contradictions(text1, text2)
                st.session_state["contradictions"] = contradictions
                st.rerun()
    
    if st.session_state.get("contradictions"):
        contradictions = st.session_state["contradictions"]
        
        if not contradictions:
            st.info("No contradictions found between the two texts.")
        else:
            st.markdown(f"### Found {len(contradictions)} Contradiction(s)")
            
            # Display as table
            contradiction_data = []
            for c in contradictions:
                severity = c.get("severity", "medium")
                color = {
                    "high": "🔴",
                    "medium": "🟡",
                    "low": "🟢"
                }.get(severity.lower(), "⚪")
                
                contradiction_data.append({
                    "Severity": f"{color} {severity.upper()}",
                    "Topic": c.get("topic", "Unknown"),
                    "Claim A": c.get("claim_a", "")[:100] + "..." if len(c.get("claim_a", "")) > 100 else c.get("claim_a", ""),
                    "Claim B": c.get("claim_b", "")[:100] + "..." if len(c.get("claim_b", "")) > 100 else c.get("claim_b", ""),
                })
            
            df = pd.DataFrame(contradiction_data)
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            # Show details in expanders
            for i, c in enumerate(contradictions, 1):
                with st.expander(f"Contradiction {i}: {c.get('topic', 'Unknown')} ({c.get('severity', 'medium')})"):
                    st.markdown(f"**Claim A:** {c.get('claim_a', '')}")
                    st.markdown(f"**Claim B:** {c.get('claim_b', '')}")
                    st.markdown(f"**Severity:** {c.get('severity', 'medium')}")
