"""Q&A Tracker page for saving and re-running questions."""
from __future__ import annotations

from datetime import datetime
from typing import List

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import SavedQuestion, SavedAnswer
from amprenta_rag.database.session import db_session


def _run_rag_query(question: str) -> str:
    """Stub for RAG query execution. Replace with real call if available."""
    # TODO: integrate actual RAG query service
    return f"Stubbed answer for: {question}"


def render_qa_tracker_page():
    st.title("🧠 Q&A Tracker")
    tab1, tab2, tab3, tab4 = st.tabs(["Ask & Save", "My Questions", "Re-run", "Export"])

    with tab1:
        st.subheader("Ask & Save")
        question = st.text_area("Question", height=100)
        if st.button("Run", type="primary"):
            if not question.strip():
                st.error("Please enter a question")
            else:
                answer = _run_rag_query(question)
                st.session_state["qa_last_question"] = question
                st.session_state["qa_last_answer"] = answer
        if st.session_state.get("qa_last_answer"):
            st.markdown("**Answer:**")
            st.write(st.session_state["qa_last_answer"])
            if st.button("Save", key="qa_save_btn"):
                with db_session() as db:
                    try:
                        q = SavedQuestion(
                            question_text=st.session_state.get("qa_last_question", question),
                            tags=None,
                            created_at=datetime.utcnow(),
                        )
                        db.add(q)
                        db.flush()
                        ans = SavedAnswer(
                            question_id=q.id,
                            answer_text=st.session_state["qa_last_answer"],
                            evidence=None,
                            model_used="stub",
                            version=1,
                            confidence_score=None,
                            created_at=datetime.utcnow(),
                        )
                        db.add(ans)
                        db.commit()
                        st.success("Saved question and answer")
                    except Exception as e:
                        db.rollback()
                        st.error(f"Failed to save: {e}")

    with tab2:
        st.subheader("My Questions")
        search = st.text_input("Search questions")
        with db_session() as db:
            query = db.query(SavedQuestion).filter(SavedQuestion.is_archived == False)
            if search:
                like = f"%{search}%"
                query = query.filter(SavedQuestion.question_text.ilike(like))
            questions = query.order_by(SavedQuestion.created_at.desc()).limit(200).all()
            if not questions:
                st.info("No questions found")
            else:
                options = {q.question_text[:120]: q for q in questions}
                selected_text = st.selectbox("Select question", list(options.keys()))
                selected = options[selected_text]
                answers = selected.answers or []
                st.markdown(f"**Question:** {selected.question_text}")
                if answers:
                    for ans in sorted(answers, key=lambda a: a.version, reverse=True):
                        with st.expander(f"Answer v{ans.version}"):
                            st.write(ans.answer_text)
                            if ans.evidence:
                                st.json(ans.evidence)
                            st.caption(f"Model: {ans.model_used or 'n/a'} | Confidence: {ans.confidence_score or 'n/a'}")
                else:
                    st.info("No answers saved for this question.")

    with tab3:
        st.subheader("Re-run")
        with db_session() as db:
            questions = db.query(SavedQuestion).order_by(SavedQuestion.created_at.desc()).limit(200).all()
            if not questions:
                st.info("No questions available")
            else:
                qmap = {q.question_text[:120]: q for q in questions}
                sel_txt = st.selectbox("Select question to re-run", list(qmap.keys()))
                sel_q = qmap[sel_txt]
                if st.button("Re-run", key="qa_rerun_btn"):
                    new_answer = _run_rag_query(sel_q.question_text)
                    st.write("**New Answer:**")
                    st.write(new_answer)
                    if st.button("Save as new version", key="qa_rerun_save"):
                        try:
                            max_ver = max([a.version for a in sel_q.answers], default=1)
                            ans = SavedAnswer(
                                question_id=sel_q.id,
                                answer_text=new_answer,
                                evidence=None,
                                model_used="stub",
                                version=max_ver + 1,
                                confidence_score=None,
                                created_at=datetime.utcnow(),
                            )
                            db.add(ans)
                            db.commit()
                            st.success("Saved new answer version")
                        except Exception as e:
                            db.rollback()
                            st.error(f"Failed to save: {e}")

    with tab4:
        st.subheader("Export")
        with db_session() as db:
            questions = db.query(SavedQuestion).order_by(SavedQuestion.created_at.desc()).limit(500).all()
            if not questions:
                st.info("No questions to export")
            else:
                selection = st.multiselect("Select questions", [q.question_text[:120] for q in questions])
                if st.button("Download CSV", key="qa_export_btn"):
                    rows: List[dict] = []
                    for q in questions:
                        if q.question_text[:120] in selection:
                            for ans in q.answers or []:
                                rows.append(
                                    {
                                        "question": q.question_text,
                                        "answer": ans.answer_text,
                                        "version": ans.version,
                                        "model": ans.model_used,
                                        "confidence": ans.confidence_score,
                                        "created_at": ans.created_at.isoformat() if ans.created_at else None,
                                    }
                                )
                    if rows:
                        df = pd.DataFrame(rows)
                        csv = df.to_csv(index=False)
                        st.download_button("Download CSV", data=csv, file_name="qa_export.csv", mime="text/csv")
                    else:
                        st.warning("No rows to export.")


if __name__ == "__main__":
    render_qa_tracker_page()
