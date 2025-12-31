"""Electronic Lab Notebook page for the Streamlit dashboard."""

from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

# Import models - use direct module import to avoid caching issues
import amprenta_rag.database.models as db_models

# Access models directly from module
LabNotebookEntry = db_models.LabNotebookEntry
LabNotebookEntryAssociation = db_models.LabNotebookEntryAssociation
Dataset = db_models.Dataset
Program = db_models.Program
Experiment = db_models.Experiment

from scripts.dashboard.db_session import db_session
from scripts.dashboard.components.annotation_panel import render_annotation_panel, render_annotation_indicator


def render_lab_notebook_page() -> None:
    """
    Render the Electronic Lab Notebook page.

    Features:
    - Create new notebook entries
    - Link entries to datasets, experiments, programs
    - Search and filter entries
    - View entry history
    - Export entries
    """
    st.header("📓 Electronic Lab Notebook")
    st.markdown("Create and manage lab notebook entries linked to your experiments, datasets, and programs.")

    # Tab selection
    tab1, tab2, tab3, tab4 = st.tabs(["New Entry", "Browse Entries", "Search", "Templates"])

    # Tab 1: Create New Entry
    with tab1:
        st.subheader("Create New Lab Notebook Entry")

        # Entry form
        entry_title = st.text_input(
            "Entry Title *",
            placeholder="e.g., Lipidomics analysis of ALS patient samples",
            help="A descriptive title for this entry",
        )

        entry_type = st.selectbox(
            "Entry Type",
            ["experiment", "observation", "protocol", "note", "analysis", "other"],
            help="Categorize this entry",
        )

        # Content editor
        entry_content = st.text_area(
            "Content *",
            height=300,
            placeholder="Enter your notes, observations, or experimental details here...",
            help="Main content of the entry",
        )

        # Tags
        tags_input = st.text_input(
            "Tags (comma-separated)",
            placeholder="e.g., lipidomics, ALS, patient samples",
            help="Add tags for easier searching",
        )

        # Link to entities
        st.markdown("### Link to Entities (Optional)")
        with db_session() as db:
            # Link to Programs
            programs = db.query(Program).order_by(Program.name).all()
            if programs:
                selected_programs = st.multiselect(
                    "Link to Programs",
                    options=[(p.id, p.name) for p in programs],
                    format_func=lambda x: x[1],
                    help="Select programs this entry relates to",
                )
            else:
                selected_programs = []
                st.info("No programs available.")

            # Link to Experiments
            experiments = db.query(Experiment).order_by(Experiment.name).all()
            if experiments:
                selected_experiments = st.multiselect(
                    "Link to Experiments",
                    options=[(e.id, e.name) for e in experiments],
                    format_func=lambda x: x[1],
                    help="Select experiments this entry relates to",
                )
            else:
                selected_experiments = []
                st.info("No experiments available.")

            # Link to Datasets
            datasets = db.query(Dataset).order_by(Dataset.name).all()
            if datasets:
                selected_datasets = st.multiselect(
                    "Link to Datasets",
                    options=[(d.id, d.name) for d in datasets],
                    format_func=lambda x: x[1],
                    help="Select datasets this entry relates to",
                )
            else:
                selected_datasets = []
                st.info("No datasets available.")

        # File attachments
        st.markdown("### Attachments")
        uploaded_files = st.file_uploader(
            "Upload files (optional)",
            type=["pdf", "txt", "csv", "tsv", "xlsx", "xls", "png", "jpg", "jpeg", "png"],
            accept_multiple_files=True,
            help="Upload files to attach to this entry",
        )

        if uploaded_files:
            st.info(f"📎 {len(uploaded_files)} file(s) selected")
            for uploaded_file in uploaded_files:
                st.caption(f"  - {uploaded_file.name} ({uploaded_file.size:,} bytes)")

        # Save button
        if st.button("💾 Save Entry", type="primary"):
            if not entry_title or not entry_content:
                st.error("Please provide both a title and content for the entry.")
            else:
                with st.spinner("Saving entry..."):
                    try:
                        with db_session() as db:
                            # Parse tags
                            tags = [tag.strip() for tag in tags_input.split(",") if tag.strip()] if tags_input else []

                            # Create entry first to get ID
                            entry = LabNotebookEntry(
                                title=entry_title,
                                content=entry_content,
                                entry_type=entry_type,
                                tags=tags if tags else None,
                            )
                            db.add(entry)
                            db.flush()  # Get the ID

                            # Handle file attachments
                            attachments_data = []
                            if uploaded_files:
                                # Create attachments directory using entry ID
                                attachments_dir = Path("data") / "eln_attachments" / str(entry.id)
                                attachments_dir.mkdir(parents=True, exist_ok=True)

                                # Save files
                                for uploaded_file in uploaded_files:
                                    file_path = attachments_dir / uploaded_file.name
                                    with open(file_path, "wb") as f:
                                        f.write(uploaded_file.getbuffer())

                                    attachments_data.append(
                                        {
                                            "filename": uploaded_file.name,
                                            "path": f"eln_attachments/{entry.id}/{uploaded_file.name}",
                                            "size": uploaded_file.size,
                                            "type": uploaded_file.type or "unknown",
                                        }
                                    )

                                # Update entry with attachments
                                entry.attachments = attachments_data
                                db.flush()

                            # Create associations
                            for program_id, _ in selected_programs:
                                assoc = LabNotebookEntryAssociation(
                                    entry_id=entry.id,
                                    entity_type="program",
                                    entity_id=program_id,
                                )
                                db.add(assoc)

                            for experiment_id, _ in selected_experiments:
                                assoc = LabNotebookEntryAssociation(
                                    entry_id=entry.id,
                                    entity_type="experiment",
                                    entity_id=experiment_id,
                                )
                                db.add(assoc)

                            for dataset_id, _ in selected_datasets:
                                assoc = LabNotebookEntryAssociation(
                                    entry_id=entry.id,
                                    entity_type="dataset",
                                    entity_id=dataset_id,
                                )
                                db.add(assoc)

                            db.commit()

                            st.success(f"✅ Entry '{entry_title}' saved successfully!")
                            st.rerun()
                    except Exception as e:
                        st.error(f"❌ Error saving entry: {str(e)}")
                        st.exception(e)

    # Tab 2: Browse Entries
    with tab2:
        st.subheader("Browse Lab Notebook Entries")

        # Filters
        col1, col2 = st.columns(2)
        with col1:
            entry_type_filter = st.selectbox(
                "Filter by Type",
                ["All", "experiment", "observation", "protocol", "note", "analysis", "other"],
            )
        with col2:
            sort_order = st.selectbox(
                "Sort by",
                ["Newest First", "Oldest First", "Title A-Z", "Title Z-A"],
            )

        with db_session() as db:
            query = db.query(LabNotebookEntry)

            if entry_type_filter != "All":
                query = query.filter(LabNotebookEntry.entry_type == entry_type_filter)

            # Apply sorting
            if sort_order == "Newest First":
                query = query.order_by(LabNotebookEntry.created_at.desc())
            elif sort_order == "Oldest First":
                query = query.order_by(LabNotebookEntry.created_at.asc())
            elif sort_order == "Title A-Z":
                query = query.order_by(LabNotebookEntry.title.asc())
            else:  # Title Z-A
                query = query.order_by(LabNotebookEntry.title.desc())

            entries = query.all()

            st.metric("Total Entries", len(entries))

            if entries:
                for entry in entries:
                    # Add annotation indicator to entry title
                    entry_title_with_annotations = f"**{entry.title}** - {entry.created_at.strftime('%Y-%m-%d %H:%M')}"
                    
                    with st.expander(entry_title_with_annotations):
                        # Show annotation panel in sidebar if this entry is selected
                        if st.session_state.get("annotation_context", {}).get("entity_id") == str(entry.id):
                            with st.sidebar:
                                render_annotation_panel(
                                    entity_type="notebook",
                                    entity_id=entry.id,
                                    position_type=st.session_state["annotation_context"].get("position_type"),
                                    position_data=st.session_state["annotation_context"].get("position_data"),
                                )
                        
                        col1, col2 = st.columns([3, 1])
                        with col1:
                            st.write(f"**Type:** {entry.entry_type or 'N/A'}")
                            if entry.tags:
                                st.write(f"**Tags:** {', '.join(entry.tags)}")

                            # Show attachments
                            if entry.attachments:
                                st.markdown("**Attachments:**")
                                for att in entry.attachments:
                                    file_path = Path("data") / att.get("path", "")
                                    if file_path.exists():
                                        with open(file_path, "rb") as f:
                                            st.download_button(
                                                label=f"📎 {att.get('filename', 'File')}",
                                                data=f.read(),
                                                file_name=att.get("filename", "file"),
                                                mime=att.get("type", "application/octet-stream"),
                                                key=f"download_{entry.id}_{att.get('filename')}",
                                            )
                                    else:
                                        st.caption(f"📎 {att.get('filename', 'File')} (file not found)")

                            st.markdown("---")
                            
                            # Split content into simulated "cells" for annotation purposes
                            content_lines = entry.content.split('\n')
                            for i, line in enumerate(content_lines):
                                if line.strip():  # Only show non-empty lines
                                    col_content, col_annotation = st.columns([10, 1])
                                    
                                    with col_content:
                                        st.markdown(line)
                                    
                                    with col_annotation:
                                        # Add annotation indicator for this "cell"
                                        if render_annotation_indicator(
                                            entity_type="notebook",
                                            entity_id=entry.id,
                                            position_type="cell",
                                            position_data={"cell_index": i},
                                            label="📝"
                                        ):
                                            # Button was clicked, annotation panel will be shown
                                            pass

                            # Edit button
                            if st.button("✏️ Edit Entry", key=f"edit_{entry.id}"):
                                st.session_state[f"editing_entry_{entry.id}"] = True
                                st.rerun()

                            # Export button
                            col_exp1, _ = st.columns(2)
                            with col_exp1:
                                # CSV export
                                entry_data = {
                                    "ID": str(entry.id),
                                    "Title": entry.title,
                                    "Content": entry.content,
                                    "Type": entry.entry_type or "",
                                    "Tags": ", ".join(entry.tags) if entry.tags else "",
                                    "Created": entry.created_at.strftime("%Y-%m-%d %H:%M:%S"),
                                    "Updated": entry.updated_at.strftime("%Y-%m-%d %H:%M:%S"),
                                }
                                df_entry = pd.DataFrame([entry_data])
                                csv = df_entry.to_csv(index=False)
                                st.download_button(
                                    label="📥 Export CSV",
                                    data=csv,
                                    file_name=f"eln_entry_{entry.id}.csv",
                                    mime="text/csv",
                                    key=f"export_csv_{entry.id}",
                                )

                        with col2:
                            st.caption(f"**ID:** `{entry.id}`")
                            st.caption(f"**Created:** {entry.created_at.strftime('%Y-%m-%d %H:%M')}")
                            st.caption(f"**Updated:** {entry.updated_at.strftime('%Y-%m-%d %H:%M')}")

                            # Show linked entities
                            linked_entities = (
                                db.query(LabNotebookEntryAssociation)
                                .filter(LabNotebookEntryAssociation.entry_id == entry.id)
                                .all()
                            )

                            if linked_entities:
                                st.markdown("**Linked Entities:**")
                                for assoc in linked_entities:
                                    st.caption(f"- {assoc.entity_type}: `{assoc.entity_id}`")

                            # Delete button
                            if st.button("🗑️ Delete", key=f"delete_{entry.id}", type="secondary"):
                                if st.session_state.get(f"confirm_delete_{entry.id}", False):
                                    try:
                                        # Delete associations
                                        for assoc in linked_entities:
                                            db.delete(assoc)

                                        # Delete attachments
                                        if entry.attachments:
                                            attachments_dir = Path("data") / "eln_attachments" / str(entry.id)
                                            if attachments_dir.exists():
                                                shutil.rmtree(attachments_dir)

                                        # Delete entry
                                        db.delete(entry)
                                        db.commit()
                                        st.success("Entry deleted successfully!")
                                        st.rerun()
                                    except Exception as e:
                                        st.error(f"Error deleting entry: {e}")
                                        db.rollback()
                                else:
                                    st.session_state[f"confirm_delete_{entry.id}"] = True
                                    st.warning("Click again to confirm deletion")
                                    st.rerun()

                            # Edit mode
                            if st.session_state.get(f"editing_entry_{entry.id}", False):
                                st.markdown("---")
                                st.markdown("### Edit Entry")

                                new_title = st.text_input("Title", value=entry.title, key=f"edit_title_{entry.id}")
                                new_content = st.text_area(
                                    "Content", value=entry.content, height=200, key=f"edit_content_{entry.id}"
                                )
                                new_type = st.selectbox(
                                    "Type",
                                    ["experiment", "observation", "protocol", "note", "analysis", "other"],
                                    index=(
                                        ["experiment", "observation", "protocol", "note", "analysis", "other"].index(
                                            entry.entry_type
                                        )
                                        if entry.entry_type
                                        else 0
                                    ),
                                    key=f"edit_type_{entry.id}",
                                )
                                new_tags = st.text_input(
                                    "Tags (comma-separated)",
                                    value=", ".join(entry.tags) if entry.tags else "",
                                    key=f"edit_tags_{entry.id}",
                                )

                                col_save, col_cancel = st.columns(2)
                                with col_save:
                                    if st.button("💾 Save Changes", key=f"save_edit_{entry.id}", type="primary"):
                                        try:
                                            entry.title = new_title
                                            entry.content = new_content
                                            entry.entry_type = new_type
                                            entry.tags = (
                                                [tag.strip() for tag in new_tags.split(",") if tag.strip()]
                                                if new_tags
                                                else None
                                            )
                                            entry.updated_at = datetime.utcnow()
                                            db.commit()
                                            st.session_state[f"editing_entry_{entry.id}"] = False
                                            st.success("Entry updated successfully!")
                                            st.rerun()
                                        except Exception as e:
                                            st.error(f"Error updating entry: {e}")
                                            db.rollback()

                                with col_cancel:
                                    if st.button("❌ Cancel", key=f"cancel_edit_{entry.id}"):
                                        st.session_state[f"editing_entry_{entry.id}"] = False
                                        st.rerun()
            else:
                st.info("No entries found. Create your first entry in the 'New Entry' tab.")

    # Tab 3: Search
    with tab3:
        st.subheader("Search Lab Notebook Entries")

        search_query = st.text_input(
            "Search", placeholder="Enter search term...", help="Search in titles, content, and tags"
        )

        if search_query:
            with db_session() as db:
                # Search in title, content, and tags
                entries = (
                    db.query(LabNotebookEntry)
                    .filter(
                        (LabNotebookEntry.title.ilike(f"%{search_query}%"))
                        | (LabNotebookEntry.content.ilike(f"%{search_query}%"))
                    )
                    .order_by(LabNotebookEntry.created_at.desc())
                    .all()
                )

                # Also check tags
                all_entries = db.query(LabNotebookEntry).all()
                tag_matches = [
                    e for e in all_entries if e.tags and any(search_query.lower() in tag.lower() for tag in e.tags)
                ]

                # Combine and deduplicate
                result_ids = {e.id for e in entries}
                for e in tag_matches:
                    if e.id not in result_ids:
                        entries.append(e)

                st.success(f"Found {len(entries)} matching entry/entries")

                if entries:
                    for entry in entries:
                        with st.expander(f"**{entry.title}** - {entry.created_at.strftime('%Y-%m-%d')}"):
                            st.write(f"**Type:** {entry.entry_type or 'N/A'}")
                            if entry.tags:
                                st.write(f"**Tags:** {', '.join(entry.tags)}")
                            st.markdown("---")
                            # Highlight search term in content (simple version)
                            content = entry.content
                            if search_query.lower() in content.lower():
                                # Simple highlighting - just show excerpt
                                idx = content.lower().find(search_query.lower())
                                start = max(0, idx - 100)
                                end = min(len(content), idx + len(search_query) + 100)
                                excerpt = content[start:end]
                                if start > 0:
                                    excerpt = "..." + excerpt
                                if end < len(content):
                                    excerpt = excerpt + "..."
                                st.markdown(excerpt)
                            else:
                                st.markdown(content[:500] + "..." if len(content) > 500 else content)
                else:
                    st.info("No entries found matching your search.")
        else:
            st.info("Enter a search term to find entries.")

    # Tab 4: Templates
    with tab4:
        st.subheader("Entry Templates")
        st.markdown("Use templates to quickly create common types of lab notebook entries.")

        template_type = st.selectbox(
            "Select Template",
            [
                "Experiment Protocol",
                "Observation Note",
                "Analysis Summary",
                "Data Collection",
                "Custom Template",
            ],
        )

        # Template definitions
        templates = {
            "Experiment Protocol": {
                "title": "Experiment Protocol: ",
                "content": """## Objective
Describe the objective of the experiment.

## Materials
- Material 1
- Material 2

## Methods
1. Step 1
2. Step 2
3. Step 3

## Expected Results
Describe expected outcomes.

## Notes
Additional notes or observations.""",
                "type": "protocol",
                "tags": "protocol, experiment",
            },
            "Observation Note": {
                "title": "Observation: ",
                "content": """## Observation
Describe what you observed.

## Context
Provide context for the observation.

## Significance
Why is this observation important?

## Follow-up
What should be done next?""",
                "type": "observation",
                "tags": "observation, note",
            },
            "Analysis Summary": {
                "title": "Analysis Summary: ",
                "content": """## Analysis Overview
Brief overview of the analysis performed.

## Methods
- Method 1
- Method 2

## Results
Key findings from the analysis.

## Conclusions
Main conclusions drawn from the analysis.

## Next Steps
Recommended follow-up actions.""",
                "type": "analysis",
                "tags": "analysis, summary",
            },
            "Data Collection": {
                "title": "Data Collection: ",
                "content": """## Data Collection Session
Date: [Date]
Location: [Location]

## Data Collected
- Dataset 1
- Dataset 2

## Collection Method
Describe how data was collected.

## Quality Checks
- Check 1: [Status]
- Check 2: [Status]

## Notes
Additional notes about the data collection.""",
                "type": "experiment",
                "tags": "data, collection",
            },
            "Custom Template": {
                "title": "",
                "content": "",
                "type": "note",
                "tags": "",
            },
        }

        template = templates.get(template_type, templates["Custom Template"])

        # Pre-fill form with template
        st.markdown("### Fill in the template:")

        template_title = st.text_input(
            "Title",
            value=template["title"],
            help="Complete the title",
        )

        template_content = st.text_area(
            "Content",
            value=template["content"],
            height=400,
            help="Edit the template content as needed",
        )

        template_entry_type = st.selectbox(
            "Entry Type",
            ["experiment", "observation", "protocol", "note", "analysis", "other"],
            index=["experiment", "observation", "protocol", "note", "analysis", "other"].index(template["type"]),
        )

        template_tags = st.text_input(
            "Tags (comma-separated)",
            value=template["tags"],
        )

        # Link to entities (same as new entry)
        st.markdown("### Link to Entities (Optional)")
        with db_session() as db:
            programs = db.query(Program).order_by(Program.name).all()
            if programs:
                selected_programs = st.multiselect(
                    "Link to Programs",
                    options=[(p.id, p.name) for p in programs],
                    format_func=lambda x: x[1],
                    help="Select programs this entry relates to",
                    key="template_programs",
                )
            else:
                selected_programs = []

            experiments = db.query(Experiment).order_by(Experiment.name).all()
            if experiments:
                selected_experiments = st.multiselect(
                    "Link to Experiments",
                    options=[(e.id, e.name) for e in experiments],
                    format_func=lambda x: x[1],
                    help="Select experiments this entry relates to",
                    key="template_experiments",
                )
            else:
                selected_experiments = []

            datasets = db.query(Dataset).order_by(Dataset.name).all()
            if datasets:
                selected_datasets = st.multiselect(
                    "Link to Datasets",
                    options=[(d.id, d.name) for d in datasets],
                    format_func=lambda x: x[1],
                    help="Select datasets this entry relates to",
                    key="template_datasets",
                )
            else:
                selected_datasets = []

        if st.button("💾 Create Entry from Template", type="primary"):
            if not template_title or not template_content:
                st.error("Please provide both a title and content for the entry.")
            else:
                with st.spinner("Creating entry from template..."):
                    try:
                        with db_session() as db:
                            # Parse tags
                            tags = (
                                [tag.strip() for tag in template_tags.split(",") if tag.strip()]
                                if template_tags
                                else []
                            )

                            # Create entry
                            entry = LabNotebookEntry(
                                title=template_title,
                                content=template_content,
                                entry_type=template_entry_type,
                                tags=tags if tags else None,
                            )
                            db.add(entry)
                            db.flush()

                            # Create associations
                            for program_id, _ in selected_programs:
                                assoc = LabNotebookEntryAssociation(
                                    entry_id=entry.id,
                                    entity_type="program",
                                    entity_id=program_id,
                                )
                                db.add(assoc)

                            for experiment_id, _ in selected_experiments:
                                assoc = LabNotebookEntryAssociation(
                                    entry_id=entry.id,
                                    entity_type="experiment",
                                    entity_id=experiment_id,
                                )
                                db.add(assoc)

                            for dataset_id, _ in selected_datasets:
                                assoc = LabNotebookEntryAssociation(
                                    entry_id=entry.id,
                                    entity_type="dataset",
                                    entity_id=dataset_id,
                                )
                                db.add(assoc)

                            db.commit()

                            st.success(f"✅ Entry '{template_title}' created from template!")
                            st.balloons()
                            st.rerun()
                    except Exception as e:
                        st.error(f"❌ Error creating entry: {str(e)}")
                        st.exception(e)
