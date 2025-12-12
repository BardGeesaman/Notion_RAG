"""Import Data page for bulk importing entities."""
from __future__ import annotations

import json
from io import StringIO

import pandas as pd
import streamlit as st

from amprenta_rag.utils.data_import import (
    validate_import_data,
    import_experiments,
    import_compounds,
    import_samples,
    get_template,
)
from scripts.dashboard.db_session import db_session
from amprenta_rag.auth.session import get_current_user


def render_import_page() -> None:
    """Render the Import Data page."""
    st.header("📥 Import Data")
    st.markdown("Bulk import experiments, compounds, or samples from CSV or JSON files.")
    
    # Entity type selection
    entity_type = st.selectbox(
        "Entity Type",
        ["Experiments", "Compounds", "Samples"],
        key="import_entity_type"
    )
    
    entity_type_lower = entity_type.lower()
    
    # Download template button
    st.markdown("### Template")
    template_df = get_template(entity_type_lower)
    if not template_df.empty:
        template_csv = template_df.to_csv(index=False)
        st.download_button(
            label=f"📄 Download {entity_type} Template (CSV)",
            data=template_csv,
            file_name=f"{entity_type_lower}_template.csv",
            mime="text/csv",
            key=f"template_{entity_type_lower}",
        )
        st.caption(f"Template includes columns: {', '.join(template_df.columns.tolist())}")
    
    st.divider()
    
    # File uploader
    st.markdown("### Upload File")
    uploaded_file = st.file_uploader(
        f"Upload {entity_type} file",
        type=["csv", "json"],
        help=f"Upload a CSV or JSON file containing {entity_type_lower} data",
        key=f"upload_{entity_type_lower}",
    )
    
    if uploaded_file is not None:
        # Read file based on extension
        try:
            if uploaded_file.name.endswith(".csv"):
                df = pd.read_csv(uploaded_file)
            elif uploaded_file.name.endswith(".json"):
                data = json.load(uploaded_file)
                if isinstance(data, list):
                    df = pd.DataFrame(data)
                elif isinstance(data, dict):
                    # If JSON is a single object, wrap in list
                    df = pd.DataFrame([data])
                else:
                    st.error("Invalid JSON format. Expected array of objects or single object.")
                    return
            else:
                st.error("Unsupported file type")
                return
            
            st.success(f"✅ File loaded: {len(df)} rows")
            
            # Show preview
            st.markdown("### Preview (First 5 Rows)")
            st.dataframe(df.head(5), use_container_width=True, hide_index=True)
            
            # Validate
            st.markdown("### Validation")
            validation_errors = validate_import_data(df, entity_type_lower)
            
            if validation_errors:
                st.error("❌ Validation Failed")
                for error in validation_errors:
                    st.error(f"• {error}")
            else:
                st.success("✅ Validation passed")
            
            st.divider()
            
            # Import section
            st.markdown("### Import")
            
            if st.button(f"Import {entity_type}", type="primary", disabled=bool(validation_errors)):
                user = get_current_user()
                user_id = user.get("id") if user else None
                
                with st.spinner(f"Importing {entity_type_lower}..."):
                    with db_session() as db:
                        if entity_type_lower == "experiment":
                            result = import_experiments(df, db, user_id)
                        elif entity_type_lower == "compound":
                            result = import_compounds(df, db)
                        elif entity_type_lower == "sample":
                            result = import_samples(df, db, user_id)
                        else:
                            st.error(f"Unknown entity type: {entity_type_lower}")
                            return
                
                # Display results
                st.markdown("### Import Results")
                
                if result.get("created", 0) > 0:
                    st.success(f"✅ Successfully imported {result['created']} {entity_type_lower}(s)")
                
                if result.get("duplicates", 0) > 0:
                    st.info(f"ℹ️ Skipped {result['duplicates']} duplicate compound(s)")
                
                if result.get("errors"):
                    st.error(f"❌ {len(result['errors'])} error(s) occurred:")
                    for error in result["errors"]:
                        st.error(f"• {error}")
                
                if result.get("created", 0) == 0 and not result.get("errors"):
                    st.warning("⚠️ No records were imported. Check your data and try again.")
        except Exception as e:
            st.error(f"Error processing file: {e}")
    
    else:
        st.info("👆 Upload a file to begin importing data.")
        st.markdown("""
        **Instructions:**
        1. Download the template CSV file above
        2. Fill in your data following the template format
        3. Upload the completed file
        4. Review the preview and validation results
        5. Click "Import" to add the data to the platform
        """)


if __name__ == "__main__":
    render_import_page()
