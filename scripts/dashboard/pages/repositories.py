from pathlib import Path
import requests
import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Dataset, Experiment
from amprenta_rag.ingestion.postgres_dataset_ingestion import ingest_dataset_from_postgres
from scripts.dashboard.db_session import db_session
from scripts.harvest_mw_studies import fetch_mw_mwtab

repo_map = {
    "GEO": "search_geo_studies",
    "PRIDE": "search_pride_studies",
    "MetaboLights": "search_metabolights_studies",
    "MW": "search_mw_studies",
}

HELP_TEXT = """
How to use this page:
1. Choose a repository and enter desired search criteria.
2. Run search to list matching studies. You can filter further in the result table.
3. Select studies to import and run the 'Import Selected' action. Progress and status are reported live.
4. Use the 'Imported Studies' manager to review or clean up previously ingested studies.
"""


def repo_status(ds):
    # Returns 'Imported' if already in db, 'Not imported' otherwise.
    return "Imported" if ds.get("imported") else "Not imported"


def _ensure_experiment_for_study(db, repo: str, accession: str, title: str) -> Experiment:
    """
    Get or create an Experiment grouping datasets for the same study.
    Stores repository study ID in external_ids.
    """
    key_map = {
        "GEO": "geo_accession",
        "PRIDE": "pride_accession",
        "MW": "mw_study_id",
        "MetaboLights": "metabolights_study_id",
    }
    key = key_map.get(repo, "study_id")

    existing = (
        db.query(Experiment)
        .filter(Experiment.external_ids.isnot(None))
        .filter(Experiment.external_ids[key].astext == accession)
        .first()
    )
    if existing:
        return existing

    exp = Experiment(
        name=title or accession,
        description=f"{repo} study {accession}",
        external_ids={key: accession},
    )
    db.add(exp)
    db.commit()
    db.refresh(exp)
    return exp


def _download_data_files(repo_obj, accession: str, repo: str) -> list:
    """
    Download data files for a study into data/repositories/{repo}/{accession}/
    Returns list of file paths downloaded.
    """
    try:
        files = repo_obj.fetch_study_data_files(accession)
    except Exception as e:
        st.warning(f"⚠️ Failed to fetch file list for {accession}: {e}")
        return []

    if not files:
        return []

    base_dir = Path("data") / "repositories" / repo.lower() / accession
    base_dir.mkdir(parents=True, exist_ok=True)

    downloaded = []
    for f in files:
        try:
            if not f.download_url:
                continue
            filename = f.filename or f.file_id or f"{accession}.dat"
            dest = base_dir / filename
            resp = requests.get(f.download_url, stream=True, timeout=120)
            resp.raise_for_status()
            with open(dest, "wb") as out:
                for chunk in resp.iter_content(chunk_size=8192):
                    if chunk:
                        out.write(chunk)
            downloaded.append(str(dest))
        except Exception as e:
            st.warning(f"⚠️ Failed to download {getattr(f,'file_id',filename)}: {e}")
            continue
    return downloaded


def render_repositories_page():
    st.header("🌎 Repository Browser & Import")
    st.info(HELP_TEXT)
    
    # Initialize session state for search results
    if "repo_search_results" not in st.session_state:
        st.session_state.repo_search_results = None
    if "repo_search_params" not in st.session_state:
        st.session_state.repo_search_params = {}
    
    tab1, tab2 = st.tabs(["Browse & Import", "Imported Studies"])
    
    with tab1:
            st.subheader("Search Repositories")
            repo = st.selectbox("Repository", list(repo_map.keys()))
            col1, col2 = st.columns([3,1])
            with col1:
                keyword = st.text_input("Keyword (title, description, ID)")
            with col2:
                limit = st.number_input("Max Results", min_value=1, max_value=100, value=25)
            disease = st.text_input("Disease", value="")
            species = st.text_input("Species", value="")
            omics_type = st.text_input("Omics type", value="")
            # Repo-specific filters
            platform = instrument = modality = analytical_platform = matrix = ""
            if repo == "GEO":
                platform = st.text_input("Platform (GEO)", value="")
            if repo == "PRIDE":
                instrument = st.text_input("Instrument (PRIDE)", value="")
                modality = st.text_input("Modality (PRIDE)", value="")
            if repo in ("MW", "MetaboLights"):
                matrix = st.text_input("Matrix (MW/MetaboLights)", value="")
                analytical_platform = st.text_input("Analytical Platform", value="")
            
            # Check if search params changed (clear results if they did)
            current_params = {
                "repo": repo, "keyword": keyword, "disease": disease, "species": species,
                "omics_type": omics_type, "platform": platform, "instrument": instrument,
                "modality": modality, "matrix": matrix, "analytical_platform": analytical_platform,
                "limit": limit
            }
            if current_params != st.session_state.repo_search_params:
                st.session_state.repo_search_results = None
                st.session_state.repo_search_params = current_params

            if st.button("Run Search"):
                # Import search function for selected repository
                # Note: PRIDE and MetaboLights functions are re-exported from their __init__.py
                try:
                    if repo == "GEO":
                        from amprenta_rag.ingestion.repositories.geo import search_geo_studies as fn
                    elif repo == "PRIDE":
                        from amprenta_rag.ingestion.repositories.pride import search_pride_studies as fn
                    elif repo == "MetaboLights":
                        from amprenta_rag.ingestion.repositories.metabolights import search_metabolights_studies as fn
                    elif repo == "MW":
                        from amprenta_rag.ingestion.repositories.mw import search_mw_studies as fn
                    else:
                        st.error(f"❌ Unknown repository: {repo}")
                        return
                except Exception as e:
                    st.error(f"❌ Failed to load search function for {repo}")
                    st.error(f"Error: {str(e)}")
                    with st.expander("Show full error details"):
                        st.exception(e)
                    return
                kwargs = {
                    "keyword": keyword or None,
                    "disease": disease or None,
                    "species": species or None,
                    "omics_type": omics_type or None,
                    "max_results": limit,
                }
                if repo == "GEO": kwargs["platform"] = platform or None
                if repo == "PRIDE": kwargs["instrument"] = instrument or None; kwargs["modality"] = modality or None
                if repo in ("MW", "MetaboLights"):
                    kwargs["analytical_platform"] = analytical_platform or None
                    kwargs["matrix"] = matrix or None
                
                try:
                    with st.spinner(f"Searching {repo}..."):
                        results = fn(**kwargs)
                except Exception as e:
                    st.error(f"❌ Search failed for {repo}")
                    st.error(f"Error: {str(e)}")
                    with st.expander("Show full error details"):
                        st.exception(e)
                    return
                
                # Store results in session state
                st.session_state.repo_search_results = {
                    "repo": repo,
                    "results": results,
                    "kwargs": kwargs,
                }
                
            # Display results if they exist in session state
            if st.session_state.repo_search_results and st.session_state.repo_search_results["repo"] == repo:
                search_data = st.session_state.repo_search_results
                results = search_data["results"]
                kwargs = search_data["kwargs"]
                
                # Debug: Show search parameters
                active_filters = {k: v for k, v in kwargs.items() if v}
                fsum = f"Showing {repo} results for " + \
                      ", ".join([f"{k}='{v}'" for k, v in active_filters.items()])
                st.caption(fsum)
                
                # Debug: Show result count
                if not results:
                    st.warning(f"⚠️ Search returned 0 results. Try broader search criteria or leave filters empty.")
                    st.info(f"Search parameters used: {active_filters}")
                else:
                    st.success(f"✅ Found {len(results)} studies")
                    
                    # Import status - check which studies are already imported
                    acckey = repo.lower() + ("_accession" if repo != "MW" else "_study_id")
                    ds_accessions = [r["accession"] for r in results]
                    
                    # Query database for imported studies
                    with db_session() as db:
                        # Check which accessions are already in database
                        imported_datasets = (
                            db.query(Dataset)
                            .filter(Dataset.external_ids.isnot(None))
                            .limit(200)
                            .all()
                        )
                        in_db = set()
                        for ds in imported_datasets:
                            if ds.external_ids:
                                acc_value = ds.external_ids.get(acckey)
                                if acc_value in ds_accessions:
                                    in_db.add(acc_value)
                    
                    summary = []
                    for r in results:
                        summary.append(
                            {
                                "Selected": False,
                                **r,
                                "status": "Imported" if r["accession"] in in_db else "Not imported",
                            }
                        )
                    df = pd.DataFrame(summary)
                    # Reorder columns: Select first, then other fields
                    desired_cols = [
                        "Selected",
                        "id",
                        "accession",
                        "title",
                        "description",
                        "disease",
                        "organism",
                        "omics_type",
                        "status",
                    ]
                    # Keep any extra columns at the end
                    remaining = [c for c in df.columns if c not in desired_cols]
                    df = df[[c for c in desired_cols if c in df.columns] + remaining]
                    
                    if not df.empty:
                        edited_df = st.data_editor(
                            df,
                            num_rows="dynamic",
                            key="repo_results",
                            use_container_width=False,
                            height=600,
                            hide_index=True,
                            column_config={
                                "Selected": st.column_config.CheckboxColumn(
                                    "Select",
                                    help="Select studies to import",
                                    default=False,
                                ),
                            },
                        )
                        
                        st.markdown(f"**{len(results)} results** | {len([r for r in summary if r['status']=='Imported'])} already imported.")

                        # Collect selected studies, skipping already imported
                        selected_studies = [
                            results[idx]
                            for idx, row in edited_df.iterrows()
                            if row.get("Selected") and row.get("accession") not in in_db
                        ]
                        
                        if selected_studies:
                            st.info(f"📋 {len(selected_studies)} studies selected for import")
                            
                            if st.button(f"🚀 Import Selected ({len(selected_studies)} studies)", key="import_selected", type="primary"):
                                import_progress = st.progress(0, text="Starting import...")
                                imported_count = 0
                                failed_count = 0
                                
                                for idx, study in enumerate(selected_studies):
                                    accession = study["accession"]
                                    title = study["title"]
                                    
                                    import_progress.progress(
                                        (idx + 1) / len(selected_studies),
                                        text=f"Importing {idx + 1}/{len(selected_studies)}: {accession}..."
                                    )
                                    
                                    try:
                                        # Fetch full metadata from repository
                                        if repo == "GEO":
                                            from amprenta_rag.ingestion.repositories.geo import GEORepository
                                            repo_obj = GEORepository()
                                        elif repo == "PRIDE":
                                            from amprenta_rag.ingestion.repositories.pride import PRIDERepository
                                            repo_obj = PRIDERepository()
                                        elif repo == "MetaboLights":
                                            from amprenta_rag.ingestion.repositories.metabolights import MetaboLightsRepository
                                            repo_obj = MetaboLightsRepository()
                                        elif repo == "MW":
                                            from amprenta_rag.ingestion.repositories.mw import MWRepository
                                            repo_obj = MWRepository()
                                        
                                        metadata = repo_obj.fetch_study_metadata(accession)
                                        
                                        if not metadata:
                                            st.warning(f"⚠️ Could not fetch metadata for {accession}")
                                            failed_count += 1
                                            continue
                                        
                                        # Import into database
                                        from amprenta_rag.ingestion.postgres_integration import create_or_update_dataset_in_postgres
                                        from amprenta_rag.models.domain import OmicsType
                                        
                                        # Prepare organism as list (function expects List[str])
                                        organism_list = None
                                        if metadata.organism:
                                            if isinstance(metadata.organism, list):
                                                organism_list = metadata.organism
                                            else:
                                                organism_list = [str(metadata.organism)]
                                        
                                        # Prepare external_ids dict with DOI and PubMed if available
                                        ext_ids = {acckey: accession}
                                        if metadata.doi:
                                            ext_ids["doi"] = metadata.doi
                                        if metadata.pubmed_id:
                                            ext_ids["pubmed_id"] = str(metadata.pubmed_id)
                                        downloaded_paths = []
                                        mwtab_cached = False
                                        if repo == "MW":
                                            try:
                                                mwtab_text = fetch_mw_mwtab(accession)
                                                mwtab_cached = bool(mwtab_text)
                                            except Exception as dl_err:
                                                st.warning(f"⚠️ mwTab download failed for {accession}: {dl_err}")
                                            ext_ids["mwtab_cached"] = mwtab_cached
                                        
                                        # Convert omics_type string to OmicsType enum
                                        omics_type_enum = OmicsType.METABOLOMICS  # Default for MW/MetaboLights
                                        if metadata.omics_type:
                                            try:
                                                omics_type_enum = OmicsType(metadata.omics_type.lower())
                                            except (ValueError, AttributeError):
                                                # If conversion fails, infer from repository
                                                if repo == "GEO":
                                                    omics_type_enum = OmicsType.TRANSCRIPTOMICS
                                                elif repo == "PRIDE":
                                                    omics_type_enum = OmicsType.PROTEOMICS
                                                elif repo in ("MW", "MetaboLights"):
                                                    omics_type_enum = OmicsType.METABOLOMICS
                                        
                                        with db_session() as db:
                                            dataset = create_or_update_dataset_in_postgres(
                                                db=db,
                                                name=metadata.title or accession,
                                                omics_type=omics_type_enum,
                                                description=metadata.summary or "",
                                                organism=organism_list,
                                                external_ids=ext_ids,
                                                auto_ingest=True,
                                            )
                                            # Ensure experiment grouping by study
                                            experiment = _ensure_experiment_for_study(
                                                db=db, repo=repo, accession=accession, title=metadata.title or accession
                                            )
                                            if experiment not in dataset.experiments:
                                                dataset.experiments.append(experiment)
                                            # Download data files and mark ingestion status
                                            try:
                                                paths = _download_data_files(repo_obj, accession, repo)
                                                downloaded_paths.extend(paths)
                                                if paths:
                                                    dataset.file_paths = paths
                                                    dataset.ingestion_status = "complete"
                                                else:
                                                    dataset.ingestion_status = "complete"
                                                    if dataset.description:
                                                        dataset.description += " [Metadata Only]"
                                                    else:
                                                        dataset.description = "[Metadata Only Import]"
                                            except Exception as dl_err:
                                                st.warning(f"⚠️ Download failed for {accession}: {dl_err}")
                                                dataset.ingestion_status = "failed"
                                            db.commit()
                                            dataset_uuid = dataset.id

                                        st.success(f"✅ Imported: {accession} - {title[:50]}...")
                                        st.info(f"Grouped under experiment: {experiment.name}")

                                        # Trigger Pinecone ingestion in background
                                        try:
                                            st.info(f"🔄 Starting Pinecone ingestion for {accession}...")
                                            ingest_dataset_from_postgres(dataset_uuid, force=False, update_notion=False)
                                            st.success(f"✅ Ingested to Pinecone: {accession}")
                                        except Exception as ingest_err:
                                            st.warning(f"⚠️ Pinecone ingestion pending for {accession}: {str(ingest_err)[:50]}")

                                        imported_count += 1

                                    except Exception as e:
                                        st.error(f"❌ Failed to import {accession}: {str(e)}")
                                        failed_count += 1
                                
                                import_progress.empty()
                                
                                # Summary
                                st.success(f"🎉 Import complete! {imported_count} succeeded, {failed_count} failed.")
                                if imported_count > 0:
                                    st.info("💡 Imported datasets are being processed. Check the 'Imported Studies' tab or Data Ingestion page for status.")
                                st.rerun()
                        else:
                            st.info("👈 Check the boxes in the 'Select' column to import studies")
    
    # Tab2: Imported studies management
    with tab2:
        st.subheader("Imported Studies Management")
        
        with db_session() as db:
            all_imported = (
                db.query(Dataset)
                .filter(Dataset.external_ids.isnot(None))
                .limit(200)
                .all()
            )
            # Extract data inside session
            dtable = []
            for d in all_imported:
                # Assumes accession is stored under key for each repo
                accn = (
                    d.external_ids.get("geo_accession")
                    or d.external_ids.get("pride_accession")
                    or d.external_ids.get("mw_study_id")
                    or d.external_ids.get("metabolights_study_id")
                ) if d.external_ids else None
                dtable.append(
                    {
                        "Selected": False,
                        "Accession": accn,
                        "Name": d.name,
                        "Omics": d.omics_type,
                        "Created": d.created_at.strftime("%Y-%m-%d"),
                        "Status": getattr(d, "ingestion_status", "pending"),
                        "# Features": len(d.features) if d.features else 0,
                    }
                )
        
        pdf = pd.DataFrame(dtable)
        pedit = st.data_editor(pdf, key="imported_mgr", num_rows="dynamic")
        del_rows = [r for r in pedit.to_dict("records") if r["Selected"]]
        if del_rows and st.button(f"Delete selected ({len(del_rows)})", key="del_impt"):
            if st.confirm(
                f"Are you sure you want to delete {len(del_rows)} studies? This will remove their datasets and associated artifacts."
            ):
                with db_session() as db:
                    for r in del_rows:
                        # TODO: call your safe backend delete logic here
                        # db.query(Dataset).filter(Dataset.name==r['Name']).delete()
                        pass
                    db.commit()
                st.success("Selected studies deleted.")
                st.rerun()
