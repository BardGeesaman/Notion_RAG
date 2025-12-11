"""Sample Inventory Management."""
import uuid
from datetime import datetime
from typing import Optional, List

import streamlit as st

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import StorageLocation, Sample, SampleTransfer, Experiment


def render_sample_inventory_page():
    st.title("🧊 Sample Inventory")

    tab1, tab2, tab3, tab4 = st.tabs([
        "Storage Locations",
        "Samples",
        "Register Sample",
        "Transfer Sample",
    ])

    with tab1:
        render_locations_tab()
    with tab2:
        render_samples_tab()
    with tab3:
        render_register_tab()
    with tab4:
        render_transfer_tab()


def _get_db():
    db_gen = get_db()
    db = next(db_gen)
    return db_gen, db


def render_locations_tab():
    st.subheader("Storage Locations")

    db_gen, db = _get_db()
    try:
        locations = db.query(StorageLocation).order_by(StorageLocation.name).all()
        parent_options = {loc.name: loc.id for loc in locations}

        with st.form("create_location"):
            name = st.text_input("Name")
            location_type = st.selectbox("Type", ["freezer", "shelf", "box", "rack", "other"], index=0)
            parent_name = st.selectbox("Parent", ["None"] + list(parent_options.keys()))
            temperature = st.text_input("Temperature", placeholder="-80C, 4C, RT")
            capacity = st.number_input("Capacity", min_value=0, step=1, value=0)
            description = st.text_area("Description", height=80)
            submitted = st.form_submit_button("Add Location", type="primary")

            if submitted:
                parent_id = parent_options.get(parent_name) if parent_name != "None" else None
                loc = StorageLocation(
                    name=name,
                    location_type=location_type,
                    parent_id=parent_id,
                    temperature=temperature or None,
                    capacity=capacity or None,
                    description=description or None,
                )
                db.add(loc)
                db.commit()
                st.success("Location created")
                st.rerun()

        if locations:
            st.markdown("### Existing Locations")
            for loc in locations:
                label = f"{loc.name} ({loc.location_type or 'n/a'})"
                if loc.parent:
                    label += f" — Parent: {loc.parent.name}"
                st.markdown(f"- {label}")
        else:
            st.info("No storage locations yet. Add one above.")
    finally:
        db_gen.close()


def render_samples_tab():
    st.subheader("Samples")

    search = st.text_input("Search", placeholder="Name or barcode")
    status_filter = st.selectbox("Status", ["All", "available", "reserved", "depleted"], index=0)

    db_gen, db = _get_db()
    try:
        query = db.query(Sample)
        if search:
            like = f"%{search}%"
            query = query.filter((Sample.name.ilike(like)) | (Sample.barcode.ilike(like)))
        if status_filter != "All":
            query = query.filter(Sample.status == status_filter)
        samples = query.order_by(Sample.created_at.desc()).limit(200).all()

        if not samples:
            st.info("No samples found.")
            return

        rows = []
        for s in samples:
            rows.append({
                "Name": s.name,
                "Barcode": s.barcode or "-",
                "Type": s.sample_type or "-",
                "Status": s.status,
                "Location": s.storage_location.name if s.storage_location else "-",
                "Position": s.position or "-",
                "Experiment": s.experiment.name if s.experiment else "-",
                "Created": s.created_at.strftime("%Y-%m-%d") if s.created_at else "",
            })
        st.dataframe(rows, use_container_width=True, hide_index=True)
    finally:
        db_gen.close()


def render_register_tab():
    st.subheader("Register Sample")

    db_gen, db = _get_db()
    try:
        locations = db.query(StorageLocation).order_by(StorageLocation.name).all()
        parent_samples = db.query(Sample).order_by(Sample.name).limit(200).all()
        experiments = db.query(Experiment).order_by(Experiment.name).limit(200).all()

        loc_options = {"(none)": None, **{loc.name: loc.id for loc in locations}}
        parent_options = {"(none)": None, **{s.name: s.id for s in parent_samples}}
        exp_options = {"(none)": None, **{e.name: e.id for e in experiments}}

        with st.form("register_sample"):
            name = st.text_input("Name*")
            sample_type = st.text_input("Sample Type")
            auto_barcode = st.checkbox("Auto-generate barcode", value=True)
            barcode = st.text_input("Barcode", disabled=auto_barcode)
            if auto_barcode:
                barcode_val = f"S-{datetime.utcnow().strftime('%Y%m%d%H%M%S')}-{uuid.uuid4().hex[:6]}"
            else:
                barcode_val = barcode

            location_name = st.selectbox("Storage Location", list(loc_options.keys()))
            position = st.text_input("Position (e.g., A1)")
            quantity = st.number_input("Quantity", min_value=0.0, step=0.1, value=0.0)
            unit = st.text_input("Unit", value="µL")
            parent_name = st.selectbox("Parent Sample", list(parent_options.keys()))
            exp_name = st.selectbox("Experiment", list(exp_options.keys()))
            notes = st.text_area("Notes", height=80)

            submitted = st.form_submit_button("Register", type="primary")
            if submitted:
                if not name:
                    st.error("Name is required")
                    return
                if not barcode_val:
                    st.error("Barcode is required")
                    return

                sample = Sample(
                    name=name,
                    sample_type=sample_type or None,
                    barcode=barcode_val,
                    storage_location_id=loc_options.get(location_name),
                    position=position or None,
                    parent_sample_id=parent_options.get(parent_name),
                    experiment_id=exp_options.get(exp_name),
                    quantity=quantity or None,
                    unit=unit or None,
                    status="available",
                    created_at=datetime.utcnow(),
                    notes=notes or None,
                )
                db.add(sample)
                db.commit()
                st.success(f"Sample '{name}' registered with barcode {barcode_val}")
                st.rerun()
    finally:
        db_gen.close()


def render_transfer_tab():
    st.subheader("Transfer Sample")

    db_gen, db = _get_db()
    try:
        samples = db.query(Sample).order_by(Sample.created_at.desc()).limit(200).all()
        locations = db.query(StorageLocation).order_by(StorageLocation.name).all()

        if not samples or not locations:
            st.info("Need at least one sample and one storage location to transfer.")
            return

        sample_options = {f"{s.name} ({s.barcode or 'no barcode'})": s.id for s in samples}
        location_options = {loc.name: loc.id for loc in locations}

        with st.form("transfer_sample"):
            sample_key = st.selectbox("Sample", list(sample_options.keys()))
            from_loc = st.selectbox("From", ["(unknown)"] + list(location_options.keys()))
            to_loc = st.selectbox("To", list(location_options.keys()))
            notes = st.text_area("Notes", height=80)
            submitted = st.form_submit_button("Transfer", type="primary")

            if submitted:
                transfer = SampleTransfer(
                    sample_id=sample_options[sample_key],
                    from_location_id=location_options.get(from_loc) if from_loc != "(unknown)" else None,
                    to_location_id=location_options[to_loc],
                    transferred_at=datetime.utcnow(),
                    notes=notes or None,
                )
                db.add(transfer)
                # Update sample location
                sample_obj = db.query(Sample).filter(Sample.id == sample_options[sample_key]).first()
                if sample_obj:
                    sample_obj.storage_location_id = location_options[to_loc]
                db.commit()
                st.success("Transfer recorded")
                st.rerun()
    finally:
        db_gen.close()


if __name__ == "__main__":
    render_sample_inventory_page()
