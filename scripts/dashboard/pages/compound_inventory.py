"""Compound Inventory Dashboard - Physical compound sample management."""

import streamlit as st
import pandas as pd
from datetime import date, datetime, timedelta
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import StorageLocation, User
from amprenta_rag.models.chemistry import Compound
from amprenta_rag.services import inventory as service


def get_compounds():
    """Fetch all compounds for dropdown."""
    with db_session() as db:
        compounds = db.query(Compound).order_by(Compound.compound_id).limit(500).all()
        return [(str(c.id), c.compound_id) for c in compounds]


def get_storage_locations():
    """Fetch storage locations for dropdown."""
    with db_session() as db:
        locations = db.query(StorageLocation).order_by(StorageLocation.name).all()
        return [(str(loc.id), loc.name) for loc in locations]


def get_users():
    """Fetch users for dropdown."""
    with db_session() as db:
        users = db.query(User).order_by(User.email).limit(100).all()
        return [(str(u.id), u.email) for u in users]


def main():
    """Main compound inventory dashboard."""
    st.title("🧪 Compound Inventory")
    st.markdown("Physical compound sample management with barcode tracking and request workflow.")

    # Initialize session state
    if "selected_sample_id" not in st.session_state:
        st.session_state.selected_sample_id = None

    # ============================================================================
    # TABS
    # ============================================================================

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📦 Compound Stocks",
        "🔬 Plates",
        "📋 Requests",
        "➕ Register Sample",
        "🔍 Barcode Lookup"
    ])

    # ============================================================================
    # TAB 1: COMPOUND STOCKS
    # ============================================================================

    with tab1:
        st.header("Compound Sample Inventory")
        
        # Filters
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            status_filter = st.selectbox(
                "Status",
                ["All", "available", "depleted", "archived", "reserved"],
                key="stock_status"
            )
        
        with col2:
            # Low Stock alert filter
            low_stock_filter = st.checkbox("⚠️ Low Stock Only (<100µL)", key="low_stock")
        
        with col3:
            # Expiring filter
            expiring_filter = st.checkbox("⏰ Expiring Soon (30 days)", key="expiring")
        
        with col4:
            # Expired filter
            expired_filter = st.checkbox("❌ Expired", key="expired")
        
        # Fetch samples based on filters
        with db_session() as db:
            if low_stock_filter:
                samples = service.get_low_stock_samples(db, threshold=100.0, limit=200)
                st.warning(f"⚠️ Showing {len(samples)} low stock samples")
            elif expiring_filter:
                samples = service.get_expiring_samples(db, days=30, limit=200)
                st.warning(f"⏰ Showing {len(samples)} samples expiring within 30 days")
            elif expired_filter:
                samples = service.get_expired_samples(db, limit=200)
                st.error(f"❌ Showing {len(samples)} expired samples")
            else:
                status_val = None if status_filter == "All" else status_filter
                samples = service.list_compound_samples(db, status=status_val, limit=200)
            
            if samples:
                # Build dataframe
                data = []
                for s in samples:
                    # Visual indicators for alerts
                    alerts = []
                    if s.quantity and s.quantity < 100 and s.status == "available":
                        alerts.append("⚠️ Low")
                    if s.expiry_date:
                        if s.expiry_date < date.today():
                            alerts.append("❌ Expired")
                        elif s.expiry_date <= date.today() + timedelta(days=30):
                            alerts.append("⏰ Expiring")
                    
                    data.append({
                        "Barcode": s.barcode or "-",
                        "Compound": s.compound.compound_id if s.compound else "-",
                        "Quantity": f"{s.quantity:.1f} {s.unit}" if s.quantity else "-",
                        "Concentration": f"{s.concentration} {s.concentration_unit}" if s.concentration else "-",
                        "Location": s.storage_location.name if s.storage_location else "-",
                        "Position": s.position or "-",
                        "Status": s.status,
                        "Expiry": str(s.expiry_date) if s.expiry_date else "-",
                        "Alerts": " ".join(alerts) if alerts else "✓",
                    })
                
                df = pd.DataFrame(data)
                
                # Color code rows based on status
                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "Alerts": st.column_config.TextColumn("Alerts", width="small"),
                    }
                )
                
                st.caption(f"Showing {len(samples)} samples")
            else:
                st.info("No compound samples found. Register samples in the 'Register Sample' tab.")

    # ============================================================================
    # TAB 2: PLATES
    # ============================================================================

    with tab2:
        st.header("Compound Plates")
        
        col1, col2 = st.columns([3, 1])
        
        with col1:
            plate_status = st.selectbox(
                "Status",
                ["All", "active", "archived", "empty"],
                key="plate_status"
            )
        
        with col2:
            if st.button("➕ New Plate", type="primary"):
                st.session_state.show_plate_form = True
        
        # New Plate Form
        if st.session_state.get("show_plate_form", False):
            with st.form("new_plate_form"):
                st.subheader("Register New Plate")
                
                plate_barcode = st.text_input("Plate Barcode *", placeholder="PLATE-2025-001")
                plate_format = st.selectbox("Format *", ["96", "384", "1536"])
                plate_type = st.selectbox("Type *", ["mother", "daughter", "screening"])
                
                locations = get_storage_locations()
                loc_options = [("", "-- Select Location --")] + locations
                plate_location = st.selectbox(
                    "Storage Location",
                    options=[l[0] for l in loc_options],
                    format_func=lambda x: dict(loc_options).get(x, x)
                )
                
                col1, col2 = st.columns(2)
                with col1:
                    if st.form_submit_button("Create Plate", type="primary"):
                        if plate_barcode:
                            try:
                                with db_session() as db:
                                    service.create_compound_plate(
                                        db=db,
                                        barcode=plate_barcode,
                                        plate_format=plate_format,
                                        plate_type=plate_type,
                                        storage_location_id=UUID(plate_location) if plate_location else None,
                                    )
                                st.success(f"✅ Plate {plate_barcode} created!")
                                st.session_state.show_plate_form = False
                                st.rerun()
                            except ValueError as e:
                                st.error(str(e))
                        else:
                            st.error("Barcode is required")
                with col2:
                    if st.form_submit_button("Cancel"):
                        st.session_state.show_plate_form = False
                        st.rerun()
        
        # Plates list
        with db_session() as db:
            status_val = None if plate_status == "All" else plate_status
            plates = service.list_plates(db, status=status_val, limit=100)
            
            if plates:
                plate_data = []
                for p in plates:
                    # Count wells
                    well_count = len(service.get_plate_contents(db, p.id))
                    plate_data.append({
                        "Barcode": p.barcode,
                        "Format": p.plate_format,
                        "Type": p.plate_type,
                        "Wells Used": well_count,
                        "Location": p.storage_location.name if p.storage_location else "-",
                        "Status": p.status,
                        "Created": p.created_at.strftime("%Y-%m-%d") if p.created_at else "-",
                    })
                
                st.dataframe(pd.DataFrame(plate_data), use_container_width=True, hide_index=True)
            else:
                st.info("No plates found.")

    # ============================================================================
    # TAB 3: REQUESTS
    # ============================================================================

    with tab3:
        st.header("Compound Requests")
        
        col1, col2 = st.columns([3, 1])
        
        with col1:
            request_view = st.radio(
                "View",
                ["Pending", "All", "My Requests"],
                horizontal=True,
                key="request_view"
            )
        
        with col2:
            if st.button("➕ New Request", type="primary"):
                st.session_state.show_request_form = True
        
        # New Request Form
        if st.session_state.get("show_request_form", False):
            with st.form("new_request_form"):
                st.subheader("Create Compound Request")
                
                compounds = get_compounds()
                comp_options = [("", "-- Select Compound --")] + compounds
                req_compound = st.selectbox(
                    "Compound *",
                    options=[c[0] for c in comp_options],
                    format_func=lambda x: dict(comp_options).get(x, x)
                )
                
                req_quantity = st.number_input("Quantity (µL) *", min_value=1.0, value=50.0)
                req_purpose = st.selectbox("Purpose", ["screening", "assay", "synthesis", "other"])
                req_priority = st.selectbox("Priority", ["normal", "low", "high", "urgent"])
                req_notes = st.text_area("Notes")
                
                col1, col2 = st.columns(2)
                with col1:
                    if st.form_submit_button("Submit Request", type="primary"):
                        if req_compound:
                            try:
                                with db_session() as db:
                                    # Get first user as requester (in real app, use current_user)
                                    users = get_users()
                                    requester_id = UUID(users[0][0]) if users else None
                                    
                                    service.create_request(
                                        db=db,
                                        requester_id=requester_id,
                                        compound_id=UUID(req_compound),
                                        requested_quantity=req_quantity,
                                        quantity_unit="µL",
                                        purpose=req_purpose,
                                        priority=req_priority,
                                        notes=req_notes if req_notes else None,
                                    )
                                st.success("✅ Request submitted!")
                                st.session_state.show_request_form = False
                                st.rerun()
                            except ValueError as e:
                                st.error(str(e))
                        else:
                            st.error("Compound is required")
                with col2:
                    if st.form_submit_button("Cancel"):
                        st.session_state.show_request_form = False
                        st.rerun()
        
        # Requests list
        with db_session() as db:
            if request_view == "Pending":
                requests = service.list_pending_requests(db, limit=100)
            else:
                requests = service.list_requests(db, limit=100)
            
            if requests:
                req_data = []
                for r in requests:
                    req_data.append({
                        "ID": str(r.id)[:8],
                        "Compound": r.compound.compound_id if r.compound else "-",
                        "Quantity": f"{r.requested_quantity} {r.quantity_unit}",
                        "Purpose": r.purpose or "-",
                        "Priority": r.priority,
                        "Status": r.status,
                        "Requester": r.requester.email if r.requester else "-",
                        "Requested": r.requested_at.strftime("%Y-%m-%d %H:%M") if r.requested_at else "-",
                    })
                
                df = pd.DataFrame(req_data)
                
                # Priority color coding would go here
                st.dataframe(df, use_container_width=True, hide_index=True)
                
                # Action buttons for pending requests
                if request_view == "Pending" and requests:
                    st.subheader("Actions")
                    selected_request = st.selectbox(
                        "Select request to action",
                        options=[str(r.id) for r in requests],
                        format_func=lambda x: f"{x[:8]}... - {next((r.compound.compound_id if r.compound else 'Unknown') for r in requests if str(r.id) == x)}"
                    )
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        if st.button("✅ Approve", type="primary"):
                            try:
                                with db_session() as db:
                                    users = get_users()
                                    approver_id = UUID(users[0][0]) if users else None
                                    service.approve_request(db, UUID(selected_request), approver_id)
                                st.success("Request approved!")
                                st.rerun()
                            except ValueError as e:
                                st.error(str(e))
                    
                    with col2:
                        if st.button("📦 Fulfill"):
                            try:
                                with db_session() as db:
                                    users = get_users()
                                    fulfiller_id = UUID(users[0][0]) if users else None
                                    service.fulfill_request(db, UUID(selected_request), fulfiller_id)
                                st.success("Request fulfilled!")
                                st.rerun()
                            except ValueError as e:
                                st.error(str(e))
                    
                    with col3:
                        reject_reason = st.text_input("Rejection reason", key="reject_reason")
                        if st.button("❌ Reject"):
                            if reject_reason:
                                try:
                                    with db_session() as db:
                                        users = get_users()
                                        rejector_id = UUID(users[0][0]) if users else None
                                        service.reject_request(db, UUID(selected_request), rejector_id, reject_reason)
                                    st.success("Request rejected")
                                    st.rerun()
                                except ValueError as e:
                                    st.error(str(e))
                            else:
                                st.warning("Please provide a rejection reason")
            else:
                st.info("No requests found.")

    # ============================================================================
    # TAB 4: REGISTER SAMPLE
    # ============================================================================

    with tab4:
        st.header("Register Compound Sample")
        
        with st.form("register_sample_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Compound & Quantity")
                
                compounds = get_compounds()
                comp_options = [("", "-- Select Compound --")] + compounds
                sample_compound = st.selectbox(
                    "Compound *",
                    options=[c[0] for c in comp_options],
                    format_func=lambda x: dict(comp_options).get(x, x),
                    key="reg_compound"
                )
                
                sample_quantity = st.number_input("Quantity (µL) *", min_value=0.1, value=100.0)
                sample_concentration = st.number_input("Concentration (mM)", min_value=0.0, value=10.0)
                sample_solvent = st.selectbox("Solvent", ["DMSO", "Water", "PBS", "Methanol", "Other"])
                sample_format = st.selectbox("Format", ["tube", "vial", "plate_well"])
            
            with col2:
                st.subheader("Storage & Tracking")
                
                sample_barcode = st.text_input("Barcode (auto-generated if empty)")
                sample_batch = st.text_input("Batch/Lot Number")
                sample_expiry = st.date_input("Expiry Date", value=date.today() + timedelta(days=365))
                
                locations = get_storage_locations()
                loc_options = [("", "-- Select Location --")] + locations
                sample_location = st.selectbox(
                    "Storage Location",
                    options=[l[0] for l in loc_options],
                    format_func=lambda x: dict(loc_options).get(x, x),
                    key="reg_location"
                )
                
                sample_position = st.text_input("Position (e.g., A01, Shelf-3)")
            
            sample_notes = st.text_area("Notes")
            
            if st.form_submit_button("Register Sample", type="primary", use_container_width=True):
                if sample_compound:
                    try:
                        with db_session() as db:
                            sample = service.create_compound_sample(
                                db=db,
                                compound_id=UUID(sample_compound),
                                quantity=sample_quantity,
                                quantity_unit="µL",
                                concentration=sample_concentration if sample_concentration > 0 else None,
                                concentration_unit="mM" if sample_concentration > 0 else None,
                                solvent=sample_solvent,
                                format=sample_format,
                                batch_lot=sample_batch if sample_batch else None,
                                expiry_date=sample_expiry,
                                storage_location_id=UUID(sample_location) if sample_location else None,
                                position=sample_position if sample_position else None,
                                barcode=sample_barcode if sample_barcode else None,
                                notes=sample_notes if sample_notes else None,
                            )
                        st.success(f"✅ Sample registered with barcode: **{sample.barcode}**")
                        st.balloons()
                    except ValueError as e:
                        st.error(f"Error: {e}")
                else:
                    st.error("Please select a compound")

    # ============================================================================
    # TAB 5: BARCODE LOOKUP
    # ============================================================================

    with tab5:
        st.header("🔍 Barcode Lookup")
        
        barcode_input = st.text_input(
            "Scan or enter barcode",
            placeholder="COMP-20250103-ABC123",
            key="barcode_lookup"
        )
        
        if barcode_input:
            with db_session() as db:
                result = service.lookup_barcode(db, barcode_input)
                
                if result["type"] == "sample":
                    sample = result["entity"]
                    st.success("✅ Found: Compound Sample")
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Barcode", sample.barcode)
                        st.metric("Compound", sample.compound.compound_id if sample.compound else "-")
                        st.metric("Quantity", f"{sample.quantity} {sample.unit}" if sample.quantity else "-")
                        st.metric("Status", sample.status)
                    
                    with col2:
                        st.metric("Concentration", f"{sample.concentration} {sample.concentration_unit}" if sample.concentration else "-")
                        st.metric("Location", sample.storage_location.name if sample.storage_location else "-")
                        st.metric("Position", sample.position or "-")
                        st.metric("Expiry", str(sample.expiry_date) if sample.expiry_date else "-")
                    
                    # Quick actions
                    st.subheader("Quick Actions")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        dispense_qty = st.number_input("Dispense (µL)", min_value=1.0, value=10.0, key="dispense_qty")
                        if st.button("💧 Dispense"):
                            try:
                                with db_session() as db:
                                    service.update_compound_sample(
                                        db, sample.id, 
                                        quantity=sample.quantity - dispense_qty
                                    )
                                st.success(f"Dispensed {dispense_qty} µL")
                                st.rerun()
                            except ValueError as e:
                                st.error(str(e))
                    
                    with col2:
                        if st.button("🗑️ Mark Depleted"):
                            try:
                                with db_session() as db:
                                    service.deplete_compound_sample(db, sample.id)
                                st.success("Sample marked as depleted")
                                st.rerun()
                            except ValueError as e:
                                st.error(str(e))
                
                elif result["type"] == "plate":
                    plate = result["entity"]
                    st.success("✅ Found: Compound Plate")
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Barcode", plate.barcode)
                        st.metric("Format", plate.plate_format)
                        st.metric("Type", plate.plate_type)
                    
                    with col2:
                        st.metric("Status", plate.status)
                        st.metric("Location", plate.storage_location.name if plate.storage_location else "-")
                    
                    # Show plate contents
                    contents = service.get_plate_contents(db, plate.id)
                    if contents:
                        st.subheader(f"Plate Contents ({len(contents)} wells)")
                        well_data = []
                        for s in contents:
                            well_data.append({
                                "Well": s.well_position or "-",
                                "Compound": s.compound.compound_id if s.compound else "-",
                                "Quantity": f"{s.quantity} {s.unit}" if s.quantity else "-",
                                "Status": s.status,
                            })
                        st.dataframe(pd.DataFrame(well_data), use_container_width=True, hide_index=True)
                    else:
                        st.info("Plate is empty")
                
                else:
                    st.warning(f"⚠️ No sample or plate found with barcode: {barcode_input}")
        else:
            st.info("Enter a barcode above to look up a sample or plate")
            
            # Show recent barcodes for convenience
            st.subheader("Recent Samples")
            with db_session() as db:
                recent = service.list_compound_samples(db, limit=10)
                if recent:
                    for s in recent:
                        st.code(s.barcode, language=None)


if __name__ == "__main__":
    main()
