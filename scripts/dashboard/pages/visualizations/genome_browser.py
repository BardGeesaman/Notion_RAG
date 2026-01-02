from __future__ import annotations

import json
import os
import requests
from typing import List, Dict

import streamlit as st
import streamlit.components.v1 as components
from scripts.dashboard.auth import require_auth


API_BASE = os.getenv("API_BASE_URL", "http://localhost:8000/api/v1")


def _get_auth_headers() -> dict:
    token = st.session_state.get("access_token", "")
    return {"Authorization": f"Bearer {token}"}


def _fetch_alignments() -> list:
    """Fetch uploaded alignment files from API."""
    try:
        response = requests.get(
            f"{API_BASE}/genomics/alignments",
            headers=_get_auth_headers(),
            params={"limit": 100, "has_index": True},  # Only indexed files work with IGV.js
            timeout=30,
        )
        if response.status_code == 200:
            return response.json()
    except Exception:
        pass
    return []


def _igv_html(reference: str, locus: str, tracks: List[Dict]) -> str:
    ref_map = {
        "hg38": {"id": "hg38", "name": "Human (hg38)"},
        "hg19": {"id": "hg19", "name": "Human (hg19)"},
        "mm10": {"id": "mm10", "name": "Mouse (mm10)"},
    }
    ref = ref_map.get(reference, ref_map["hg38"])
    tracks_json = json.dumps(tracks)
    locus_js = f"locus: {json.dumps(locus)}," if locus else ""
    return f"""
    <div id="igv-container" style="width:100%; height:700px; border:1px solid #ddd;"></div>
    <script src="https://cdn.jsdelivr.net/npm/igv@2.15.5/dist/igv.min.js"></script>
    <script>
      const options = {{
        genome: "{ref['id']}",
        {locus_js}
        tracks: {tracks_json}
      }};
      igv.createBrowser(document.getElementById("igv-container"), options);
    </script>
    """


def _build_default_tracks(ref: str) -> list:
    """Build default annotation tracks for reference genome."""
    if ref == "hg38":
        return [{
            "name": "RefSeq Genes",
            "type": "annotation",
            "format": "bed",
            "url": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.bed.gz",
            "indexURL": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.bed.gz.tbi",
            "displayMode": "EXPANDED",
        }]
    return []  # Other genomes - user can add custom tracks


def render() -> None:
    require_auth()
    st.header("IGV.js Genome Browser")
    st.caption("Interactive genome visualization with BAM/CRAM alignment support.")

    # Reference genome
    ref = st.selectbox("Reference genome", ["hg38", "hg19", "mm10"], index=0)
    gene = st.text_input("Gene / locus (e.g., BRCA1 or chr17:43,040,000-43,125,000)", value="")
    
    # Alignment file selector
    st.subheader("Alignment Tracks")
    alignments = _fetch_alignments()
    
    selected_alignments = []
    if alignments:
        alignment_options = {a["id"]: f"{a['filename']} ({a['file_format']})" for a in alignments}
        selected_ids = st.multiselect(
            "Select BAM/CRAM files",
            options=list(alignment_options.keys()),
            format_func=lambda x: alignment_options[x],
        )
        selected_alignments = [a for a in alignments if a["id"] in selected_ids]
    else:
        st.info("No indexed alignment files found. Upload BAM/CRAM with index in Alignments page.")
    
    # [P1] CRAM reference URL input
    cram_reference_url = st.text_input(
        "CRAM Reference FASTA URL (required for CRAM files)",
        placeholder="https://example.com/reference.fa",
        help="CRAM files need a reference genome FASTA to decode. Leave blank for BAM files.",
    )
    
    # Custom URL tracks (existing)
    st.subheader("Custom URL Tracks")
    bed_url = st.text_input("BED track URL (optional)")
    vcf_url = st.text_input("VCF track URL (optional)")
    
    # Build tracks
    tracks = _build_default_tracks(ref)
    
    # Add alignment tracks
    for alignment in selected_alignments:
        track = {
            "name": alignment["filename"],
            "type": "alignment",
            "format": alignment["file_format"].lower(),
            "url": f"{API_BASE}/genomics/alignments/{alignment['id']}/file",
            "indexURL": f"{API_BASE}/genomics/alignments/{alignment['id']}/index",
            "height": 300,
            "colorBy": "strand",
        }
        # Add reference for CRAM files
        if alignment["file_format"] == "CRAM" and cram_reference_url:
            track["referenceUrl"] = cram_reference_url
        tracks.append(track)
    
    # Add custom URL tracks
    if bed_url:
        tracks.append({"name": "Custom BED", "type": "annotation", "format": "bed", "url": bed_url})
    if vcf_url:
        tracks.append({"name": "Custom VCF", "type": "variant", "format": "vcf", "url": vcf_url})
    
    locus = gene.strip() if gene.strip() else ""
    
    if st.button("Load Genome Browser", type="primary"):
        html = _igv_html(ref, locus, tracks)
        components.html(html, height=750, scrolling=True)
    else:
        st.info("Configure tracks and click Load Genome Browser.")

