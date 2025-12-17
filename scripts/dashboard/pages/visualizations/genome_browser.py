from __future__ import annotations

import json
from typing import List, Dict

import streamlit as st
import streamlit.components.v1 as components


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


def render() -> None:
    st.header("IGV.js Genome Browser")
    st.caption("Interactive genome visualization with IGV.js (demo tracks).")

    ref = st.selectbox("Reference genome", ["hg38", "hg19", "mm10"], index=0)
    gene = st.text_input("Gene / locus (e.g., BRCA1 or chr17:43,040,000-43,125,000)", value="")
    bed_url = st.text_input("BED track URL (optional)")
    vcf_url = st.text_input("VCF track URL (optional)")

    default_tracks = [
        {
            "name": "RefSeq Genes",
            "type": "annotation",
            "format": "bed",
            "url": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.bed.gz",
            "indexURL": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.bed.gz.tbi",
            "displayMode": "EXPANDED",
        }
    ]

    tracks: List[Dict] = list(default_tracks)

    if bed_url:
        tracks.append(
            {
                "name": "Custom BED",
                "type": "annotation",
                "format": "bed",
                "url": bed_url,
            }
        )
    if vcf_url:
        tracks.append(
            {
                "name": "Custom VCF",
                "type": "variant",
                "format": "vcf",
                "url": vcf_url,
                "indexURL": vcf_url + ".tbi" if not vcf_url.endswith(".tbi") else None,
            }
        )

    locus = gene.strip() if gene.strip() else ""

    if st.button("Load Genome Browser", type="primary"):
        html = _igv_html(ref, locus, tracks)
        components.html(html, height=750, scrolling=True)
    else:
        st.info("Set reference and optional tracks, then click Load Genome Browser.")

