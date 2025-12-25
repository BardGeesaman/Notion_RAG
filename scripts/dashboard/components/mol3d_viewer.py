"""Streamlit 3D molecule viewer component (py3Dmol optional)."""

from __future__ import annotations

import uuid
from io import BytesIO
from typing import List, Optional

import streamlit as st

try:
    import py3Dmol  # type: ignore

    PY3DMOL_AVAILABLE = True
except Exception:  # noqa: BLE001
    PY3DMOL_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import Draw

    RDKIT_AVAILABLE = True
except Exception:  # noqa: BLE001
    RDKIT_AVAILABLE = False

from amprenta_rag.chemistry.conformers import (
    RDKIT_AVAILABLE as RDKIT_CONFORMERS_AVAILABLE,
    align_molecules_to_reference,
    conformer_to_pdb,
    conformer_energies,
    generate_conformers,
    optimize_conformer,
)


def _build_py3dmol_html(pdb_string: str, style: str, width: int, height: int) -> str:
    """Build a standalone HTML snippet to render a single PDB in 3Dmol.js."""
    stl = (style or "stick").lower()
    style_obj = (
        {"stick": {}}
        if stl == "stick"
        else {"sphere": {}}
        if stl == "sphere"
        else {"cartoon": {}}
        if stl == "cartoon"
        else {"line": {}}
    )

    dom_id = f"mol3d_{uuid.uuid4().hex}"
    pdb_esc = (pdb_string or "").replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$")
    html = f"""
<div id="{dom_id}" style="width:{int(width)}px;height:{int(height)}px;position:relative;"></div>
<script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
<script>
  (function() {{
    const viewer = $3Dmol.createViewer(document.getElementById("{dom_id}"), {{ backgroundColor: "0xeeeeee" }});
    const pdb = `{pdb_esc}`;
    viewer.addModel(pdb, "pdb");
    viewer.setStyle({{}}, {style_obj});
    viewer.zoomTo();
    viewer.render();
  }})();
</script>
"""
    return html


def _build_py3dmol_html_multi(pdb_strings: List[str], colors: List[str], style: str, width: int, height: int) -> str:
    stl = (style or "stick").lower()
    base_style = "stick" if stl == "stick" else "sphere" if stl == "sphere" else "line"

    dom_id = f"mol3d_multi_{uuid.uuid4().hex}"
    models = []
    for pdb in pdb_strings:
        pdb_esc = (pdb or "").replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$")
        models.append(f"`{pdb_esc}`")
    colors_js = [c if c.startswith("#") else f"#{c}" for c in colors]

    html = f"""
<div id="{dom_id}" style="width:{int(width)}px;height:{int(height)}px;position:relative;"></div>
<script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
<script>
  (function() {{
    const viewer = $3Dmol.createViewer(document.getElementById("{dom_id}"), {{ backgroundColor: "0xeeeeee" }});
    const pdbs = [{",".join(models)}];
    const colors = {colors_js};
    for (let i = 0; i < pdbs.length; i++) {{
      const m = viewer.addModel(pdbs[i], "pdb");
      const c = colors[i % colors.length] || "#1f77b4";
      const style = {{}};
      style["{base_style}"] = {{ color: c }};
      viewer.setStyle({{model: i}}, style);
    }}
    viewer.zoomTo();
    viewer.render();
  }})();
</script>
"""
    return html


def _render_2d_fallback(smiles: str) -> None:
    if not RDKIT_AVAILABLE:
        st.info("3D rendering unavailable (py3Dmol failed) and RDKit not installed for 2D fallback.")
        return
    m = Chem.MolFromSmiles(str(smiles or ""))
    if m is None:
        st.error("Invalid SMILES.")
        return
    img = Draw.MolToImage(m, size=(400, 300))
    buf = BytesIO()
    img.save(buf, format="PNG")
    st.image(buf.getvalue(), caption="2D structure (fallback)", use_container_width=False)


def render_molecule_3d(smiles: str, style: str = "stick", width: int = 400, height: int = 400) -> None:
    """Render a single molecule as lowest-energy conformer."""
    if not RDKIT_CONFORMERS_AVAILABLE:
        st.info("RDKit not available; cannot generate 3D conformers.")
        _render_2d_fallback(smiles)
        return

    try:
        confs = generate_conformers(smiles, n_conformers=5, method="ETKDG")
        energies = []
        for m in confs:
            optimize_conformer(m, force_field="MMFF")
            energies.append(conformer_energies(m, prefer="MMFF")[0])
        best_idx = int(min(range(len(energies)), key=lambda i: energies[i]))
        pdb = conformer_to_pdb(confs[best_idx], conf_id=0)
    except Exception as e:  # noqa: BLE001
        st.error(f"3D conformer generation failed: {e}")
        _render_2d_fallback(smiles)
        return

    try:
        if PY3DMOL_AVAILABLE:
            view = py3Dmol.view(width=int(width), height=int(height))
            view.addModel(pdb, "pdb")
            view.setStyle({style: {}})
            view.setBackgroundColor("0xeeeeee")
            view.zoomTo()
            html = view._make_html()  # noqa: SLF001
        else:
            html = _build_py3dmol_html(pdb, style, width, height)
        st.components.v1.html(html, height=int(height) + 20, scrolling=False)
    except Exception:
        _render_2d_fallback(smiles)


def render_conformers_3d(smiles: str, n_conformers: int = 5, style: str = "stick") -> None:
    """Render multiple conformers, colored by conformer index."""
    if not RDKIT_CONFORMERS_AVAILABLE:
        st.info("RDKit not available; cannot generate 3D conformers.")
        _render_2d_fallback(smiles)
        return

    try:
        mols = generate_conformers(smiles, n_conformers=int(n_conformers), method="ETKDG")
        pdbs = []
        for m in mols:
            optimize_conformer(m, force_field="MMFF")
            pdbs.append(conformer_to_pdb(m, conf_id=0))
    except Exception as e:  # noqa: BLE001
        st.error(f"Conformer generation failed: {e}")
        _render_2d_fallback(smiles)
        return

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
    html = _build_py3dmol_html_multi(pdbs, colors=colors, style=style, width=520, height=420)
    st.components.v1.html(html, height=440, scrolling=False)


def render_overlay_3d(smiles_list: List[str], colors: Optional[List[str]] = None, reference_idx: int = 0) -> None:
    """Render multiple molecules aligned (best-effort) and overlaid."""
    if not RDKIT_CONFORMERS_AVAILABLE:
        st.info("RDKit not available; cannot generate 3D overlays.")
        return

    if not smiles_list:
        st.error("No SMILES provided.")
        return

    mols = []
    pdbs: List[str] = []
    for smi in smiles_list:
        try:
            ms = generate_conformers(smi, n_conformers=1, method="ETKDG")
            m = ms[0]
            optimize_conformer(m, force_field="MMFF")
            mols.append(m)
        except Exception as e:  # noqa: BLE001
            st.error(f"Invalid/failed SMILES '{smi}': {e}")
            return

    try:
        align_molecules_to_reference(mols, reference_idx=int(reference_idx))
    except Exception:
        pass

    for m in mols:
        pdbs.append(conformer_to_pdb(m, conf_id=0))

    used_colors = colors or ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]
    html = _build_py3dmol_html_multi(pdbs, colors=used_colors, style="stick", width=520, height=420)
    st.components.v1.html(html, height=440, scrolling=False)


def render_protein_3d(
    pdb_string: str,
    style: str = "cartoon",
    width: int = 520,
    height: int = 420,
    model_format: str = "pdb",
) -> None:
    """Render a protein (or any structure file) from a PDB/mmCIF string."""
    if not pdb_string or not isinstance(pdb_string, str):
        st.error("No structure content provided.")
        return

    fmt = str(model_format or "pdb").lower()
    stl = (style or "cartoon").lower()

    try:
        if PY3DMOL_AVAILABLE:
            view = py3Dmol.view(width=int(width), height=int(height))
            view.addModel(pdb_string, fmt)
            if stl == "cartoon":
                view.setStyle({"cartoon": {"color": "spectrum"}})
            else:
                view.setStyle({stl: {}})
            view.setBackgroundColor("0xeeeeee")
            view.zoomTo()
            html = view._make_html()  # noqa: SLF001
        else:
            # Fall back to CDN-based 3Dmol; reuse builder for PDB.
            html = _build_py3dmol_html(pdb_string, style=style, width=width, height=height)
        st.components.v1.html(html, height=int(height) + 20, scrolling=False)
    except Exception as e:  # noqa: BLE001
        st.error(f"Protein 3D rendering failed: {e}")


__all__ = [
    "render_molecule_3d",
    "render_conformers_3d",
    "render_overlay_3d",
    "render_protein_3d",
    "_build_py3dmol_html",
]


