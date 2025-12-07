# amprenta_rag/metadata/classify_literature.py

from __future__ import annotations

import json
import re
from typing import Any, Dict, List, Optional

import requests

# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.openai_client import (get_default_models,
                                                get_openai_client)
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

def notion_headers() -> Dict[str, str]:
    """DEPRECATED: Notion support removed. Returns empty headers dict."""
    logger.debug("[METADATA][CLASSIFY-LITERATURE] notion_headers() deprecated - Notion support removed")
    return {}


# ----------------- Notion helpers ----------------- #


def _notion_base_url() -> str:
    return get_config().notion.base_url


def _query_literature_pages(limit: Optional[int] = None) -> List[Dict[str, Any]]:
    """
    Fetch all pages from the Literature DB (or up to `limit` pages if provided).
    """
    cfg = get_config().notion
    lit_db_id = cfg.lit_db_id
    if not lit_db_id:
        raise RuntimeError("NOTION_LIT_DB_ID not configured in config.py")

    url = f"{_notion_base_url()}/databases/{lit_db_id}/query"
    payload: Dict[str, Any] = {"page_size": 100}
    pages: List[Dict[str, Any]] = []
    next_cursor: Optional[str] = None

    while True:
        body = payload.copy()
        if next_cursor:
            body["start_cursor"] = next_cursor

        resp = requests.post(url, headers=notion_headers(), data=json.dumps(body))
        resp.raise_for_status()
        data = resp.json()

        batch = data.get("results", [])
        pages.extend(batch)

        if limit is not None and len(pages) >= limit:
            return pages[:limit]

        if not data.get("has_more"):
            break
        next_cursor = data.get("next_cursor")

    return pages


def _get_title(props: Dict[str, Any]) -> str:
    t = props.get("Title", {}).get("title", []) or []
    return t[0].get("plain_text", "").strip() if t else ""


def _get_abstract(props: Dict[str, Any]) -> str:
    r = props.get("Abstract", {}).get("rich_text", []) or []
    return r[0].get("plain_text", "").strip() if r else ""


def _get_select_name(props: Dict[str, Any], name: str) -> Optional[str]:
    sel = props.get(name, {}).get("select")
    if sel and sel.get("name"):
        return sel["name"]
    return None


def _get_multi_names(props: Dict[str, Any], name: str) -> List[str]:
    ms = props.get(name, {}).get("multi_select", []) or []
    return [x["name"] for x in ms if x.get("name")]


def _needs_classification(props: Dict[str, Any]) -> bool:
    """
    Decide if a Literature page should be classified.

    Default rule: if Disease AND Targets AND Lipid Signatures are all empty,
    we classify. You can override this with force=True at the call site.
    """
    diseases = _get_multi_names(props, "Disease")
    targets = _get_multi_names(props, "Targets")
    lipid_sigs = _get_multi_names(props, "Lipid Signatures")
    return not (diseases or targets or lipid_sigs)


def _update_literature_properties(page_id: str, updates: Dict[str, Any]) -> None:
    """
    Patch Notion page with new select/multi-select/number values.
    `updates` should map property name -> Notion property object.
    """
    url = f"{_notion_base_url()}/pages/{page_id}"
    payload = {"properties": updates}
    resp = requests.patch(url, headers=notion_headers(), data=json.dumps(payload))
    if resp.status_code >= 300:
        print(f"⚠️ Failed to update Literature page {page_id}: {resp.text}")


# ----------------- OpenAI classifier ----------------- #


def _build_classification_prompt(title: str, abstract: str) -> str:
    """
    Build a text prompt summarizing the document for classification.
    For now we use Title + Abstract; if you want, we can later pull RAG chunks.
    """
    text = f"Title: {title}\n\nAbstract: {abstract}\n"
    return text[:8000]  # just in case


def _classify_literature_doc(title: str, abstract: str) -> Dict[str, Any]:
    """
    Call OpenAI to classify a literature document into semantic + lipid metadata.

    Returns a dict with keys:
      - disease (list of strings)
      - targets (list of strings)
      - modality (list of strings)
      - stage (string or null)
      - model_systems (list of strings)
      - biomarker_role (list of strings)
      - importance (int 1-5 or null)
      - lipids_raw (list of strings)
      - lipid_signatures (list of strings)
      - lipid_signature_role (list of strings)
      - phenotype_axes (list of strings)
      - matrix (list of strings)
      - treatment_arms (list of strings)
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    content_text = _build_classification_prompt(title, abstract)

    system_prompt = """You are helping classify ALS/neurodegeneration research papers into structured metadata.

You must respond with a *single JSON object* with these keys:

- "disease": list of strings from {"ALS","AD","PD","MS","FTD","Other"}.

- "targets": list of gene/protein/complex names that could be targeted by a drug or genetic intervention
  (e.g. ["SPTLC1","SPTLC2","SPTLC3","DEGS1","CERS2","ORMDL3","ASAH1","TNF","IL17RA"]).
  DO NOT include lipid/metabolite names like "ceramide", "sphingomyelin", or "hexosylceramide" here.
  Lipid names MUST go only in "lipids_raw".

- "modality": list from {"SmallMolecule","ASO","GeneTherapy","Biologic","CellTherapy","Other"}.

- "stage": one of {"Preclinical","Phase1","Phase2","Phase3","Translational","Mechanistic","Epidemiology","Unknown"}.

- "model_systems": list of models/matrices (e.g. ["HumanCSF","HumanPlasma","WobblerMouse","iPSCNeurons","Mouse","Rat"]).

- "biomarker_role": list from {"Diagnostic","Prognostic","PD","Stratification","None"}.

- "importance": integer 1-5 (5 = central to ALS ceramide/sphingolipid story, 1 = peripheral or tangential).

- "lipids_raw": list of lipid species names as they appear (e.g. ["Cer(d18:1/16:0)","SM(d18:1/24:1)","C16:0 ceramide"]).
  If the paper discusses ceramide, sphingomyelin, glycosphingolipids, etc., put those here, NOT in "targets".

- "lipid_signatures": list of short names for lipid panels/signatures, or [] if not clearly defined.

- "lipid_signature_role": list from {"Diagnostic","Prognostic","PD","Stratification","Mechanistic"}.

- "phenotype_axes": list of phenotype or readout axes (e.g. ["ALSFRS-R_slope","NfL_baseline","PC1_CeramidePanel"]).

- "matrix": list from {"CSF","Plasma","Serum","Brain","SpinalCord","MotorCortex","Skin","Fibroblasts","Other"}.

- "treatment_arms": list of treatment arms or conditions if relevant.

If the paper does not clearly specify a field, use an empty list [] or "Unknown" rather than guessing wildly.

Focus on ceramide/sphingolipid biology and neurodegeneration (especially ALS) when assigning importance.
"""

    user_prompt = f"""Classify the following document:

{content_text}

Return ONLY the JSON object, nothing else."""

    resp = client.chat.completions.create(
        model=chat_model,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        temperature=0.1,
    )

    raw = resp.choices[0].message.content.strip()  # type: ignore[union-attr]

    # Extract JSON (robust to extra text, though we asked for none)
    match = re.search(r"\{.*\}", raw, flags=re.DOTALL)
    if match:
        raw = match.group(0)

    try:
        data = json.loads(raw)
        if not isinstance(data, dict):
            raise ValueError("Classification response is not a JSON object")
        return data
    except Exception as e:
        print("⚠️ Failed to parse classification JSON:", e)
        print("Raw content was:", raw[:500])
        # Fallback: empty metadata
        return {
            "disease": [],
            "targets": [],
            "modality": [],
            "stage": "Unknown",
            "model_systems": [],
            "biomarker_role": [],
            "importance": None,
            "lipids_raw": [],
            "lipid_signatures": [],
            "lipid_signature_role": [],
            "phenotype_axes": [],
            "matrix": [],
            "treatment_arms": [],
        }


# ----------------- Target / lipid post-processing ----------------- #

LIPID_LIKE_TERMS = {
    "ceramide",
    "ceramides",
    "sphingomyelin",
    "sphingomyelins",
    "hexcer",
    "hexosylceramide",
    "hexosylceramides",
    "glycosphingolipid",
    "glycosphingolipids",
    "sphingolipid",
    "sphingolipids",
    "dhcer",
    "dihydroceramide",
    "sm",
    "glccer",
    "glucosylceramide",
    "galcer",
    "galactosylceramide",
    "laccer",
    "lactosylceramide",
}


def _filter_targets_for_molecular_entities(raw_targets: List[str]) -> List[str]:
    """
    Keep only molecular targets (genes/proteins/enzymes/receptors) and drop
    generic lipid/metabolite terms like 'ceramide' or 'sphingomyelin'.

    Heuristics:
      - Drop if lowercased term is in LIPID_LIKE_TERMS.
      - Keep short all-caps-with-digits tokens (e.g. 'SPTLC1','DEGS1','CERS2').
      - Keep other names that don't look like generic lipids (e.g. 'ORMDL3').
    """
    cleaned: List[str] = []
    for t in raw_targets:
        if not t:
            continue
        t_str = str(t).strip()
        if not t_str:
            continue

        low = t_str.lower()
        if low in LIPID_LIKE_TERMS:
            # treat as lipid, not target
            continue

        # Likely gene symbol: short, uppercase, with digits/hyphens allowed
        if len(t_str) <= 16 and re.match(r"^[A-Z0-9\-]+$", t_str):
            cleaned.append(t_str)
            continue

        # Otherwise, keep if it's not obviously a generic lipid term
        cleaned.append(t_str)

    return cleaned


# ----------------- Build Notion updates ----------------- #


def _build_notion_updates_from_classification(cls: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert classifier output into Notion property updates.
    """

    def ms(name: str) -> Dict[str, Any]:
        vals = cls.get(name, []) or []
        return {"multi_select": [{"name": str(v)} for v in vals if v]}

    def sel(name: str) -> Dict[str, Any]:
        val = cls.get(name)

        # Handle list or single value from the classifier
        if isinstance(val, list):
            val = val[0] if val else None

        if not val or val == "Unknown":
            return {"select": None}

        # Ensure it's a string
        val_str = str(val)
        return {"select": {"name": val_str}}

    # Build Targets with lipid-stripping logic
    raw_targets = cls.get("targets", []) or []
    cleaned_targets = _filter_targets_for_molecular_entities(raw_targets)
    targets_prop = {"multi_select": [{"name": t} for t in cleaned_targets]}

    updates: Dict[str, Any] = {
        "Disease": ms("disease"),
        "Targets": targets_prop,
        "Modality": ms("modality"),
        "Stage": sel("stage"),
        "Model Systems": ms("model_systems"),
        "Biomarker Role": ms("biomarker_role"),
        "Lipid Species (raw)": ms("lipids_raw"),
        "Lipid Signatures": ms("lipid_signatures"),
        "Lipid Signature Role": ms("lipid_signature_role"),
        "Phenotype Axes": ms("phenotype_axes"),
        "Matrix": ms("matrix"),
        "Treatment Arms": ms("treatment_arms"),
    }

    importance = cls.get("importance")
    if isinstance(importance, int):
        updates["Importance"] = {"number": importance}

    return updates


# ----------------- Public entrypoint ----------------- #


def classify_and_update_all_literature(
    limit: Optional[int] = None,
    force: bool = False,
) -> None:
    """
    For each Literature DB page, call OpenAI once to classify semantic + lipid
    metadata, and write back to Notion.

    By default, only pages that appear unclassified (no Disease, no Targets,
    no Lipid Signatures) are processed.

    If force=True, all pages are classified and updated regardless of existing metadata.

    limit: if provided, only classify up to this many pages (for testing).
    """
    pages = _query_literature_pages(limit=limit)
    total = len(pages)
    print(f"📚 Found {total} literature page(s).")

    classified = 0
    skipped = 0

    for idx, page in enumerate(pages, 1):
        props = page.get("properties", {}) or {}
        title = _get_title(props)
        abstract = _get_abstract(props)
        page_id = page["id"]

        if not title and not abstract:
            print(f"\n[{idx}/{total}] (no title/abstract) – skipping.")
            skipped += 1
            continue

        if (not force) and (not _needs_classification(props)):
            print(
                f"\n[{idx}/{total}] {title[:80]} – already has semantic metadata, skipping."
            )
            skipped += 1
            continue

        print(f"\n[{idx}/{total}] Classifying: {title[:80]}")

        cls = _classify_literature_doc(title, abstract)
        updates = _build_notion_updates_from_classification(cls)
        _update_literature_properties(page_id, updates)

        classified += 1

    print("\n=====================================")
    print(f"✅ Classification complete.")
    print(f"   Classified: {classified}")
    print(f"   Skipped:    {skipped}")
    print("=====================================\n")
