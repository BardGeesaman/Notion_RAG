# LEGACY / TRANSITIONAL MODULE
#
# This file contains the original Zotero ingestion implementation.
# The active ingestion pipeline now lives in:
#   - amprenta_rag.ingestion.zotero_ingest.ingest_zotero_item
#   - amprenta_rag.ingestion.zotero_api
#   - amprenta_rag.ingestion.text_extraction
#   - amprenta_rag.ingestion.notion_pages
#   - amprenta_rag.ingestion.metadata_semantic
#   - amprenta_rag.ingestion.pinecone_utils
#
# No other modules should import from zotero_item.py going forward.
# It is kept only as a reference until all logic has been fully reviewed and migrated.

from datetime import datetime, timezone
from typing import Dict, Any, List, Optional

import hashlib
import re
import textwrap

import requests  # still used for Notion calls
# no PdfReader / BytesIO here anymore

from amprenta_rag.config import get_config
from amprenta_rag.clients.openai_client import get_openai_client, get_default_models
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.ingestion.zotero_api import (
    ZoteroItem,
    fetch_zotero_item,
    fetch_zotero_attachments,
    fetch_zotero_notes,
    download_zotero_file,
)

from amprenta_rag.ingestion.text_extraction import (
    extract_text_from_attachment_bytes,
    html_to_text,
    is_boilerplate,
)

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------





# ---------------------------------------------------------------------------
# Zotero helpers
# ---------------------------------------------------------------------------
def sanitize_metadata(meta: Dict[str, Any]) -> Dict[str, Any]:
    """
    Ensure all metadata values are allowed by Pinecone:
      - string
      - number (int/float)
      - boolean
      - list of strings

    Drop keys with None or empty lists.
    Convert non-string, non-primitive list items to strings.
    """
    cleaned: Dict[str, Any] = {}

    for key, value in meta.items():
        if value is None:
            # Pinecone does not allow null metadata
            continue

        # Primitive types are fine
        if isinstance(value, (str, bool, int, float)):
            cleaned[key] = value
            continue

        # Lists need to be list of strings
        if isinstance(value, list):
            string_list = [str(v) for v in value if v is not None]
            if string_list:
                cleaned[key] = string_list
            continue

        # Anything else (e.g. dicts) – convert to string for now
        cleaned[key] = str(value)

    return cleaned


# ---------------------------------------------------------------------------
# Text extraction from attachments & notes
# ---------------------------------------------------------------------------


def extract_text_from_pdf_bytes(data: bytes) -> str:
    reader = PdfReader(BytesIO(data))
    text_parts: List[str] = []
    for page in reader.pages:
        t = page.extract_text() or ""
        if t:
            text_parts.append(t)
    return "\n".join(text_parts)


def extract_text_from_attachment_bytes(
    content_type: Optional[str],
    data: bytes,
    filename: str = "",
) -> str:
    """
    Convert different attachment types to text.

    - PDFs: use PdfReader
    - text/*: decode as UTF-8
    - everything else: skip
    """
    ct = (content_type or "").lower()

    if ct == "application/pdf":
        return extract_text_from_pdf_bytes(data)

    if ct.startswith("text/"):
        try:
            return data.decode("utf-8", errors="ignore")
        except Exception as e:
            print(f"⚠️ Failed to decode text attachment {filename}: {e}")
            return ""

    print(f"ℹ️ Skipping non-text attachment {filename!r} ({content_type})")
    return ""


def html_to_text(html: str) -> str:
    """
    Very simple HTML → text conversion for Zotero notes.
    """
    # Strip comments
    html = re.sub(r"<!--.*?-->", " ", html, flags=re.DOTALL)
    # <br> → newline
    html = re.sub(r"<br\s*/?>", "\n", html, flags=re.IGNORECASE)
    # Remove tags
    text = re.sub(r"<[^>]+>", " ", html)
    # Normalize whitespace
    text = re.sub(r"\s+", " ", text)
    return text.strip()


# ---------------------------------------------------------------------------
# Embeddings & chunking
# ---------------------------------------------------------------------------


def _embed_texts(texts: List[str], batch_size: int = 64) -> List[List[float]]:
    """
    Embed texts in batches using the configured embedding model.
    """
    if not texts:
        return []

    client = get_openai_client()
    _, embed_model = get_default_models()
    all_embeddings: List[List[float]] = []

    for i in range(0, len(texts), batch_size):
        batch = texts[i : i + batch_size]
        resp = client.embeddings.create(
            model=embed_model,
            input=batch,
        )
        all_embeddings.extend(d.embedding for d in resp.data)  # type: ignore[attr-defined]

    return all_embeddings


def _chunk_text(text: str, max_chars: int = 2000) -> List[str]:
    """
    Paragraph-aware, character-limited chunking.
    """
    text = text.strip()
    if not text:
        return []

    paras = [p.strip() for p in text.split("\n") if p.strip()]
    chunks: List[str] = []
    buf = ""

    for p in paras:
        if buf and len(buf) + len(p) + 1 > max_chars:
            chunks.append(buf.strip())
            buf = p
        else:
            if buf:
                buf += "\n" + p
            else:
                buf = p

    if buf:
        chunks.append(buf.strip())
    return chunks


_BOILERPLATE_PATTERNS = [
    "all rights reserved",
    "no responsibility or liability",
    "creativecommons.org/licenses",
    "this is an open access article",
    "distributed under the terms of the creative commons",
    "medrxiv preprint",
    "biorxiv preprint",
    "author/funder",
    "supplementary table",
    "supplementary figure",
    "view this article online",
]


def is_boilerplate(text: str) -> bool:
    t = text.lower()
    return any(pat in t for pat in _BOILERPLATE_PATTERNS)


# ---------------------------------------------------------------------------
# Notion helpers (Literature & RAG chunk pages)
# ---------------------------------------------------------------------------


def _require_notion_db_id(name: str, value: Optional[str]) -> str:
    if not value:
        raise RuntimeError(f"Notion {name} DB ID is not configured in config.py")
    return value


def _ensure_literature_page(item: ZoteroItem, parent_type: str) -> str:
    """
    Find or create a Literature DB page for this Zotero item.

    Matches your schema:
      - Title (title)
      - Source Type (select)  -> "Paper"
      - Embedding Status (select) -> "Not Embedded"
      - Zotero Item Key (rich_text)
      - Abstract, Journal, DOI (url), Year, URL, Tags (optional)
    """
    cfg = get_config().notion
    lit_db_id = _require_notion_db_id("LIT", cfg.lit_db_id)

    # 1) Try find by Zotero Item Key
    query_url = f"{cfg.base_url}/databases/{lit_db_id}/query"
    payload = {
        "page_size": 1,
        "filter": {
            "property": "Zotero Item Key",
            "rich_text": {"equals": item.key},
        },
    }

    resp = requests.post(query_url, headers=notion_headers(), json=payload)
    if resp.status_code >= 300:
        logger.warning(
            "⚠️ Notion Literature query error for item %s: %s %s",
            item.key,
            resp.status_code,
            resp.text,
        )
        resp.raise_for_status()

    results = resp.json().get("results", [])
    if results:
        return results[0]["id"]

    # 2) Create new
    logger.info("Creating new Literature page for Zotero item %s", item.key)

    year_val = None
    if item.date:
        m = re.search(r"(19|20)\d{2}", item.date)
        if m:
            year_val = int(m.group(0))

    props: Dict[str, Any] = {
        "Title": {"title": [{"text": {"content": item.title or "Untitled"}}]},
        "Source Type": {"select": {"name": "Paper"}},
        "Embedding Status": {"select": {"name": "Not Embedded"}},
        "Zotero Item Key": {"rich_text": [{"text": {"content": item.key}}]},
    }

    if item.abstract:
        props["Abstract"] = {
            "rich_text": [{"text": {"content": item.abstract}}],
        }
    if item.journal:
        props["Journal"] = {
            "rich_text": [{"text": {"content": item.journal}}],
        }
    if item.doi:
        doi_str = item.doi.strip()
        if doi_str and not doi_str.startswith("http"):
            doi_str = "https://doi.org/" + doi_str
        props["DOI"] = {"url": doi_str}
    if year_val is not None:
        props["Year"] = {"number": year_val}
    if item.url:
        props["URL"] = {"url": item.url}
    if item.tags:
        props["Tags"] = {"multi_select": [{"name": t} for t in item.tags]}

    create_payload = {
        "parent": {"database_id": lit_db_id},
        "properties": props,
    }

    resp = requests.post(
        f"{cfg.base_url}/pages",
        headers=notion_headers(),
        json=create_payload,
    )
    if resp.status_code >= 300:
        logger.error(
            "⚠️ Notion Literature create error for item %s: %s %s",
            item.key,
            resp.status_code,
            resp.text,
        )
        resp.raise_for_status()

    return resp.json()["id"]


def update_literature_page(parent_page_id: str, when_iso: str) -> None:
    """
    Set Embedding Status = 'Embedded' and update Last Embedded.
    """
    import json

    cfg = get_config().notion
    payload = {
        "properties": {
            "Embedding Status": {"select": {"name": "Embedded"}},
            "Last Embedded": {"date": {"start": when_iso}},
        },
    }

    resp = requests.patch(
        f"{cfg.base_url}/pages/{parent_page_id}",
        headers=notion_headers(),
        data=json.dumps(payload),
    )
    if resp.status_code >= 300:
        print("⚠️ Notion Literature update error:", resp.text)


def create_rag_chunk_page(
    chunk_id: str,
    chunk_text: str,
    parent_type: str,
    parent_page_id: str,
    order: int,
    when_iso: str,
) -> str:
    """
    Create a chunk page in the RAG DB, matching your original schema:

      - Chunk ID (title)
      - Chunk Text (rich_text)
      - Parent Type (select)
      - Parent Item (relation)
      - Order (number)
      - Embedding Vector ID (rich_text)
      - Last Embedded (date)
    """
    import json

    cfg = get_config().notion
    rag_db_id = _require_notion_db_id("RAG", cfg.rag_db_id)

    props: Dict[str, Any] = {
        "Chunk ID": {"title": [{"text": {"content": chunk_id}}]},
        "Chunk Text": {"rich_text": [{"text": {"content": chunk_text[:1900]}}]},
        "Parent Type": {"select": {"name": parent_type}},
        "Parent Item": {"relation": [{"id": parent_page_id}]},
        "Order": {"number": order},
        "Embedding Vector ID": {"rich_text": [{"text": {"content": chunk_id}}]},
        "Last Embedded": {"date": {"start": when_iso}},
    }

    payload = {
        "parent": {"database_id": rag_db_id},
        "properties": props,
    }

    resp = requests.post(
        f"{cfg.base_url}/pages",
        headers=notion_headers(),
        data=json.dumps(payload),
    )

    if resp.status_code >= 300:
        print("⚠️ Notion RAG chunk create error:", resp.status_code, resp.text)
        resp.raise_for_status()

    page_id = resp.json()["id"].replace("-", "")
    return page_id

# ---------------------------------------------------------------------------
# Semantic + lipid metadata helpers (Notion → Pinecone)
# ---------------------------------------------------------------------------

def _fetch_notion_page(page_id: str) -> Dict[str, Any]:
    cfg = get_config().notion
    url = f"{cfg.base_url}/pages/{page_id}"
    resp = requests.get(url, headers=notion_headers())
    resp.raise_for_status()
    return resp.json()


def _get_select_name(props: Dict[str, Any], name: str) -> Optional[str]:
    sel = props.get(name, {}).get("select")
    if sel and sel.get("name"):
        return sel["name"]
    return None


def _get_multi_names(props: Dict[str, Any], name: str) -> List[str]:
    ms = props.get(name, {}).get("multi_select", []) or []
    return [x["name"] for x in ms if x.get("name")]

def get_literature_semantic_metadata(parent_page_id: str, item: ZoteroItem) -> Dict[str, Any]:
    """
    Read semantic + lipid metadata from the Literature DB page and build a
    doc-level + lipid-level metadata dict that we will mirror into Pinecone.

    Expects (but does not require) these Notion properties on the Literature DB:
      - Disease (multi-select)
      - Targets (multi-select)
      - Modality (multi-select)
      - Stage (select)
      - Model Systems (multi-select)
      - Biomarker Role (multi-select)
      - Importance (number)
      - Lipid Species (raw) (multi-select or text)
      - Canonical Lipid Species (relation) – optional
      - Lipid Signatures (multi-select)
      - Lipid Signature Role (multi-select)
      - Phenotype Axes (multi-select)
      - Matrix (multi-select)
      - Treatment Arms (multi-select)
    """
    page = _fetch_notion_page(parent_page_id)
    props = page.get("properties", {}) or {}

    # --- Basic semantic fields ---
    diseases = _get_multi_names(props, "Disease")
    targets = _get_multi_names(props, "Targets")
    modality = _get_multi_names(props, "Modality")
    stage = _get_select_name(props, "Stage")
    model_systems = _get_multi_names(props, "Model Systems")
    biomarker_role = _get_multi_names(props, "Biomarker Role")

    importance = props.get("Importance", {}).get("number")

    phenotype_axes = _get_multi_names(props, "Phenotype Axes")
    matrix = _get_multi_names(props, "Matrix")
    treatment_arms = _get_multi_names(props, "Treatment Arms")

    # --- Lipid-related fields ---
    lipids_raw = _get_multi_names(props, "Lipid Species (raw)")

    # Canonical lipids via relation (store Notion page IDs for now)
    canonical_rel = props.get("Canonical Lipid Species", {}).get("relation", []) or []
    canonical_lipids = [r.get("id", "").replace("-", "") for r in canonical_rel if r.get("id")]

    lipid_signatures = _get_multi_names(props, "Lipid Signatures")
    lipid_signature_role = _get_multi_names(props, "Lipid Signature Role")

    # --- Simple doc_type inference from Zotero item_type and tags ---
    doc_type = "JournalArticle"
    if item.item_type.lower() == "book":
        doc_type = "Book"
    elif item.item_type.lower() == "thesis":
        doc_type = "Thesis"
    elif any("review" in t.lower() for t in item.tags):
        doc_type = "Review"
    elif any("trial" in t.lower() for t in item.tags):
        doc_type = "ClinicalTrialDesign"

    # Extract year from date string if available
    year_val = None
    if item.date:
        m = re.search(r"(19|20)\d{2}", item.date)
        if m:
            year_val = int(m.group(0))

    doc_meta: Dict[str, Any] = {
        "doc_id": f"ZOTERO:{item.key}",
        "doc_source": "Zotero",
        "doc_source_subtype": "Attachment",  # notes will override if needed
        "doc_type": doc_type,
        "diseases": diseases,
        "targets": targets,
        "modality": modality,
        "stage": stage,
        "model_systems": model_systems,
        "biomarker_role": biomarker_role,
        "year": year_val,
        "journal": item.journal,
        "importance": importance,
        "manual_tags": item.tags,
    }

    lipid_meta: Dict[str, Any] = {
        "lipids_raw": lipids_raw,
        "lipids": canonical_lipids,
        "lipid_classes": [],  # can be filled later by deref'ing Lipid Species DB
        "lipid_signatures": lipid_signatures,
        "lipid_signature_role": lipid_signature_role,
        "phenotype_axes": phenotype_axes,
        "matrix": matrix,
        "treatment_arms": treatment_arms,
    }

    return {**doc_meta, **lipid_meta}

# ---------------------------------------------------------------------------
# Pinecone idempotency helpers
# ---------------------------------------------------------------------------


def _query_pinecone_by_filter(filter_obj: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Helper: query Pinecone by metadata filter using a dummy vector of the right dimension.
    We rely on the index being 3072-dim (your current index).
    """
    index = get_pinecone_index()
    cfg = get_config()

    dummy = [0.0] * 3072  # matches your index dimension
    res = index.query(
        vector=dummy,
        top_k=10,
        include_metadata=True,
        namespace=cfg.pinecone.namespace,
        filter=filter_obj,
    )

    matches = getattr(res, "matches", None)
    if matches is None:
        matches = res.get("matches", [])
    return list(matches)


def attachment_already_ingested(
    item_key: str,
    attachment_key: str,
    current_md5: Optional[str] = None,
) -> bool:
    """
    Return True if this (item_key, attachment_key) is already ingested AND
    md5 matches. If md5 changed, delete old vectors and return False.
    """
    matches = _query_pinecone_by_filter(
        {
            "zotero_item_key": {"$eq": item_key},
            "attachment_key": {"$eq": attachment_key},
        }
    )
    if not matches:
        return False

    if not current_md5:
        return True

    # Look at first match metadata
    m0 = matches[0]
    meta = getattr(m0, "metadata", None) or m0.get("metadata", {})
    stored_md5 = meta.get("attachment_md5")

    index = get_pinecone_index()
    cfg = get_config()

    if (not stored_md5) or (stored_md5 != current_md5):
        print(
            f"⚠️ Attachment {attachment_key} for item {item_key} has changed "
            f"(stored md5={stored_md5!r}, current md5={current_md5!r}); "
            f"deleting old vectors and re-ingesting."
        )
        index.delete(
            filter={
                "zotero_item_key": {"$eq": item_key},
                "attachment_key": {"$eq": attachment_key},
            },
            namespace=cfg.pinecone.namespace,
        )
        return False

    return True


def note_already_ingested(
    item_key: str,
    note_key: str,
    note_hash: str,
) -> bool:
    """
    Same pattern as attachments: check if note with this hash is already in Pinecone.
    If hash changed, delete old vectors.
    """
    matches = _query_pinecone_by_filter(
        {
            "zotero_item_key": {"$eq": item_key},
            "note_key": {"$eq": note_key},
        }
    )
    if not matches:
        return False

    m0 = matches[0]
    meta = getattr(m0, "metadata", None) or m0.get("metadata", {})
    stored_hash = meta.get("note_hash")

    index = get_pinecone_index()
    cfg = get_config()

    if (not stored_hash) or (stored_hash != note_hash):
        print(
            f"⚠️ Note {note_key} for item {item_key} has changed "
            f"(stored={stored_hash!r}, current={note_hash!r}); "
            f"deleting old vectors and re-ingesting."
        )
        index.delete(
            filter={
                "zotero_item_key": {"$eq": item_key},
                "note_key": {"$eq": note_key},
            },
            namespace=cfg.pinecone.namespace,
        )
        return False

    return True


# ---------------------------------------------------------------------------
# Main ingestion: attachments + notes
# ---------------------------------------------------------------------------


def ingest_zotero_item(
    item_key: str,
    parent_type: str = "Literature",
    force: bool = False,
) -> None:
    """
    Ingest a single Zotero item into Notion + Pinecone, with
    attachment-level and note-level idempotency.

    For each attachment OR note on the Zotero item:
      - Skip if already ingested with same content (md5 / hash) unless force=True
      - Attachments: download file, extract text (PDF or text/*)
      - Notes: HTML → text
      - Chunk + embed
      - Create RAG chunk pages in Notion
      - Upsert vectors with rich metadata

    If all attachments/notes are up to date, this is a no-op.
    """
    logger.info("📥 Ingesting Zotero item %s", item_key)
    cfg = get_config()
    index = get_pinecone_index()

    # 1) Item metadata
    item = fetch_zotero_item(item_key)

    # 2) Ensure Literature page exists
    parent_page_id = _ensure_literature_page(item, parent_type)

    # 2b) Read semantic + lipid metadata from the Literature page
    base_meta = get_literature_semantic_metadata(parent_page_id, item)

    # 3) Fetch children
    attachments = fetch_zotero_attachments(item_key)
    notes = fetch_zotero_notes(item_key)

    if not attachments and not notes:
        logger.info("   ℹ️ No attachments or notes for %s; nothing to ingest.", item_key)
        return

    print(
        f"📄 Found {len(attachments)} attachment(s) and {len(notes)} note(s) "
        f"for Zotero item {item_key}."
    )

    now = datetime.now(timezone.utc).isoformat()
    any_ingested = False

    # ---------------- Attachments ---------------- #
    for att in attachments:
        att_key = att.get("key")
        if not att_key:
            continue

        filename = att.get("filename") or ""
        content_type = att.get("contentType")
        att_md5 = att.get("md5")

        print(f"\n📎 Attachment {att_key} ({filename or 'unnamed'})")

        # 🔥 Skip Zotero HTML Snapshots (web page dumps)
        title = att.get("title", "") or ""
        if "snapshot" in filename.lower() or "snapshot" in title.lower():
            print("   ⏭️ Skipping Snapshot attachment.")
            continue

        if content_type and content_type.lower().startswith("text/html"):
            print("   ⏭️ Skipping HTML snapshot attachment.")
            continue

        # Idempotency (skip unless force=True)
        if not force and attachment_already_ingested(
            item_key=item_key,
            attachment_key=att_key,
            current_md5=att_md5,
        ):
            print("   ⏭️ Already ingested with matching md5; skipping.")
            continue

        # Download file and extract text
        try:
            data = download_zotero_file(att_key)
        except Exception as e:
            print(f"   ❌ Failed to download attachment {att_key}: {e}")
            continue

        text = extract_text_from_attachment_bytes(content_type, data, filename=filename)
        if not text or len(text.strip()) < 100:
            print("   ⏭️ Attachment text is empty or very short; skipping.")
            continue

        header = item.title + "\n"
        if item.journal:
            header += f"{item.journal}\n"
        if item.doi:
            header += f"DOI: {item.doi}\n"
        header += f"[Attachment: {filename}]\n\n"
        full_text = header + text

        chunks = _chunk_text(full_text)
        chunks = [c for c in chunks if not is_boilerplate(c)]
        if not chunks:
            print("   ⏭️ No usable chunks produced; skipping.")
            continue

        print(f"   ✂️ Generated {len(chunks)} chunk(s) for attachment {att_key}.")
        print("   🧠 Embedding chunks with OpenAI...")
        embeddings = _embed_texts(chunks)

        vectors: List[Dict[str, Any]] = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{parent_page_id}_{att_key}_chunk_{order:03d}"

            chunk_page_id = create_rag_chunk_page(
                chunk_id=chunk_id,
                chunk_text=chunk,
                parent_type=parent_type,
                parent_page_id=parent_page_id,
                order=order,
                when_iso=now,
            )

            snippet = textwrap.shorten(chunk, width=300)

            meta: Dict[str, Any] = {
                **base_meta,  # doc-level + lipid-level metadata from Notion
                # chunk-level:
                "chunk_id": chunk_id,
                "chunk_index": order,
                "snippet": snippet,
                # link-back fields:
                "zotero_item_key": item.key,
                "attachment_key": att_key,
                "attachment_filename": filename,
                "attachment_content_type": content_type,
                "attachment_md5": att_md5,
                "notion_chunk_page_id": chunk_page_id,
                "source": "Literature",
                "source_type": parent_type,
                "title": item.title,
                "zotero_tags": item.tags,
                "item_type": item.item_type,
            }
            if item.doi:
                meta["doi"] = item.doi
            if item.url:
                meta["url"] = item.url
            if item.date:
                meta["date"] = item.date

            vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(meta),
                }
            )

        if vectors:
            print(f"   📡 Upserting {len(vectors)} vectors into Pinecone...")
            index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
            any_ingested = True

    # ---------------- Notes ---------------- #
    for note in notes:
        note_key = note.get("key")
        if not note_key:
            continue

        raw_html = note.get("note") or ""
        text = html_to_text(raw_html)
        if not text or len(text.strip()) < 50:
            continue

        # Hash of note text for idempotency
        note_hash = hashlib.md5(text.encode("utf-8", errors="ignore")).hexdigest()
        print(f"\n📝 Note {note_key} (hash={note_hash})")

        if not force and note_already_ingested(
            item_key=item_key,
            note_key=note_key,
            note_hash=note_hash,
        ):
            print("   ⏭️ Note already ingested with same hash; skipping.")
            continue

        header = f"{item.title}\n[Zotero Note]\n\n"
        full_text = header + text

        chunks = _chunk_text(full_text)
        chunks = [c for c in chunks if not is_boilerplate(c)]
        if not chunks:
            print("   ⏭️ No usable chunks for note; skipping.")
            continue

        print(f"   ✂️ Generated {len(chunks)} chunk(s) for note {note_key}.")
        print("   🧠 Embedding note chunks with OpenAI...")
        embeddings = _embed_texts(chunks)

        vectors: List[Dict[str, Any]] = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{parent_page_id}_{note_key}_note_{order:03d}"

            chunk_page_id = create_rag_chunk_page(
                chunk_id=chunk_id,
                chunk_text=chunk,
                parent_type=parent_type,
                parent_page_id=parent_page_id,
                order=order,
                when_iso=now,
            )

            snippet = textwrap.shorten(chunk, width=300)

            meta: Dict[str, Any] = {
                **base_meta,
                "chunk_id": chunk_id,
                "chunk_index": order,
                "snippet": snippet,
                "zotero_item_key": item.key,
                "note_key": note_key,
                "note_hash": note_hash,
                "notion_chunk_page_id": chunk_page_id,
                "source": "Literature",
                "source_type": parent_type,
                "title": item.title,
                "zotero_tags": item.tags,
                "item_type": item.item_type,
            }
            if item.doi:
                meta["doi"] = item.doi
            if item.url:
                meta["url"] = item.url
            if item.date:
                meta["date"] = item.date

            vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(meta),
                }
            )

        if vectors:
            print(f"   📡 Upserting {len(vectors)} note vectors into Pinecone...")
            index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
            any_ingested = True

    # ---------------- Final status ---------------- #
    if any_ingested:
        print("📝 Updating Literature page status in Notion...")
        update_literature_page(parent_page_id, now)
        print("✅ Ingestion complete.")
    else:
        print("ℹ️ No new or changed attachments/notes; nothing ingested.")