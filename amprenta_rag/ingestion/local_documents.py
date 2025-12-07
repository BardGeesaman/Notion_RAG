import datetime
import textwrap
from pathlib import Path
from typing import Dict
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Email, Literature
from amprenta_rag.ingestion.pinecone_utils import get_pinecone_index, sanitize_metadata
from amprenta_rag.ingestion.postgres_rag_chunk import create_rag_chunk_in_postgres
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.ingestion.text_extraction import extract_text_from_pdf_bytes


def ingest_local_document(
    file_path: Path,
    doc_type: str,
    title: str,
    metadata: Dict,
) -> UUID:
    ext = file_path.suffix.lower()
    # Determine mime type (better error handling for .pdf, .txt, .md)
    if ext == ".pdf":
        text = extract_text_from_pdf_bytes(file_path.read_bytes())
    elif ext in {".txt", ".md"}:
        text = file_path.read_text(encoding="utf-8", errors="replace")
    else:
        raise ValueError(f"Unsupported file type: {ext}")
    if not text or len(text.strip()) < 50:
        raise ValueError("Document text is empty or too short after extraction.")
    # Compose semantic_metadata
    semantic_metadata = dict(metadata)
    # Save to Postgres
    db = next(get_db())
    rec = None
    if doc_type.lower() == "literature":
        rec = Literature(
            title=title,
            source_type="LocalUpload",
            semantic_metadata=semantic_metadata,
            created_at=datetime.datetime.utcnow(),
            updated_at=datetime.datetime.utcnow(),
            embedding_status="Ingesting",
        )
        db.add(rec)
        db.commit()
        db.refresh(rec)
    else:
        rec = Email(
            title=title,
            content=text,
            item_type=doc_type,
            semantic_metadata=semantic_metadata,
            created_at=datetime.datetime.utcnow(),
            updated_at=datetime.datetime.utcnow(),
            embedding_status="Ingesting",
        )
        db.add(rec)
        db.commit()
        db.refresh(rec)
    # Chunk, embed, and create RAG chunks:
    chunks = chunk_text(text)
    chunks = [c for c in chunks if c.strip()]
    if not chunks:
        raise ValueError("No substantial chunks generated from extracted text.")
    embeddings = embed_texts(chunks)
    vectors = []
    for idx, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_uuid = create_rag_chunk_in_postgres(
            chunk_id=f"{rec.id}_local_chunk_{idx:03d}",
            chunk_text=chunk,
            source_type="Literature" if doc_type.lower() == "literature" else doc_type.capitalize(),
            source_id=rec.id,
            source_name=title,
            chunk_index=idx,
            snippet=textwrap.shorten(chunk, width=300),
            chunk_metadata=semantic_metadata,
            literature_id=rec.id if doc_type.lower() == "literature" else None,
            email_id=rec.id if doc_type.lower() != "literature" else None,
        )
        vectors.append(
            {
                "id": f"{rec.id}_local_chunk_{idx:03d}",
                "values": emb,
                "metadata": sanitize_metadata(
                    {
                        **semantic_metadata,
                        "postgres_chunk_id": str(chunk_uuid),
                        "source_id": str(rec.id),
                        "doc_type": doc_type,
                        "title": title,
                    }
                ),
            }
        )
    pinecone_index = get_pinecone_index()
    pinecone_index.upsert(vectors=vectors)
    # Mark embedding/ingestion as complete
    rec.embedding_status = "Embedded"
    rec.last_ingested_at = datetime.datetime.utcnow()
    db.add(rec)
    db.commit()
    return rec.id
