"""Connectivity scoring service (Signature -> LINCS signatures)."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Dict, Optional
from uuid import UUID

from amprenta_rag.connectivity.scorer import compute_connectivity_score
from amprenta_rag.database.models import ConnectivityScore, LINCSSignature, Signature, SignatureComponent
from amprenta_rag.database.session import db_session


def _direction_sign(direction: Optional[str]) -> float:
    if not direction:
        return 1.0
    d = direction.strip().lower()
    if "down" in d or "decrease" in d or "↓" in d or d in ("-", "neg", "negative"):
        return -1.0
    if "up" in d or "increase" in d or "↑" in d or d in ("+", "pos", "positive"):
        return 1.0
    return 1.0


def _extract_entrez_from_external_ids(external_ids: object) -> Optional[int]:
    if not isinstance(external_ids, dict):
        return None
    for k in ("entrez_id", "entrez", "ncbi_gene_id", "gene_id", "ncbigene"):
        v = external_ids.get(k)
        if v is None:
            continue
        try:
            return int(v)
        except Exception:
            continue
    return None


class ConnectivityService:
    """Compute connectivity scores for an internal Signature against LINCS signatures."""

    def compute_scores(self, signature_id: UUID, *, score_type: str = "weighted_cosine") -> int:
        """Compute and store connectivity scores. Returns count stored."""
        with db_session() as db:
            sig = db.query(Signature).filter_by(id=signature_id).first()
            if not sig:
                raise ValueError("Signature not found")

            # Load components + features
            comps = db.query(SignatureComponent).filter_by(signature_id=signature_id).all()

            query_vec: Dict[int, float] = {}
            for c in comps:
                feat = getattr(c, "feature", None)
                if not feat:
                    continue
                if getattr(feat, "feature_type", None) != "gene":
                    continue
                entrez = _extract_entrez_from_external_ids(getattr(feat, "external_ids", None))
                if not entrez:
                    continue
                sign = _direction_sign(getattr(c, "direction", None))
                w = float(getattr(c, "weight", 1.0) or 1.0)
                query_vec[entrez] = query_vec.get(entrez, 0.0) + (sign * w)

            if not query_vec:
                raise ValueError("No gene Entrez IDs available for signature (missing Feature.external_ids.entrez_id)")

            # Clear existing scores for this signature+type
            db.query(ConnectivityScore).filter_by(query_signature_id=signature_id, score_type=score_type).delete()

            now = datetime.now(timezone.utc)
            count = 0
            batch = []

            for lincs_sig in db.query(LINCSSignature).yield_per(200):
                raw = getattr(lincs_sig, "gene_expression", None) or {}
                lincs_vec: Dict[int, float] = {}
                if isinstance(raw, dict):
                    for k, v in raw.items():
                        try:
                            lincs_vec[int(k)] = float(v)
                        except Exception:
                            continue

                score = compute_connectivity_score(query_vec, lincs_vec)
                batch.append(
                    ConnectivityScore(
                        query_signature_id=signature_id,
                        lincs_signature_id=lincs_sig.id,
                        score=float(score),
                        score_type=score_type,
                        p_value=None,
                        computed_at=now,
                    )
                )
                count += 1

                if len(batch) >= 1000:
                    db.bulk_save_objects(batch)
                    batch = []

            if batch:
                db.bulk_save_objects(batch)

            return count


__all__ = ["ConnectivityService"]


