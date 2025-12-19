from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT, repo_rate_limit
from amprenta_rag.ingestion.repositories.base import RepositoryInterface
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.repository import DataFile, StudyMetadata

logger = get_logger(__name__)

API_BASE = "https://www.ebi.ac.uk/biostudies/api/v1/"
DEFAULT_TIMEOUT = 30


def _headers() -> Dict[str, str]:
    return {"User-Agent": REPOSITORY_USER_AGENT}


def _safe_get(url: str, params: Optional[Dict[str, Any]] = None) -> Optional[Dict[str, Any]]:
    repo_rate_limit()
    try:
        resp = requests.get(url, params=params or {}, headers=_headers(), timeout=DEFAULT_TIMEOUT)
        resp.raise_for_status()
        return resp.json()
    except Exception as exc:
        logger.warning("[REPO][ARRAYEXPRESS] Request failed for %s: %r", url, exc)
        return None


class ArrayExpressRepository(RepositoryInterface):
    """ArrayExpress/EBI BioStudies repository implementation."""

    def get_repository_name(self) -> str:
        return "ArrayExpress"

    def get_omics_type(self) -> str:
        return "transcriptomics"

    def search_studies(
        self,
        keywords: List[str],
        filters: Optional[Dict[str, Any]] = None,
        max_results: int = 100,
    ) -> List[str]:
        query = " ".join(keywords) if keywords else ""
        if filters:
            # Basic filter concat (e.g., disease/organism) as free text
            for k, v in filters.items():
                if v:
                    query += f" {k}:{v}"
        params = {"query": query.strip(), "pageSize": max_results}
        url = f"{API_BASE}studies"
        data = _safe_get(url, params=params)
        if not data:
            return []
        studies = data.get("studies") or data.get("entries") or []
        results: List[str] = []
        for s in studies:
            acc = s.get("accno") or s.get("id") or s.get("accession")
            if acc:
                results.append(acc)
            if len(results) >= max_results:
                break
        logger.info("[REPO][ARRAYEXPRESS] search '%s' -> %d results", query, len(results))
        return results

    def fetch_study_metadata(self, study_id: str) -> Optional[StudyMetadata]:
        url = f"{API_BASE}studies/{study_id}"
        data = _safe_get(url)
        if not data:
            return None

        title = data.get("title") or data.get("name") or study_id
        description = data.get("description") or ""
        organism = None
        attrs = data.get("attributes", {})
        if isinstance(attrs, dict):
            organism = attrs.get("Organism") or attrs.get("organism")

        files = self.fetch_study_data_files(study_id, data=data)

        return StudyMetadata(
            study_id=study_id,
            repository=self.get_repository_name(),
            title=title,
            description=description,
            organism=organism,
            data_files=files,
            raw_metadata=data,
        )

    def fetch_study_data_files(
        self,
        study_id: str,
        file_types: Optional[List[str]] = None,
        data: Optional[Dict[str, Any]] = None,
    ) -> List[DataFile]:
        # If metadata already fetched, reuse it
        info = data or _safe_get(f"{API_BASE}studies/{study_id}") or {}
        files = info.get("files") or []
        results: List[DataFile] = []

        for f in files:
            fname = f.get("name") or f.get("fileName") or ""
            if file_types and fname:
                if not any(fname.lower().endswith(ext.lower()) for ext in file_types):
                    continue
            download_url = f.get("url") or f.get("link")
            size = f.get("size")
            file_id = f.get("accno") or fname
            results.append(
                DataFile(
                    file_id=str(file_id),
                    file_name=fname,
                    file_type=f.get("type") or "",
                    download_url=download_url,
                    size=size,
                )
            )
        return results

    def download_data_file(self, study_id: str, file_info: DataFile, dest_dir: str) -> bool:
        if not file_info.download_url:
            return False
        repo_rate_limit()
        dest_path = Path(dest_dir) / (file_info.file_name or file_info.file_id)
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with requests.get(file_info.download_url, stream=True, headers=_headers(), timeout=300) as r:
                r.raise_for_status()
                with open(dest_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            logger.info(
                "[REPO][ARRAYEXPRESS] Downloaded %s -> %s",
                file_info.file_id,
                dest_path,
            )
            return True
        except Exception as exc:
            logger.warning(
                "[REPO][ARRAYEXPRESS] Failed to download %s: %r",
                file_info.file_id,
                exc,
            )
            return False


__all__ = ["ArrayExpressRepository"]
