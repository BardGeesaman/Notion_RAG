"""AI Extraction Tools API endpoints.

Provides REST API for:
- OCR text extraction from images/PDFs
- Web page content scraping
- Entity normalization (genes, compounds, diseases)
"""

from __future__ import annotations

import ipaddress
import socket
from typing import Optional

from fastapi import APIRouter, Depends, File, HTTPException, Query, UploadFile, status
from pydantic import BaseModel

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.database.models import User

router = APIRouter(prefix="/extraction", tags=["extraction"])


# --- Constants (P1 Security) ---
MAX_FILE_SIZE = 10_000_000  # 10MB
ALLOWED_CONTENT_TYPES = [
    "application/pdf",
    "image/png",
    "image/jpeg",
    "image/tiff",
    "image/gif",
]

# SSRF Prevention
BLOCKED_SCHEMES = {"file", "ftp", "gopher"}
PRIVATE_IP_RANGES = [
    # IPv4
    ipaddress.ip_network("10.0.0.0/8"),
    ipaddress.ip_network("172.16.0.0/12"),
    ipaddress.ip_network("192.168.0.0/16"),
    ipaddress.ip_network("127.0.0.0/8"),
    ipaddress.ip_network("169.254.0.0/16"),  # Link-local
    ipaddress.ip_network("0.0.0.0/8"),       # "This" network
    # IPv6
    ipaddress.ip_network("::1/128"),         # Loopback
    ipaddress.ip_network("fe80::/10"),       # Link-local
    ipaddress.ip_network("fc00::/7"),        # Unique local
]

# OCR Language validation (must match installed tessdata)
SUPPORTED_OCR_LANGUAGES = {
    "eng", "fra", "deu", "spa", "ita", "por", "nld", "rus",
    "chi_sim", "chi_tra", "jpn", "kor", "ara", "hin", "tha",
    "pol", "ces", "dan", "fin", "hun", "nor", "swe", "tur",
}


def _is_safe_url(url: str) -> tuple[bool, str]:
    """Validate URL is not internal/private (SSRF prevention)."""
    from urllib.parse import urlparse
    
    parsed = urlparse(url)
    
    # Block dangerous schemes
    if parsed.scheme.lower() in BLOCKED_SCHEMES:
        return False, f"Blocked scheme: {parsed.scheme}"
    
    # Must have hostname
    hostname = parsed.hostname
    if not hostname:
        return False, "Invalid URL: no hostname"
    
    # Block localhost variants
    if hostname.lower() in ("localhost", "0.0.0.0"):
        return False, "Blocked: localhost"
    
    # Check if hostname is already an IP address
    try:
        ip = ipaddress.ip_address(hostname)
        for network in PRIVATE_IP_RANGES:
            if ip in network:
                return False, "Blocked: private IP range"
    except ValueError:
        # Not an IP address, try to resolve hostname
        try:
            ip = ipaddress.ip_address(socket.gethostbyname(hostname))
            for network in PRIVATE_IP_RANGES:
                if ip in network:
                    return False, "Blocked: private IP range"
        except (socket.gaierror, ValueError):
            pass  # Let WebScraper handle DNS errors
    
    return True, ""


# --- Schemas ---

class OCRResponse(BaseModel):
    """OCR extraction result."""
    success: bool
    text: str
    word_count: int
    language: str
    error: Optional[str] = None


class ScrapeRequest(BaseModel):
    """Web scraping request."""
    url: str


class ScrapeResponse(BaseModel):
    """Web scraping result."""
    success: bool
    url: str
    title: Optional[str]
    author: Optional[str]
    content: str
    word_count: int
    error: Optional[str] = None


class NormalizeRequest(BaseModel):
    """Entity normalization request."""
    entity_type: str  # gene, compound, disease
    name: str


class GeneResult(BaseModel):
    """Normalized gene information."""
    symbol: str
    name: Optional[str]
    entrez_id: Optional[str]
    uniprot_id: Optional[str]
    organism: Optional[str]


class CompoundResult(BaseModel):
    """Normalized compound information."""
    name: str
    cid: Optional[int]
    inchi_key: Optional[str]
    smiles: Optional[str]
    molecular_formula: Optional[str]
    molecular_weight: Optional[float]


class NormalizeResponse(BaseModel):
    """Entity normalization result."""
    success: bool
    entity_type: str
    query: str
    result: Optional[dict] = None
    error: Optional[str] = None


# --- Endpoints ---

@router.post("/ocr", response_model=OCRResponse)
async def extract_ocr(
    file: UploadFile = File(...),
    language: str = Query("eng", description="Tesseract language code"),
    current_user: User = Depends(get_current_user),
) -> OCRResponse:
    """
    Extract text from uploaded image or PDF using OCR.
    
    P1 Security:
    - File size limit: 10MB
    - Allowed types: PDF, PNG, JPEG, TIFF, GIF
    """
    # Validate language code
    if language not in SUPPORTED_OCR_LANGUAGES:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Unsupported language: {language}. Supported: {sorted(SUPPORTED_OCR_LANGUAGES)}",
        )
    
    # P1: Validate content type
    if file.content_type not in ALLOWED_CONTENT_TYPES:
        raise HTTPException(
            status_code=status.HTTP_415_UNSUPPORTED_MEDIA_TYPE,
            detail=f"Unsupported file type: {file.content_type}. Allowed: {ALLOWED_CONTENT_TYPES}",
        )
    
    # Read file content
    content = await file.read()
    
    # P1: Validate file size
    if len(content) > MAX_FILE_SIZE:
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail=f"File too large. Maximum size: {MAX_FILE_SIZE // 1_000_000}MB",
        )
    
    try:
        from amprenta_rag.extraction.ocr_service import OCRService
        
        ocr = OCRService(language=language)
        
        # Handle PDF vs image
        if file.content_type == "application/pdf":
            text = ocr.extract_from_scanned_pdf(content)
        else:
            text = ocr.extract_from_image(content)
        
        return OCRResponse(
            success=True,
            text=text,
            word_count=len(text.split()),
            language=language,
        )
    except ImportError as e:
        return OCRResponse(
            success=False,
            text="",
            word_count=0,
            language=language,
            error=str(e),
        )
    except Exception as e:
        return OCRResponse(
            success=False,
            text="",
            word_count=0,
            language=language,
            error=f"OCR failed: {str(e)}",
        )


@router.post("/scrape", response_model=ScrapeResponse)
def scrape_url(
    request: ScrapeRequest,
    current_user: User = Depends(get_current_user),
) -> ScrapeResponse:
    """
    Scrape content from a URL.
    
    Rate limiting: Inherits WebScraper's 0.5s per-domain limit.
    """
    # SSRF Prevention
    is_safe, error = _is_safe_url(request.url)
    if not is_safe:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"URL blocked for security: {error}",
        )
    
    try:
        from amprenta_rag.extraction.web_scraper import WebScraper
        
        scraper = WebScraper()
        result = scraper.extract_from_url(request.url)
        
        return ScrapeResponse(
            success=result.success,
            url=result.url,
            title=result.title,
            author=result.author,
            content=result.content,
            word_count=result.word_count,
            error=result.error,
        )
    except Exception as e:
        return ScrapeResponse(
            success=False,
            url=request.url,
            title=None,
            author=None,
            content="",
            word_count=0,
            error=f"Scraping failed: {str(e)}",
        )


@router.post("/normalize", response_model=NormalizeResponse)
def normalize_entity(
    request: NormalizeRequest,
    current_user: User = Depends(get_current_user),
) -> NormalizeResponse:
    """
    Normalize entity name to standard identifiers.
    
    Supported entity types: gene, compound, disease
    """
    entity_type = request.entity_type.lower()
    
    if entity_type not in ["gene", "compound", "disease"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type: {entity_type}. Supported: gene, compound, disease",
        )
    
    try:
        from amprenta_rag.extraction.entity_normalizer import EntityNormalizer
        
        normalizer = EntityNormalizer()
        
        if entity_type == "gene":
            result = normalizer.normalize_gene(request.name)
            if result:
                return NormalizeResponse(
                    success=True,
                    entity_type=entity_type,
                    query=request.name,
                    result={
                        "symbol": result.symbol,
                        "name": result.name,
                        "entrez_id": result.entrez_id,
                        "uniprot_id": result.uniprot_id,
                        "organism": result.organism,
                    },
                )
        
        elif entity_type == "compound":
            result = normalizer.normalize_compound(request.name)
            if result:
                return NormalizeResponse(
                    success=True,
                    entity_type=entity_type,
                    query=request.name,
                    result={
                        "name": result.name,
                        "cid": result.cid,
                        "inchi_key": result.inchi_key,
                        "smiles": result.smiles,
                        "molecular_formula": result.molecular_formula,
                        "molecular_weight": result.molecular_weight,
                    },
                )
        
        elif entity_type == "disease":
            result = normalizer.normalize_disease(request.name)
            if result:
                return NormalizeResponse(
                    success=True,
                    entity_type=entity_type,
                    query=request.name,
                    result={
                        "name": result.name,
                        "mesh_id": result.mesh_id,
                        "doid": result.doid,
                    },
                )
        
        # Not found
        return NormalizeResponse(
            success=False,
            entity_type=entity_type,
            query=request.name,
            error=f"Could not normalize {entity_type}: {request.name}",
        )
    
    except Exception as e:
        return NormalizeResponse(
            success=False,
            entity_type=entity_type,
            query=request.name,
            error=f"Normalization failed: {str(e)}",
        )