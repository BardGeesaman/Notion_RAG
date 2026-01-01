"""OCR service for extracting text from scanned documents and images."""

from __future__ import annotations

import io
import logging
import shutil
from typing import List, Optional

from PIL import Image

logger = logging.getLogger(__name__)

# Check if tesseract is available (P1 FIX: graceful degradation)
TESSERACT_AVAILABLE = shutil.which('tesseract') is not None

if TESSERACT_AVAILABLE:
    import pytesseract
else:
    pytesseract = None  # type: ignore


class OCRService:
    """Extract text from scanned documents and images using Tesseract OCR."""

    def __init__(self, language: str = "eng"):
        """
        Initialize OCR service.
        
        Args:
            language: Tesseract language code (default: English)
            
        Raises:
            ImportError: If tesseract-ocr is not installed
        """
        if not TESSERACT_AVAILABLE:
            raise ImportError(
                "tesseract-ocr not installed. OCR functionality unavailable.\n"
                "Install with:\n"
                "  - Linux: apt-get install tesseract-ocr\n"
                "  - macOS: brew install tesseract\n"
                "  - Windows: Download from https://github.com/UB-Mannheim/tesseract/wiki"
            )
        self.language = language

    def extract_from_image(self, image_bytes: bytes) -> str:
        """
        Extract text from an image using OCR.
        
        Args:
            image_bytes: Image file content as bytes
            
        Returns:
            Extracted text string
        """
        try:
            image = Image.open(io.BytesIO(image_bytes))
            text = pytesseract.image_to_string(image, lang=self.language)
            logger.info(f"OCR extracted {len(text)} chars from image")
            return text.strip()
        except Exception as e:
            logger.error(f"OCR extraction failed: {e}")
            return ""

    def extract_from_scanned_pdf(self, pdf_bytes: bytes) -> str:
        """
        Extract text from a scanned PDF using OCR.
        
        Converts PDF pages to images then runs OCR on each.
        
        Args:
            pdf_bytes: PDF file content as bytes
            
        Returns:
            Extracted text from all pages concatenated
        """
        try:
            from pdf2image import convert_from_bytes
            
            # Convert PDF pages to images
            images = convert_from_bytes(pdf_bytes, dpi=300)
            logger.info(f"Converting {len(images)} PDF pages to images for OCR")
            
            # OCR each page
            texts: List[str] = []
            for i, image in enumerate(images):
                page_text = pytesseract.image_to_string(image, lang=self.language)
                texts.append(page_text)
                logger.debug(f"Page {i+1}: extracted {len(page_text)} chars")
            
            full_text = "\n\n".join(texts)
            logger.info(f"OCR extracted {len(full_text)} total chars from PDF")
            return full_text.strip()
            
        except ImportError:
            logger.error("pdf2image not installed. Install with: pip install pdf2image")
            logger.error("Also requires poppler-utils: apt-get install poppler-utils (Linux) or brew install poppler (macOS)")
            return ""
        except Exception as e:
            logger.error(f"PDF OCR extraction failed: {e}")
            return ""

    def is_scanned_pdf(self, pdf_bytes: bytes, min_text_ratio: float = 0.01) -> bool:
        """
        Detect if a PDF is scanned (contains mostly images, little embedded text).
        
        Args:
            pdf_bytes: PDF file content as bytes
            min_text_ratio: Minimum text-to-page ratio to consider non-scanned
            
        Returns:
            True if PDF appears to be scanned
        """
        try:
            from pypdf import PdfReader
            
            reader = PdfReader(io.BytesIO(pdf_bytes))
            total_pages = len(reader.pages)
            total_text = 0
            
            for page in reader.pages:
                text = page.extract_text() or ""
                total_text += len(text)
            
            # Heuristic: less than ~100 chars per page suggests scanned
            avg_chars_per_page = total_text / total_pages if total_pages > 0 else 0
            is_scanned = avg_chars_per_page < 100
            
            logger.debug(
                f"PDF scan detection: {total_pages} pages, {total_text} chars, "
                f"avg {avg_chars_per_page:.0f} chars/page, is_scanned={is_scanned}"
            )
            return is_scanned
            
        except Exception as e:
            logger.error(f"PDF scan detection failed: {e}")
            return False  # Assume not scanned on error


def is_ocr_available() -> bool:
    """Check if OCR functionality is available."""
    return TESSERACT_AVAILABLE


def get_ocr_service(language: str = "eng") -> Optional[OCRService]:
    """
    Get OCR service if available.
    
    Returns:
        OCRService instance or None if tesseract not installed
    """
    if not TESSERACT_AVAILABLE:
        logger.warning("OCR not available: tesseract-ocr not installed")
        return None
    return OCRService(language=language)
