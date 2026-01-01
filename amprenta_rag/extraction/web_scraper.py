"""Web scraper service for extracting content from URLs."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Optional
from urllib.parse import urlparse

import requests
from bs4 import BeautifulSoup
from readability import Document

logger = logging.getLogger(__name__)

# Default request timeout
REQUEST_TIMEOUT = 30

# Rate limiting: minimum seconds between requests to same domain
MIN_REQUEST_INTERVAL = 0.5


@dataclass
class ScrapedContent:
    """Result of web page scraping."""
    url: str
    title: Optional[str]
    content: str  # Main text content
    html: Optional[str]  # Raw HTML (optional)
    author: Optional[str]
    published_date: Optional[str]
    word_count: int
    success: bool
    error: Optional[str] = None


class WebScraper:
    """Extract content from web pages."""

    def __init__(self, user_agent: str = "Amprenta-RAG/1.0"):
        """
        Initialize web scraper.
        
        Args:
            user_agent: User agent string for requests
        """
        self.user_agent = user_agent
        self._last_request_time: dict[str, float] = {}  # domain -> timestamp

    def _rate_limit(self, domain: str) -> None:
        """Apply rate limiting per domain."""
        if domain in self._last_request_time:
            elapsed = time.time() - self._last_request_time[domain]
            if elapsed < MIN_REQUEST_INTERVAL:
                time.sleep(MIN_REQUEST_INTERVAL - elapsed)
        self._last_request_time[domain] = time.time()

    def extract_from_url(self, url: str, include_html: bool = False) -> ScrapedContent:
        """
        Fetch and extract main content from a URL.
        
        Args:
            url: Web page URL to scrape
            include_html: Whether to include raw HTML in result
            
        Returns:
            ScrapedContent with extracted text and metadata
        """
        try:
            # Parse domain for rate limiting
            parsed = urlparse(url)
            domain = parsed.netloc
            self._rate_limit(domain)
            
            # Fetch page
            logger.info(f"Fetching URL: {url}")
            response = requests.get(
                url,
                headers={"User-Agent": self.user_agent},
                timeout=REQUEST_TIMEOUT
            )
            response.raise_for_status()
            
            return self.extract_from_html(
                response.text, 
                url=url, 
                include_html=include_html
            )
            
        except requests.RequestException as e:
            logger.error(f"Failed to fetch {url}: {e}")
            return ScrapedContent(
                url=url,
                title=None,
                content="",
                html=None,
                author=None,
                published_date=None,
                word_count=0,
                success=False,
                error=str(e)
            )

    def extract_from_html(
        self, 
        html: str, 
        url: str = "", 
        include_html: bool = False
    ) -> ScrapedContent:
        """
        Extract main content from raw HTML.
        
        Uses readability-lxml for main content extraction.
        
        Args:
            html: Raw HTML string
            url: Original URL (for metadata)
            include_html: Whether to include raw HTML in result
            
        Returns:
            ScrapedContent with extracted text and metadata
        """
        try:
            # Use readability to extract main content
            doc = Document(html)
            title = doc.title()
            content_html = doc.summary()
            
            # Parse with BeautifulSoup for clean text
            soup = BeautifulSoup(content_html, "html.parser")
            content_text = soup.get_text(separator="\n", strip=True)
            
            # Extract metadata from original HTML
            original_soup = BeautifulSoup(html, "html.parser")
            author = self._extract_author(original_soup)
            published_date = self._extract_date(original_soup)
            
            word_count = len(content_text.split())
            
            logger.info(f"Extracted {word_count} words from {url or 'HTML'}")
            
            return ScrapedContent(
                url=url,
                title=title,
                content=content_text,
                html=content_html if include_html else None,
                author=author,
                published_date=published_date,
                word_count=word_count,
                success=True
            )
            
        except Exception as e:
            logger.error(f"Failed to extract content: {e}")
            return ScrapedContent(
                url=url,
                title=None,
                content="",
                html=None,
                author=None,
                published_date=None,
                word_count=0,
                success=False,
                error=str(e)
            )

    def _extract_author(self, soup: BeautifulSoup) -> Optional[str]:
        """Extract author from HTML metadata."""
        # Try meta tags
        for meta in soup.find_all("meta"):
            if meta.get("name", "").lower() in ["author", "article:author"]:
                return meta.get("content")
            if meta.get("property", "").lower() in ["author", "article:author"]:
                return meta.get("content")
        return None

    def _extract_date(self, soup: BeautifulSoup) -> Optional[str]:
        """Extract publication date from HTML metadata."""
        # Try meta tags
        for meta in soup.find_all("meta"):
            prop = meta.get("property", "").lower()
            name = meta.get("name", "").lower()
            if prop in ["article:published_time", "datePublished"] or \
               name in ["date", "pubdate", "publishdate"]:
                return meta.get("content")
        
        # Try time element
        time_elem = soup.find("time")
        if time_elem and time_elem.get("datetime"):
            return time_elem.get("datetime")
        
        return None
