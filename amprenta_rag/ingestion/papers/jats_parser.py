"""
JATS XML parser for extracting full-text sections from PMC Open Access articles.

JATS (Journal Article Tag Suite) is the XML format used by PubMed Central
for full-text articles. This parser extracts sections like Abstract, Methods,
Results, Discussion, etc.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import List, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class PaperSection:
    """
    A section of a scientific paper.

    Attributes:
        title: Section title (e.g., "Introduction", "Methods")
        content: Section text content
        order: Display order (lower numbers first)
    """

    title: str
    content: str
    order: int = 0


@dataclass
class PaperContent:
    """
    Full content of a scientific paper.

    Attributes:
        sections: List of paper sections in order
        raw_text: Complete paper text (all sections concatenated)
    """

    sections: List[PaperSection] = field(default_factory=list)
    raw_text: str = ""


# Standard section types and their display order
SECTION_ORDER = {
    "abstract": 0,
    "introduction": 1,
    "intro": 1,
    "background": 2,
    "methods": 3,
    "materials": 3,
    "materials and methods": 3,
    "results": 4,
    "discussion": 5,
    "conclusions": 6,
    "conclusion": 6,
    "acknowledgments": 7,
    "acknowledgements": 7,
    "references": 8,
}


def parse_jats_xml(xml_string: str) -> Optional[PaperContent]:
    """
    Parse JATS XML to extract paper sections.

    Args:
        xml_string: JATS XML content as string

    Returns:
        PaperContent object with extracted sections, or None if parsing fails

    Example:
        >>> xml = '''<article>
        ...   <front><article-meta><abstract><p>Abstract text</p></abstract></article-meta></front>
        ...   <body><sec><title>Introduction</title><p>Intro text</p></sec></body>
        ... </article>'''
        >>> content = parse_jats_xml(xml)
        >>> print(content.sections[0].title)
        Abstract
    """
    try:
        root = ET.fromstring(xml_string)
        sections = []

        # Extract abstract (in <front> section)
        abstract_section = _extract_abstract(root)
        if abstract_section:
            sections.append(abstract_section)

        # Extract body sections
        body = root.find(".//body")
        if body is not None:
            body_sections = _extract_body_sections(body)
            sections.extend(body_sections)

        # Sort sections by order
        sections.sort(key=lambda s: s.order)

        # Create raw text (concatenate all sections)
        raw_text = "\n\n".join(
            f"{section.title}\n{section.content}" for section in sections
        )

        logger.info("[JATS] Parsed %d sections from JATS XML", len(sections))
        return PaperContent(sections=sections, raw_text=raw_text)

    except ET.ParseError as e:
        logger.error("[JATS] Failed to parse JATS XML: %r", e)
        return None
    except Exception as e:
        logger.error("[JATS] Unexpected error parsing JATS XML: %r", e)
        return None


def _extract_abstract(root: ET.Element) -> Optional[PaperSection]:
    """
    Extract abstract section from JATS XML.

    Args:
        root: Root XML element

    Returns:
        PaperSection for abstract, or None if not found
    """
    abstract = root.find(".//abstract")
    if abstract is None:
        return None

    # Extract all paragraph text
    paragraphs = []
    for elem in abstract.iter():
        if elem.tag == "p" and elem.text:
            paragraphs.append(elem.text.strip())

    if not paragraphs:
        return None

    content = "\n\n".join(paragraphs)
    return PaperSection(
        title="Abstract", content=content, order=SECTION_ORDER.get("abstract", 0)
    )


def _extract_body_sections(body: ET.Element) -> List[PaperSection]:
    """
    Extract all sections from paper body.

    Args:
        body: Body XML element

    Returns:
        List of PaperSection objects
    """
    sections = []

    for sec in body.findall(".//sec"):
        section = _parse_section(sec)
        if section:
            sections.append(section)

    return sections


def _parse_section(sec: ET.Element) -> Optional[PaperSection]:
    """
    Parse a single <sec> element into a PaperSection.

    Args:
        sec: Section XML element

    Returns:
        PaperSection object, or None if section has no content
    """
    # Get section title
    title_elem = sec.find("title")
    if title_elem is not None and title_elem.text:
        title = title_elem.text.strip()
    else:
        # Try to get from sec-type attribute
        sec_type = sec.get("sec-type", "")
        title = sec_type.replace("-", " ").title() if sec_type else "Section"

    # Extract all paragraph text
    paragraphs = []
    for elem in sec.iter():
        if elem.tag == "p" and elem.text:
            paragraphs.append(elem.text.strip())

    if not paragraphs:
        # Section has no content
        return None

    content = "\n\n".join(paragraphs)

    # Determine section order based on title
    title_lower = title.lower()
    order = SECTION_ORDER.get(title_lower, 99)  # Unknown sections go last

    return PaperSection(title=title, content=content, order=order)


def _get_text_recursive(element: ET.Element) -> str:
    """
    Recursively extract all text from an XML element.

    Args:
        element: XML element

    Returns:
        Concatenated text content
    """
    text_parts = []

    if element.text:
        text_parts.append(element.text.strip())

    for child in element:
        child_text = _get_text_recursive(child)
        if child_text:
            text_parts.append(child_text)
        if child.tail:
            text_parts.append(child.tail.strip())

    return " ".join(text_parts)

