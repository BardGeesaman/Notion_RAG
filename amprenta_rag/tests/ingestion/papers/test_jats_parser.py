"""
Unit tests for JATS XML parser.

Tests parsing of PMC JATS XML format with various section structures.
"""

from __future__ import annotations

import pytest

from amprenta_rag.ingestion.papers.jats_parser import PaperSection, parse_jats_xml


# Sample JATS XML for testing
SAMPLE_JATS_XML = """<?xml version="1.0" encoding="UTF-8"?>
<article xmlns:xlink="http://www.w3.org/1999/xlink">
    <front>
        <article-meta>
            <abstract>
                <p>This is the abstract of the paper. It provides a brief overview.</p>
                <p>Second abstract paragraph with more details.</p>
            </abstract>
        </article-meta>
    </front>
    <body>
        <sec sec-type="intro">
            <title>Introduction</title>
            <p>Introduction paragraph explaining the background.</p>
            <p>Second introduction paragraph.</p>
        </sec>
        <sec sec-type="methods">
            <title>Materials and Methods</title>
            <p>Description of experimental methods.</p>
        </sec>
        <sec sec-type="results">
            <title>Results</title>
            <p>Presentation of experimental results.</p>
            <p>Additional results paragraph.</p>
        </sec>
        <sec sec-type="discussion">
            <title>Discussion</title>
            <p>Interpretation of results.</p>
        </sec>
        <sec sec-type="conclusions">
            <title>Conclusions</title>
            <p>Summary and conclusions.</p>
        </sec>
    </body>
</article>
"""

MINIMAL_JATS_XML = """<?xml version="1.0" encoding="UTF-8"?>
<article>
    <front>
        <article-meta>
            <abstract>
                <p>Minimal abstract.</p>
            </abstract>
        </article-meta>
    </front>
    <body>
        <sec>
            <title>Introduction</title>
            <p>Minimal introduction.</p>
        </sec>
    </body>
</article>
"""

NO_ABSTRACT_JATS_XML = """<?xml version="1.0" encoding="UTF-8"?>
<article>
    <front>
        <article-meta>
        </article-meta>
    </front>
    <body>
        <sec>
            <title>Introduction</title>
            <p>Introduction without abstract.</p>
        </sec>
    </body>
</article>
"""


class TestJATSParser:
    """Tests for JATS XML parser."""

    def test_parse_complete_jats_xml(self):
        """Test parsing complete JATS XML with all sections."""
        content = parse_jats_xml(SAMPLE_JATS_XML)

        assert content is not None
        assert len(content.sections) == 6
        
        # Check section titles
        section_titles = [s.title for s in content.sections]
        assert "Abstract" in section_titles
        assert "Introduction" in section_titles
        assert "Materials and Methods" in section_titles
        assert "Results" in section_titles
        assert "Discussion" in section_titles
        assert "Conclusions" in section_titles

    def test_sections_ordered_correctly(self):
        """Test that sections are ordered correctly."""
        content = parse_jats_xml(SAMPLE_JATS_XML)

        assert content is not None
        assert content.sections[0].title == "Abstract"
        assert content.sections[1].title == "Introduction"
        # Methods comes after intro
        assert "method" in content.sections[2].title.lower()
        # Results come after methods
        assert content.sections[3].title == "Results"
        # Discussion comes after results
        assert content.sections[4].title == "Discussion"
        # Conclusions come last
        assert content.sections[5].title == "Conclusions"

    def test_abstract_extraction(self):
        """Test abstract is extracted correctly."""
        content = parse_jats_xml(SAMPLE_JATS_XML)

        assert content is not None
        abstract = content.sections[0]
        assert abstract.title == "Abstract"
        assert "brief overview" in abstract.content
        assert "Second abstract paragraph" in abstract.content

    def test_section_content_extraction(self):
        """Test section content includes all paragraphs."""
        content = parse_jats_xml(SAMPLE_JATS_XML)

        assert content is not None
        intro = next(s for s in content.sections if s.title == "Introduction")
        assert "Introduction paragraph" in intro.content
        assert "Second introduction paragraph" in intro.content

    def test_raw_text_concatenation(self):
        """Test raw_text includes all sections."""
        content = parse_jats_xml(SAMPLE_JATS_XML)

        assert content is not None
        assert len(content.raw_text) > 0
        assert "Abstract" in content.raw_text
        assert "Introduction" in content.raw_text
        assert "Results" in content.raw_text

    def test_minimal_jats_xml(self):
        """Test parsing minimal JATS XML."""
        content = parse_jats_xml(MINIMAL_JATS_XML)

        assert content is not None
        assert len(content.sections) == 2
        assert content.sections[0].title == "Abstract"
        assert content.sections[1].title == "Introduction"

    def test_no_abstract_jats_xml(self):
        """Test parsing JATS XML without abstract."""
        content = parse_jats_xml(NO_ABSTRACT_JATS_XML)

        assert content is not None
        assert len(content.sections) == 1
        assert content.sections[0].title == "Introduction"

    def test_invalid_xml_returns_none(self):
        """Test that invalid XML returns None."""
        invalid_xml = "This is not valid XML <<<<>"
        content = parse_jats_xml(invalid_xml)

        assert content is None

    def test_empty_sections_ignored(self):
        """Test that sections with no paragraphs are ignored."""
        xml_empty_section = """<?xml version="1.0" encoding="UTF-8"?>
        <article>
            <body>
                <sec>
                    <title>Empty Section</title>
                </sec>
                <sec>
                    <title>With Content</title>
                    <p>Some content here.</p>
                </sec>
            </body>
        </article>
        """
        content = parse_jats_xml(xml_empty_section)

        assert content is not None
        # Only section with content should be included
        assert len(content.sections) == 1
        assert content.sections[0].title == "With Content"

    def test_paper_section_dataclass(self):
        """Test PaperSection dataclass."""
        section = PaperSection(
            title="Test Section", content="Test content", order=5
        )

        assert section.title == "Test Section"
        assert section.content == "Test content"
        assert section.order == 5

