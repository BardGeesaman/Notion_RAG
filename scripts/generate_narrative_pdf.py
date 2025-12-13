import sys
from pathlib import Path
import markdown
from fpdf import FPDF
import re

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

class PDF(FPDF):
    def header(self):
        self.set_font('helvetica', 'B', 12)
        self.cell(0, 10, 'Amprenta Narrative Walkthrough', border=False, align='C', new_x="LMARGIN", new_y="NEXT")
        self.ln(5)

    def footer(self):
        self.set_y(-15)
        self.set_font('helvetica', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}', align='C')

def clean_text_for_pdf(text: str) -> str:
    # Strip emojis and other non-latin-1 characters that standard PDF fonts handle poorly
    # This is a robust way to ensure PDF generation works without custom fonts
    return text.encode('latin-1', 'ignore').decode('latin-1')

def markdown_to_pdf_fpdf2(md_content: str, output_path: Path):
    # Pre-process markdown to remove unsupported characters
    cleaned_md = clean_text_for_pdf(md_content)
    
    # Convert Markdown to HTML
    html = markdown.markdown(cleaned_md)
    
    pdf = PDF()
    pdf.add_page()
    pdf.set_font("helvetica", size=12)
    
    # Write HTML
    pdf.write_html(html)
    
    pdf.output(output_path)

def main():
    source_file = Path("docs/NARRATIVE_WALKTHROUGH.md")
    output_file = Path("docs/NARRATIVE_WALKTHROUGH.pdf")

    if not source_file.exists():
        print(f"Error: Source file {source_file} not found.")
        return

    print(f"Reading {source_file}...")
    content = source_file.read_text(encoding="utf-8")

    print(f"Converting to PDF using fpdf2...")
    try:
        markdown_to_pdf_fpdf2(content, output_file)
        print(f"Success! PDF saved to {output_file}")
    except Exception as e:
        print(f"Error converting to PDF: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
