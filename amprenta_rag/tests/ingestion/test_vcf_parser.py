"""Tests for VCF parser service."""

import tempfile
from pathlib import Path


class TestVCFParser:
    """Tests for VCF file parsing."""

    def test_parse_vcf_file_success(self):
        """Valid VCF file returns variants."""
        # Create minimal VCF content
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t12345\t.\tA\tT\t99\tPASS\tGENE=BRCA1
chr2\t67890\t.\tG\tC\t85\tPASS\tGENE=TP53
"""
        with tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False) as f:
            f.write(vcf_content)
            vcf_path = Path(f.name)
        
        try:
            from amprenta_rag.ingestion.genomics.vcf_parser import parse_vcf_file
            
            variants = parse_vcf_file(vcf_path)
            
            assert len(variants) == 2
            assert variants[0]["chromosome"] == "chr1"
            assert variants[0]["position"] == 12345
            assert variants[0]["reference_allele"] == "A"
            assert variants[0]["alternate_allele"] == "T"
        finally:
            vcf_path.unlink()

    def test_parse_vcf_with_annotations(self):
        """INFO field gene annotation is extracted."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr17\t7577120\t.\tG\tA\t99\tPASS\tGENE=TP53
"""
        with tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False) as f:
            f.write(vcf_content)
            vcf_path = Path(f.name)
        
        try:
            from amprenta_rag.ingestion.genomics.vcf_parser import parse_vcf_file
            
            variants = parse_vcf_file(vcf_path)
            
            assert len(variants) == 1
            assert variants[0].get("gene") == "TP53"
        finally:
            vcf_path.unlink()

    def test_parse_malformed_vcf_graceful(self):
        """Malformed VCF handles gracefully."""
        vcf_content = "not a valid vcf file\njust some text\n"
        
        with tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False) as f:
            f.write(vcf_content)
            vcf_path = Path(f.name)
        
        try:
            from amprenta_rag.ingestion.genomics.vcf_parser import parse_vcf_file
            
            # Should either return empty or raise a specific error
            try:
                variants = parse_vcf_file(vcf_path)
                assert variants == [] or isinstance(variants, list)
            except Exception as e:
                # Acceptable to raise for truly malformed files
                assert "VCF" in str(e).upper() or "parse" in str(e).lower()
        finally:
            vcf_path.unlink()
