"""Tests for BAM/CRAM parser."""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock


class TestBAMParser:
    """Tests for BAM file parsing functions."""

    def test_check_index_exists_finds_bai(self):
        """check_index_exists finds .bai index files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bam_path = Path(tmpdir) / "test.bam"
            bai_path = Path(tmpdir) / "test.bam.bai"
            
            # Create dummy files
            bam_path.touch()
            bai_path.touch()
            
            from amprenta_rag.ingestion.genomics.bam_parser import check_index_exists
            
            has_index, index_path = check_index_exists(bam_path)
            
            assert has_index is True
            assert index_path is not None

    def test_check_index_exists_returns_false_when_missing(self):
        """check_index_exists returns False when no index."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bam_path = Path(tmpdir) / "test.bam"
            bam_path.touch()
            
            from amprenta_rag.ingestion.genomics.bam_parser import check_index_exists
            
            has_index, index_path = check_index_exists(bam_path)
            
            assert has_index is False
            assert index_path is None

    @patch("amprenta_rag.ingestion.genomics.bam_parser.pysam")
    def test_parse_bam_header_extracts_references(self, mock_pysam):
        """parse_bam_header extracts reference sequences."""
        mock_bam = MagicMock()
        mock_bam.header.to_dict.return_value = {
            "SQ": [
                {"SN": "chr1", "LN": 249250621},
                {"SN": "chr2", "LN": 243199373},
            ],
            "RG": [{"ID": "sample1", "SM": "Sample1", "PL": "ILLUMINA"}],
        }
        mock_pysam.AlignmentFile.return_value.__enter__.return_value = mock_bam
        
        from amprenta_rag.ingestion.genomics.bam_parser import parse_bam_header
        
        result = parse_bam_header(Path("/fake/path.bam"))
        
        assert result["num_references"] == 2
        assert len(result["reference_sequences"]) == 2
        assert result["reference_sequences"][0]["name"] == "chr1"

    @patch("amprenta_rag.ingestion.genomics.bam_parser.pysam")
    def test_get_alignment_stats_calculates_metrics(self, mock_pysam):
        """get_alignment_stats calculates read statistics."""
        # Create mock reads
        mock_reads = []
        for i in range(100):
            read = MagicMock()
            read.is_unmapped = (i % 10 == 0)  # 10% unmapped
            read.is_duplicate = (i % 20 == 0)  # 5% duplicates
            mock_reads.append(read)
        
        mock_bam = MagicMock()
        mock_bam.fetch.return_value = iter(mock_reads)
        mock_pysam.AlignmentFile.return_value.__enter__.return_value = mock_bam
        
        from amprenta_rag.ingestion.genomics.bam_parser import get_alignment_stats
        
        stats = get_alignment_stats(Path("/fake/path.bam"))
        
        assert stats.total_reads == 100
        assert stats.unmapped_reads == 10
        assert stats.mapped_reads == 90
        assert stats.duplicate_rate == 0.05
