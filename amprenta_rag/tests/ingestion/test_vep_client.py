"""Tests for VEP REST API client."""

from unittest.mock import patch, MagicMock
from amprenta_rag.ingestion.genomics.vep_client import (
    annotate_variant,
    annotate_variants_batch,
    VEPAnnotation,
    _parse_vep_response,
)


class TestVEPAnnotation:
    """Test VEPAnnotation Pydantic model."""
    
    def test_vep_annotation_model_creation(self):
        """Test VEPAnnotation model with all fields."""
        annotation = VEPAnnotation(
            consequence="missense_variant",
            impact="HIGH",
            symbol="TP53",
            gene_id="ENSG00000141510",
            sift_prediction="deleterious",
            sift_score=0.01,
            polyphen_prediction="probably_damaging",
            polyphen_score=0.95,
        )
        assert annotation.impact == "HIGH"
        assert annotation.sift_score == 0.01
    
    def test_vep_annotation_optional_fields(self):
        """Test VEPAnnotation with minimal fields."""
        annotation = VEPAnnotation(consequence="synonymous_variant")
        assert annotation.consequence == "synonymous_variant"
        assert annotation.impact is None


class TestAnnotateVariant:
    """Test single variant annotation."""
    
    @patch("amprenta_rag.ingestion.genomics.vep_client._vep_request")
    def test_annotate_variant_success(self, mock_vep_request):
        """Test successful single variant annotation."""
        mock_response = MagicMock()
        mock_response.json.return_value = [{
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{
                "impact": "MODERATE",
                "gene_symbol": "BRCA1",
                "sift_prediction": "tolerated",
            }]
        }]
        mock_vep_request.return_value = mock_response
        
        result = annotate_variant("17", 43082403, "G", "A")
        
        assert result is not None
        assert result.consequence == "missense_variant"
        assert result.impact == "MODERATE"
        assert result.symbol == "BRCA1"
    
    @patch("amprenta_rag.ingestion.genomics.vep_client._vep_request")
    def test_annotate_variant_api_error(self, mock_vep_request):
        """Test handling of API errors."""
        mock_vep_request.side_effect = Exception("API Error")
        
        result = annotate_variant("1", 12345, "A", "T")
        
        assert result is None


class TestAnnotateVariantsBatch:
    """Test batch variant annotation."""
    
    @patch("amprenta_rag.ingestion.genomics.vep_client._vep_request")
    def test_batch_annotation_success(self, mock_vep_request):
        """Test successful batch annotation."""
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {"most_severe_consequence": "missense_variant", "transcript_consequences": [{"impact": "HIGH"}]},
            {"most_severe_consequence": "synonymous_variant", "transcript_consequences": [{"impact": "LOW"}]},
        ]
        mock_vep_request.return_value = mock_response
        
        variants = [
            {"chromosome": "1", "position": 100, "ref": "A", "alt": "T"},
            {"chromosome": "2", "position": 200, "ref": "C", "alt": "G"},
        ]
        
        results = annotate_variants_batch(variants)
        
        assert len(results) == 2
        assert results[0].impact == "HIGH"
        assert results[1].impact == "LOW"


class TestParseVEPResponse:
    """Test VEP response parsing."""
    
    def test_parse_with_transcript_consequences(self):
        """Test parsing response with transcript consequences."""
        data = {
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{
                "impact": "MODERATE",
                "gene_symbol": "TP53",
                "gene_id": "ENSG00000141510",
                "transcript_id": "ENST00000269305",
                "biotype": "protein_coding",
                "amino_acids": "R/H",
                "protein_start": 273,
                "sift_prediction": "deleterious",
                "sift_score": 0.01,
                "polyphen_prediction": "probably_damaging",
                "polyphen_score": 0.999,
            }]
        }
        
        result = _parse_vep_response(data)
        
        assert result.consequence == "missense_variant"
        assert result.symbol == "TP53"
        assert result.sift_score == 0.01
    
    def test_parse_intergenic_variant(self):
        """Test parsing intergenic variant."""
        data = {
            "most_severe_consequence": "intergenic_variant",
            "intergenic_consequences": [{"impact": "MODIFIER"}]
        }
        
        result = _parse_vep_response(data)
        
        assert result.consequence == "intergenic_variant"