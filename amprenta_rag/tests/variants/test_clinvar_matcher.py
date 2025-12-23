from __future__ import annotations

from uuid import UUID

from amprenta_rag.database.models import Variant
from amprenta_rag.variants.clinvar_matcher import match_variants_to_clinvar


def test_match_variants_to_clinvar_matches_on_chr_pos_ref_alt():
    v = Variant(
        variant_set_id=UUID("00000000-0000-0000-0000-000000000001"),
        chromosome="chr1",
        position=123,
        ref_allele="a",
        alt_allele="t",
    )
    v.id = UUID("00000000-0000-0000-0000-000000000010")

    lookup = {
        ("1", 123, "A", "T"): {
            "clinvar_id": "42",
            "clinical_significance": "Pathogenic",
            "review_status": "criteria provided, single submitter",
            "condition": "Example condition",
        }
    }

    anns = match_variants_to_clinvar([v], lookup)
    assert len(anns) == 1
    a = anns[0]
    assert a.variant_id == v.id
    assert a.annotation_source == "clinvar"
    assert a.clinvar_id == "42"
    assert a.clinical_significance == "Pathogenic"


