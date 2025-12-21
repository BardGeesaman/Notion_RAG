from amprenta_rag.domain.chemistry import (
    BiochemicalResultRow,
    ChemistryCampaignMetadata,
    ChemistryHitRow,
)


def test_chemistry_campaign_metadata_defaults():
    meta = ChemistryCampaignMetadata(campaign_id="c1", name="Campaign", description=None)
    assert meta.campaign_id == "c1"
    assert meta.metadata is None


def test_chemistry_hit_row_allows_optional():
    hit = ChemistryHitRow(hit_id="h1", compound_id="cmp", score=0.5)
    assert hit.score == 0.5
    assert hit.metadata is None


def test_biochemical_result_row_accepts_any_value():
    row = BiochemicalResultRow(result_id="r1", compound_id="cmp", campaign_id="c1", value={"x": 1})
    assert row.value == {"x": 1}

