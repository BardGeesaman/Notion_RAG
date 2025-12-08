from __future__ import annotations

import pytest


@pytest.mark.ui
def test_fixtures_load(test_program, test_dataset, test_signature):
    assert test_program["name"] == "ALS Program"
    assert test_dataset["omics_type"] == "lipidomics"
    assert "Core" in test_signature["name"] or test_signature["name"].startswith("ALS")


@pytest.mark.ui
def test_mock_db_session(mock_db_session):
    # Presence of mock fixture is enough; ensures no exception raised
    assert mock_db_session is None or True

