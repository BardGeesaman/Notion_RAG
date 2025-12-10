import sqlite3

import numpy as np
import pytest

from amprenta_rag.chemistry.compound_linking import (
    compute_signature_reversal_score,
    find_compounds_reversing_signature,
    find_signatures_affected_by_compound,
    link_compound_to_program,
    link_compound_to_signature,
)
from amprenta_rag.chemistry.database import initialize_chemistry_database
from amprenta_rag.chemistry.schema import COMPOUNDS_TABLE


@pytest.fixture
def temp_db(tmp_path):
    db_path = tmp_path / "chem.db"
    initialize_chemistry_database(db_path=db_path)
    return db_path


def _insert_compound(conn: sqlite3.Connection, cid: str):
    conn.execute(
        """
        INSERT INTO compounds (compound_id, smiles, molecular_weight)
        VALUES (?, 'C', 100.0)
        """,
        (cid,),
    )
    conn.commit()


def test_link_compound_to_signature(temp_db):
    conn = sqlite3.connect(str(temp_db))
    _insert_compound(conn, "C1")
    conn.close()

    link_compound_to_signature("C1", "S1", effect_type="reverses", correlation=-0.7, p_value=0.01, db_path=temp_db)

    conn = sqlite3.connect(str(temp_db))
    row = conn.execute(
        "SELECT effect_type, correlation, p_value FROM compound_signature WHERE compound_id=? AND signature_id=?",
        ("C1", "S1"),
    ).fetchone()
    conn.close()

    assert row is not None
    assert row[0] == "reverses"
    assert pytest.approx(row[1]) == -0.7
    assert pytest.approx(row[2]) == 0.01


def test_find_compounds_reversing_signature(temp_db):
    conn = sqlite3.connect(str(temp_db))
    _insert_compound(conn, "C1")
    _insert_compound(conn, "C2")
    conn.execute(COMPOUNDS_TABLE)
    conn.commit()
    conn.close()

    link_compound_to_signature("C1", "S1", effect_type="reverses", correlation=-0.6, p_value=0.02, db_path=temp_db)
    link_compound_to_signature("C2", "S1", effect_type="mimics", correlation=0.4, p_value=0.05, db_path=temp_db)

    results = find_compounds_reversing_signature("S1", min_correlation=-0.5, db_path=temp_db)
    assert ("C1", -0.6, 0.02) in results
    assert all(r[0] != "C2" for r in results)


def test_find_signatures_affected_by_compound(temp_db):
    conn = sqlite3.connect(str(temp_db))
    _insert_compound(conn, "C1")
    conn.close()

    link_compound_to_signature("C1", "S1", effect_type="reverses", correlation=-0.6, db_path=temp_db)
    link_compound_to_signature("C1", "S2", effect_type="mimics", correlation=0.3, db_path=temp_db)

    links = find_signatures_affected_by_compound("C1", db_path=temp_db)
    sig_ids = {l.signature_id for l in links}
    assert sig_ids == {"S1", "S2"}


def test_compute_signature_reversal_score():
    compound_features = {"A": -1.0, "B": -0.5, "C": -0.2, "D": 0.1}
    signature_features = {"A": 1.0, "B": 0.6, "C": 0.3, "D": -0.1}

    corr, pval, effect = compute_signature_reversal_score(compound_features, signature_features)
    assert effect in ("reverses", "mimics")
    assert not np.isnan(corr)
    assert not np.isnan(pval)


def test_link_compound_to_program(temp_db):
    conn = sqlite3.connect(str(temp_db))
    _insert_compound(conn, "C1")
    conn.close()

    link_compound_to_program("C1", "P1", role="LEAD", db_path=temp_db)

    conn = sqlite3.connect(str(temp_db))
    row = conn.execute(
        "SELECT role FROM compound_program WHERE compound_id=? AND program_id=?",
        ("C1", "P1"),
    ).fetchone()
    conn.close()

    assert row is not None
    assert row[0] == "LEAD"

