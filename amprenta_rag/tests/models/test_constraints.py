import uuid

import pytest
from sqlalchemy.exc import IntegrityError

from amprenta_rag.database.base import get_session_local
from amprenta_rag.database.models import Compound, Program, Team, User


def test_user_email_unique():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        email = f"{uuid.uuid4().hex[:10]}@example.com"
        u1 = User(
            id=uuid.uuid4(),
            username=f"user_{uuid.uuid4().hex[:8]}",
            email=email,
            password_hash="x",
        )
        u2 = User(
            id=uuid.uuid4(),
            username=f"user_{uuid.uuid4().hex[:8]}",
            email=email,  # duplicate
            password_hash="x",
        )

        db.add(u1)
        db.flush()

        db.add(u2)
        with pytest.raises(IntegrityError):
            db.flush()
    finally:
        db.rollback()
        db.close()


def test_compound_compound_id_unique():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        cid = f"CMP-{uuid.uuid4().hex[:10]}"
        c1 = Compound(id=uuid.uuid4(), compound_id=cid, smiles="CCO")
        c2 = Compound(id=uuid.uuid4(), compound_id=cid, smiles="CCO")  # duplicate

        db.add(c1)
        db.flush()

        db.add(c2)
        with pytest.raises(IntegrityError):
            db.flush()
    finally:
        db.rollback()
        db.close()


def test_team_name_not_null():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        t = Team(id=uuid.uuid4(), name=None)  # type: ignore[arg-type]
        db.add(t)
        with pytest.raises(IntegrityError):
            db.flush()
    finally:
        db.rollback()
        db.close()


def test_nullable_fields_accept_none():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        program = Program(
            id=uuid.uuid4(),
            name=f"Program {uuid.uuid4().hex[:8]}",
            description=None,
            disease=None,
            external_ids=None,
        )
        compound = Compound(
            id=uuid.uuid4(),
            compound_id=f"CMP-{uuid.uuid4().hex[:10]}",
            smiles="CCO",
            inchi_key=None,
            molecular_weight=None,
            external_ids=None,
        )

        db.add_all([program, compound])
        db.flush()
    finally:
        db.rollback()
        db.close()


