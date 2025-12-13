import uuid

from amprenta_rag.database.base import get_session_local
from amprenta_rag.database.models import (
    Compound,
    Experiment,
    HTSCampaign,
    HTSResult,
    LabNotebookEntry,
    LabNotebookEntryAssociation,
    Program,
    Project,
    Team,
    TeamMember,
    User,
)


def test_user_team_relationship():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        user = User(
            id=uuid.uuid4(),
            username=f"user_{uuid.uuid4().hex[:8]}",
            email=f"{uuid.uuid4().hex[:8]}@example.com",
            password_hash="x",
        )
        team = Team(
            id=uuid.uuid4(),
            name=f"Team {uuid.uuid4().hex[:8]}",
        )
        membership = TeamMember(
            id=uuid.uuid4(),
            user_id=user.id,
            team_id=team.id,
            role="member",
        )

        db.add_all([user, team, membership])
        db.flush()

        assert membership.user == user
        assert membership.team == team
        assert membership in team.members
        assert membership in user.team_memberships

        db.rollback()
    finally:
        db.close()


def test_compound_program_many_to_many():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        program = Program(
            id=uuid.uuid4(),
            name=f"Program {uuid.uuid4().hex[:8]}",
        )
        compound = Compound(
            id=uuid.uuid4(),
            compound_id=f"CMP-{uuid.uuid4().hex[:10]}",
            smiles="CCO",
        )

        db.add_all([program, compound])
        db.flush()

        compound.programs.append(program)
        db.flush()

        assert program in compound.programs
        assert compound in program.compounds

        db.rollback()
    finally:
        db.close()


def test_experiment_project_relationship():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        team = Team(
            id=uuid.uuid4(),
            name=f"Team {uuid.uuid4().hex[:8]}",
        )
        project = Project(
            id=uuid.uuid4(),
            team_id=team.id,
            name=f"Project {uuid.uuid4().hex[:8]}",
            is_public=False,
        )
        experiment = Experiment(
            id=uuid.uuid4(),
            name=f"Experiment {uuid.uuid4().hex[:8]}",
            project_id=project.id,
        )

        db.add_all([team, project, experiment])
        db.flush()

        assert experiment.project == project
        assert experiment in project.experiments

        db.rollback()
    finally:
        db.close()


def test_hts_result_compound_relationship():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        compound = Compound(
            id=uuid.uuid4(),
            compound_id=f"CMP-{uuid.uuid4().hex[:10]}",
            smiles="CCO",
        )
        campaign = HTSCampaign(
            id=uuid.uuid4(),
            campaign_id=f"HTS-{uuid.uuid4().hex[:10]}",
            campaign_name="Test Campaign",
        )
        result = HTSResult(
            id=uuid.uuid4(),
            result_id=f"RES-{uuid.uuid4().hex[:10]}",
            campaign_id=campaign.id,
            compound_id=compound.id,
            hit_flag=True,
        )

        db.add_all([compound, campaign, result])
        db.flush()

        assert result.compound == compound
        assert result in compound.hts_results
        assert result.campaign == campaign
        assert result in campaign.results

        db.rollback()
    finally:
        db.close()


def test_lab_notebook_associations():
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        entry = LabNotebookEntry(
            id=uuid.uuid4(),
            title="Test Entry",
            content="Hello",
        )
        db.add(entry)
        db.flush()

        assoc = LabNotebookEntryAssociation(
            entry_id=entry.id,
            entity_type="experiment",
            entity_id=uuid.uuid4(),
        )
        db.add(assoc)
        db.flush()

        assert assoc.entry == entry
        assert assoc in entry.linked_entities

        db.rollback()
    finally:
        db.close()


