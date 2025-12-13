from unittest.mock import MagicMock

from amprenta_rag.auth import permissions as perms


USER_ID = "00000000-0000-0000-0000-000000000001"
PROJECT_ID = "00000000-0000-0000-0000-000000000002"
TEAM_ID = "00000000-0000-0000-0000-000000000003"


def _make_db_with_queries(project_obj=None, membership_obj=None, memberships_list=None):
    """
    Helper to build a db mock with separate query objects for Project and TeamMember.
    """
    db = MagicMock()
    q_project = MagicMock()
    q_member = MagicMock()

    def _query_side_effect(model):
        if model is perms.Project:
            return q_project
        if model is perms.TeamMember:
            return q_member
        # Used only by get_user_teams -> TeamMember, but keep safe fallback
        return MagicMock()

    db.query.side_effect = _query_side_effect

    if memberships_list is not None:
        q_member.filter.return_value.all.return_value = memberships_list

    q_project.filter.return_value.first.return_value = project_obj
    q_member.filter.return_value.first.return_value = membership_obj

    return db, q_project, q_member


def test_get_user_teams():
    team1 = MagicMock(spec=perms.Team)
    team2 = MagicMock(spec=perms.Team)
    m1 = MagicMock()
    m1.team = team1
    m2 = MagicMock()
    m2.team = team2

    db, _, q_member = _make_db_with_queries(memberships_list=[m1, m2])
    teams = perms.get_user_teams(USER_ID, db)

    assert teams == [team1, team2]
    assert q_member.filter.called is True


def test_can_view_public_project():
    project = MagicMock()
    project.is_public = True
    project.team_id = perms.UUID(TEAM_ID)

    db, _, q_member = _make_db_with_queries(project_obj=project, membership_obj=None)
    assert perms.can_view_project(USER_ID, PROJECT_ID, db) is True
    # For public project, membership check should not be required
    assert q_member.filter.called is False


def test_can_view_private_project_as_member():
    project = MagicMock()
    project.is_public = False
    project.team_id = perms.UUID(TEAM_ID)

    membership = MagicMock()
    db, _, q_member = _make_db_with_queries(project_obj=project, membership_obj=membership)
    assert perms.can_view_project(USER_ID, PROJECT_ID, db) is True
    assert q_member.filter.called is True


def test_cannot_view_private_project_as_nonmember():
    project = MagicMock()
    project.is_public = False
    project.team_id = perms.UUID(TEAM_ID)

    db, _, q_member = _make_db_with_queries(project_obj=project, membership_obj=None)
    assert perms.can_view_project(USER_ID, PROJECT_ID, db) is False
    assert q_member.filter.called is True


def test_can_edit_project_as_admin():
    project = MagicMock()
    project.team_id = perms.UUID(TEAM_ID)

    membership = MagicMock()
    membership.role = "admin"

    db, _, q_member = _make_db_with_queries(project_obj=project, membership_obj=membership)
    assert perms.can_edit_project(USER_ID, PROJECT_ID, db) is True
    assert q_member.filter.called is True


def test_cannot_edit_project_as_viewer():
    project = MagicMock()
    project.team_id = perms.UUID(TEAM_ID)

    db, _, q_member = _make_db_with_queries(project_obj=project, membership_obj=None)
    assert perms.can_edit_project(USER_ID, PROJECT_ID, db) is False
    assert q_member.filter.called is True


