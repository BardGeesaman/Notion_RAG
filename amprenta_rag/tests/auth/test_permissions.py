from __future__ import annotations

from dataclasses import dataclass
from uuid import uuid4


from amprenta_rag.auth import permissions as perm

@dataclass
class _Membership:
    team_id: object
    team: object = None
    role: str = "member"


@dataclass
class _Project:
    id: object
    team_id: object
    is_public: bool = False


class _Query:
    def __init__(self, all_rows=None, first_row=None):
        self._all_rows = list(all_rows or [])
        self._first_row = first_row

    def filter(self, *_args, **_kwargs):
        return self

    def all(self):
        return list(self._all_rows)

    def first(self):
        return self._first_row


class _DB:
    """
    Deterministic fake DB:
    - TeamMember queries: .all() returns memberships; .first() returns membership_first
    - Project queries: successive query(Project) calls return from project_queries list
    - Project first() (for can_view/can_edit) comes from project_first
    """

    def __init__(
        self,
        memberships=None,
        team_member_first=None,
        project_queries=None,
        project_first=None,
        raise_on_query: bool = False,
    ):
        self._memberships = memberships or []
        self._team_member_first = team_member_first
        self._project_queries = list(project_queries or [])
        self._project_first = project_first
        self._raise_on_query = raise_on_query

    def query(self, model):
        if self._raise_on_query:
            raise RuntimeError("Boom")

        if model is perm.TeamMember:
            return _Query(all_rows=self._memberships, first_row=self._team_member_first)

        if model is perm.Project:
            # If a "first" is needed (view/edit), perm calls .first() immediately.
            if self._project_first is not None:
                return _Query(first_row=self._project_first)
            # Otherwise (get_user_projects), perm calls .all(); pop successive lists.
            rows = self._project_queries.pop(0) if self._project_queries else []
            return _Query(all_rows=rows)

        # Unused in this module
        return _Query()


def test_get_user_teams_success():
    t1 = object()
    db = _DB(memberships=[_Membership(team_id=uuid4(), team=t1)])
    teams = perm.get_user_teams(str(uuid4()), db)
    assert teams == [t1]


def test_get_user_teams_error():
    db = _DB(raise_on_query=True)
    assert perm.get_user_teams(str(uuid4()), db) == []


def test_get_user_projects_team_and_public_deduped():
    team_id = uuid4()
    memberships = [_Membership(team_id=team_id)]
    shared_id = uuid4()
    team_projects = [_Project(id=shared_id, team_id=team_id, is_public=False)]
    public_projects = [_Project(id=shared_id, team_id=team_id, is_public=True)]

    # Two Project queries: team_projects then public_projects.
    db = _DB(memberships=memberships, project_queries=[team_projects, public_projects])
    projects = perm.get_user_projects(str(uuid4()), db)
    assert len(projects) == 1
    assert projects[0].id == shared_id


def test_get_user_projects_only_public_when_no_team_ids():
    public_projects = [_Project(id=uuid4(), team_id=uuid4(), is_public=True)]
    # Only one Project query should happen (public projects).
    db = _DB(memberships=[], project_queries=[public_projects])
    projects = perm.get_user_projects(str(uuid4()), db)
    assert [p.id for p in projects] == [public_projects[0].id]


def test_get_user_projects_error():
    db = _DB(raise_on_query=True)
    assert perm.get_user_projects(str(uuid4()), db) == []


def test_can_view_project_not_found():
    db = _DB(project_first=None)
    assert perm.can_view_project(str(uuid4()), str(uuid4()), db) is False


def test_can_view_project_public():
    proj = _Project(id=uuid4(), team_id=uuid4(), is_public=True)
    db = _DB(project_first=proj)
    assert perm.can_view_project(str(uuid4()), str(proj.id), db) is True


def test_can_view_project_private_member_vs_denied():
    proj = _Project(id=uuid4(), team_id=uuid4(), is_public=False)

    db_yes = _DB(project_first=proj, team_member_first=_Membership(team_id=proj.team_id))
    assert perm.can_view_project(str(uuid4()), str(proj.id), db_yes) is True

    db_no = _DB(project_first=proj, team_member_first=None)
    assert perm.can_view_project(str(uuid4()), str(proj.id), db_no) is False


def test_can_view_project_error():
    db = _DB(raise_on_query=True)
    assert perm.can_view_project(str(uuid4()), str(uuid4()), db) is False


def test_can_edit_project_owner_admin_only():
    proj = _Project(id=uuid4(), team_id=uuid4(), is_public=False)

    db_yes = _DB(project_first=proj, team_member_first=_Membership(team_id=proj.team_id, role="owner"))
    assert perm.can_edit_project(str(uuid4()), str(proj.id), db_yes) is True

    db_no = _DB(project_first=proj, team_member_first=None)
    assert perm.can_edit_project(str(uuid4()), str(proj.id), db_no) is False


def test_can_edit_project_not_found():
    db = _DB(project_first=None)
    assert perm.can_edit_project(str(uuid4()), str(uuid4()), db) is False


def test_can_edit_project_error():
    db = _DB(raise_on_query=True)
    assert perm.can_edit_project(str(uuid4()), str(uuid4()), db) is False


def test_get_team_role_success_none_error():
    team_id = uuid4()
    db_yes = _DB(team_member_first=_Membership(team_id=team_id, role="admin"))
    assert perm.get_team_role(str(uuid4()), str(team_id), db_yes) == "admin"

    db_none = _DB(team_member_first=None)
    assert perm.get_team_role(str(uuid4()), str(uuid4()), db_none) is None

    db_err = _DB(raise_on_query=True)
    assert perm.get_team_role(str(uuid4()), str(uuid4()), db_err) is None
