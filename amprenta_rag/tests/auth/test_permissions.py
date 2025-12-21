from __future__ import annotations
import pytest
from uuid import uuid4, UUID
from amprenta_rag.auth import permissions as perms

# Mock Models
class Team:
    def __init__(self, id, name="Team"):
        self.id = id
        self.name = name

class TeamMember:
    def __init__(self, user_id, team_id, role="member", team=None):
        self.user_id = user_id
        self.team_id = team_id
        self.role = role
        self.team = team

class Project:
    def __init__(self, id, team_id, is_public=False):
        self.id = id
        self.team_id = team_id
        self.is_public = is_public

class FakeQuery:
    def __init__(self, results=None):
        self.results = results if results is not None else []
        self.filter_called = False

    def filter(self, *args, **kwargs):
        self.filter_called = True
        return self

    def all(self):
        return self.results

    def first(self):
        return self.results[0] if self.results else None

class FakeDB:
    def __init__(self):
        self.queries = {}

    def query(self, model):
        return self.queries.get(model, FakeQuery())

    def register_query(self, model, results):
        self.queries[model] = FakeQuery(results)

def test_get_user_teams():
    uid = uuid4()
    tid = uuid4()
    team = Team(tid)
    tm = TeamMember(uid, tid, team=team)
    
    db = FakeDB()
    db.register_query(perms.TeamMember, [tm])
    
    teams = perms.get_user_teams(str(uid), db)
    assert len(teams) == 1
    assert teams[0].id == tid

def test_get_user_teams_error():
    class BrokenDB:
        def query(self, model): raise RuntimeError("Fail")
        
    teams = perms.get_user_teams(str(uuid4()), BrokenDB())
    assert teams == []

def test_get_user_projects():
    uid = uuid4()
    tid = uuid4()
    pid1 = uuid4() # Team project
    pid2 = uuid4() # Public project
    
    tm = TeamMember(uid, tid)
    p1 = Project(pid1, tid, is_public=False)
    p2 = Project(pid2, uuid4(), is_public=True)
    
    db = FakeDB()
    db.register_query(perms.TeamMember, [tm])
    # Mocking chained queries is tricky with simple fakes, 
    # but the logic relies on .in_() for team projects.
    # We can simulate the results directly for the two Project queries (team + public).
    
    # This is a limitation of simple FakeDB. For this test, let's patch the db.query call
    # to return specific results based on call order or inspection, or just simplify.
    
    # Let's mock the query calls specifically for this flow
    class SmartDB:
        def __init__(self):
            self.call_count = 0
            
        def query(self, model):
            if model == perms.TeamMember:
                return FakeQuery([tm])
            if model == perms.Project:
                self.call_count += 1
                if self.call_count == 1: # Team projects
                    return FakeQuery([p1])
                else: # Public projects
                    return FakeQuery([p2])
            return FakeQuery()
            
    projects = perms.get_user_projects(str(uid), SmartDB())
    assert len(projects) == 2
    ids = {p.id for p in projects}
    assert pid1 in ids
    assert pid2 in ids

def test_can_view_project_public():
    pid = uuid4()
    p = Project(pid, uuid4(), is_public=True)
    db = FakeDB()
    db.register_query(perms.Project, [p])
    
    assert perms.can_view_project(str(uuid4()), str(pid), db) is True

def test_can_view_project_private_member():
    uid = uuid4()
    tid = uuid4()
    pid = uuid4()
    p = Project(pid, tid, is_public=False)
    tm = TeamMember(uid, tid)
    
    class SmartDB:
        def query(self, model):
            if model == perms.Project:
                return FakeQuery([p])
            if model == perms.TeamMember:
                return FakeQuery([tm])
            return FakeQuery()
            
    assert perms.can_view_project(str(uid), str(pid), SmartDB()) is True

def test_can_view_project_private_non_member():
    uid = uuid4()
    pid = uuid4()
    p = Project(pid, uuid4(), is_public=False)
    
    class SmartDB:
        def query(self, model):
            if model == perms.Project:
                return FakeQuery([p])
            if model == perms.TeamMember:
                return FakeQuery([]) # No membership found
            return FakeQuery()
            
    assert perms.can_view_project(str(uid), str(pid), SmartDB()) is False

def test_can_edit_project_owner():
    uid = uuid4()
    tid = uuid4()
    pid = uuid4()
    p = Project(pid, tid)
    tm = TeamMember(uid, tid, role="owner")
    
    class SmartDB:
        def query(self, model):
            if model == perms.Project:
                return FakeQuery([p])
            if model == perms.TeamMember:
                return FakeQuery([tm])
            return FakeQuery()
            
    assert perms.can_edit_project(str(uid), str(pid), SmartDB()) is True

def test_can_edit_project_member():
    uid = uuid4()
    tid = uuid4()
    pid = uuid4()
    p = Project(pid, tid)
    tm = TeamMember(uid, tid, role="member") # Not owner/admin
    
    class SmartDB:
        def query(self, model):
            if model == perms.Project:
                return FakeQuery([p])
            if model == perms.TeamMember:
                return FakeQuery([]) # Filter for owner/admin returns empty
            return FakeQuery()
            
    assert perms.can_edit_project(str(uid), str(pid), SmartDB()) is False

def test_get_team_role():
    uid = uuid4()
    tid = uuid4()
    tm = TeamMember(uid, tid, role="admin")
    
    db = FakeDB()
    db.register_query(perms.TeamMember, [tm])
    
    assert perms.get_team_role(str(uid), str(tid), db) == "admin"

def test_get_team_role_none():
    db = FakeDB() # Returns empty
    assert perms.get_team_role(str(uuid4()), str(uuid4()), db) is None
