"""Teams and Projects management page."""
import streamlit as st
from uuid import UUID
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Team, TeamMember, Project
from amprenta_rag.auth.session import get_current_user
from amprenta_rag.auth.permissions import get_user_teams, get_user_projects, get_team_role


def render_teams_page():
    st.title("👥 Teams & Projects")
    st.markdown("Manage teams and projects")

    tab1, tab2, tab3, tab4 = st.tabs(["My Teams", "My Projects", "Create Team", "Create Project"])

    with tab1:
        render_my_teams_tab()

    with tab2:
        render_my_projects_tab()

    with tab3:
        render_create_team_tab()

    with tab4:
        render_create_project_tab()


def render_my_teams_tab():
    st.subheader("My Teams")
    user = get_current_user()
    
    if not user:
        st.error("Please log in")
        return
    
    db_gen = get_db()
    db = next(db_gen)
    try:
        teams = get_user_teams(user.get("id"), db)
        
        if not teams:
            st.info("You are not a member of any teams yet.")
            return
        
        for team in teams:
            role = get_team_role(user.get("id"), str(team.id), db)
            role_badge = {"owner": "👑", "admin": "🔧", "member": "👤", "viewer": "👁️"}.get(role or "", "❓")
            
            with st.expander(f"{role_badge} **{team.name}** ({role or 'unknown'})"):
                st.markdown(f"**Description:** {team.description or 'No description'}")
                st.caption(f"Created: {team.created_at.strftime('%Y-%m-%d') if team.created_at else 'Unknown'}")
                
                # Show team projects
                team_projects = db.query(Project).filter(Project.team_id == team.id).all()
                if team_projects:
                    st.markdown(f"**Projects:** {len(team_projects)}")
                    for proj in team_projects:
                        st.markdown(f"- {proj.name} {'(Public)' if proj.is_public else ''}")
    finally:
        db_gen.close()


def render_my_projects_tab():
    st.subheader("My Projects")
    user = get_current_user()
    
    if not user:
        st.error("Please log in")
        return
    
    db_gen = get_db()
    db = next(db_gen)
    try:
        projects = get_user_projects(user.get("id"), db)
        
        if not projects:
            st.info("You don't have access to any projects yet.")
            return
        
        # Group by team
        projects_by_team = {}
        for proj in projects:
            team = db.query(Team).filter(Team.id == proj.team_id).first()
            team_name = team.name if team else "Unknown Team"
            if team_name not in projects_by_team:
                projects_by_team[team_name] = []
            projects_by_team[team_name].append(proj)
        
        for team_name, team_projects in projects_by_team.items():
            st.markdown(f"### {team_name}")
            for proj in team_projects:
                public_badge = "🌐 Public" if proj.is_public else "🔒 Private"
                with st.expander(f"{public_badge} **{proj.name}**"):
                    st.markdown(f"**Description:** {proj.description or 'No description'}")
                    st.caption(f"Created: {proj.created_at.strftime('%Y-%m-%d') if proj.created_at else 'Unknown'}")
    finally:
        db_gen.close()


def render_create_team_tab():
    st.subheader("Create Team")
    user = get_current_user()
    
    if not user:
        st.error("Please log in")
        return
    
    # Admin only
    if user.get("role") != "admin":
        st.error("Only administrators can create teams.")
        return
    
    with st.form("create_team"):
        name = st.text_input("Team Name*", max_chars=255)
        description = st.text_area("Description")
        
        submitted = st.form_submit_button("Create Team", type="primary")
        
        if submitted:
            if not name:
                st.error("Team name is required")
                return
            
            db_gen = get_db()
            db = next(db_gen)
            try:
                # Check for duplicate name
                existing = db.query(Team).filter(Team.name == name).first()
                if existing:
                    st.error("A team with this name already exists")
                    return
                
                new_team = Team(
                    name=name,
                    description=description if description else None,
                )
                db.add(new_team)
                db.commit()
                
                # Add creator as owner
                membership = TeamMember(
                    team_id=new_team.id,
                    user_id=user.get("id") if user.get("id") != "test" else None,
                    role="owner",
                )
                db.add(membership)
                db.commit()
                
                st.success(f"Team '{name}' created! You have been added as owner.")
                st.rerun()
            finally:
                db_gen.close()


def render_create_project_tab():
    st.subheader("Create Project")
    user = get_current_user()
    
    if not user:
        st.error("Please log in")
        return
    
    db_gen = get_db()
    db = next(db_gen)
    try:
        # Get teams user belongs to
        teams = get_user_teams(user.get("id"), db)
        
        if not teams:
            st.warning("You must be a member of a team to create a project. Create a team first or ask to be added to one.")
            return
        
        team_options = {t.name: t.id for t in teams}
        
        with st.form("create_project"):
            selected_team = st.selectbox("Team*", list(team_options.keys()))
            name = st.text_input("Project Name*", max_chars=255)
            description = st.text_area("Description")
            is_public = st.checkbox("Make project public", value=False)
            
            submitted = st.form_submit_button("Create Project", type="primary")
            
            if submitted:
                if not name:
                    st.error("Project name is required")
                    return
                
                team_id = team_options[selected_team]
                
                new_project = Project(
                    team_id=team_id,
                    name=name,
                    description=description if description else None,
                    is_public=is_public,
                )
                db.add(new_project)
                db.commit()
                
                st.success(f"Project '{name}' created!")
                st.rerun()
    finally:
        db_gen.close()


if __name__ == "__main__":
    render_teams_page()
