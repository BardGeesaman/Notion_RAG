"""Project Workspaces Hub Dashboard."""

import streamlit as st
import httpx
import os
from typing import Dict, Optional


def render_workspaces_page():
    """Render the Project Workspaces Hub page."""
    st.title("🏢 Project Workspaces")
    st.markdown("Collaborate on projects and access shared resources")
    
    # Create 3-column layout
    col1, col2, col3 = st.columns([1, 2, 1])
    
    with col1:
        st.subheader("Teams")
        _render_teams_section()
    
    with col2:
        st.subheader("Projects")
        selected_team = st.session_state.get("selected_team_id")
        if selected_team:
            _render_projects_section(selected_team)
        else:
            st.info("Select a team to view projects")
    
    with col3:
        st.subheader("Shared With Me")
        _render_shared_entities_section()


def _render_teams_section():
    """Render the teams section in the left sidebar."""
    try:
        # Fetch user's teams
        teams_response = _api_get("/api/v1/teams/my-teams")
        
        if teams_response and teams_response.get("teams"):
            teams = teams_response["teams"]
            
            # Team selection
            team_options = ["Select a team..."] + [
                f"{team['name']} ({team['member_count']} members)" 
                for team in teams
            ]
            
            selected_index = st.selectbox(
                "Choose team:",
                range(len(team_options)),
                format_func=lambda x: team_options[x],
                key="team_selector"
            )
            
            if selected_index > 0:
                selected_team = teams[selected_index - 1]
                st.session_state.selected_team_id = selected_team["id"]
                st.session_state.selected_team_name = selected_team["name"]
                
                # Show team info
                st.info(f"**{selected_team['name']}**\n\n"
                       f"Role: {selected_team['user_role']}\n"
                       f"Members: {selected_team['member_count']}")
                
                if selected_team.get("description"):
                    st.caption(selected_team["description"])
            else:
                st.session_state.selected_team_id = None
                st.session_state.selected_team_name = None
        else:
            st.info("You're not a member of any teams yet.")
    
    except Exception as e:
        st.error(f"Failed to load teams: {str(e)}")
    
    # Create Team section
    st.markdown("---")
    with st.expander("➕ Create New Team"):
        with st.form("create_team_form"):
            team_name = st.text_input("Team Name", placeholder="Enter team name")
            team_description = st.text_area(
                "Description (optional)", 
                placeholder="Brief description of the team"
            )
            
            if st.form_submit_button("Create Team"):
                if team_name.strip():
                    _create_team(team_name.strip(), team_description.strip())
                else:
                    st.error("Team name is required")


def _render_projects_section(team_id: str):
    """Render the projects section for the selected team."""
    try:
        # Fetch projects for the team
        projects_response = _api_get(f"/api/v1/projects?team_id={team_id}")
        
        if projects_response and projects_response.get("projects"):
            projects = projects_response["projects"]
            
            # Display projects in a grid
            cols_per_row = 2
            for i in range(0, len(projects), cols_per_row):
                cols = st.columns(cols_per_row)
                
                for j, project in enumerate(projects[i:i+cols_per_row]):
                    with cols[j]:
                        _render_project_card(project)
        else:
            st.info("No projects found for this team.")
            
    except Exception as e:
        st.error(f"Failed to load projects: {str(e)}")
    
    # Create Project button
    st.markdown("---")
    with st.expander("➕ Create New Project"):
        with st.form("create_project_form"):
            project_name = st.text_input("Project Name", placeholder="Enter project name")
            project_description = st.text_area(
                "Description", 
                placeholder="Brief description of the project"
            )
            is_public = st.checkbox("Public project", value=False)
            
            if st.form_submit_button("Create Project"):
                if project_name.strip():
                    _create_project(team_id, project_name.strip(), project_description.strip(), is_public)
                else:
                    st.error("Project name is required")


def _render_project_card(project: Dict):
    """Render a project card."""
    with st.container():
        st.markdown(
            f"""
            <div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin: 8px 0;">
                <h4 style="margin: 0 0 8px 0;">{project['name']}</h4>
                <p style="margin: 0 0 12px 0; color: #666; font-size: 14px;">
                    {project.get('description', 'No description')}
                </p>
                <div style="font-size: 12px; color: #888;">
                    <span>📊 Datasets: {project.get('dataset_count', 0)}</span><br>
                    <span>🧪 Experiments: {project.get('experiment_count', 0)}</span><br>
                    <span>🔬 Compounds: {project.get('compound_count', 0)}</span>
                </div>
            </div>
            """,
            unsafe_allow_html=True
        )
        
        if st.button("Open Project", key=f"open_project_{project['id']}"):
            st.switch_page("pages/Projects.py")  # Navigate to project detail


def _render_shared_entities_section():
    """Render the shared entities section."""
    try:
        # Fetch entities shared with current user
        shares_response = _api_get("/api/v1/sharing/my-shares")
        
        if shares_response and shares_response.get("shares"):
            shares = shares_response["shares"]
            
            # Group shares by entity type
            grouped_shares = {}
            for share in shares:
                entity_type = share["entity_type"]
                if entity_type not in grouped_shares:
                    grouped_shares[entity_type] = []
                grouped_shares[entity_type].append(share)
            
            # Display grouped shares
            for entity_type, type_shares in grouped_shares.items():
                st.markdown(f"**{entity_type.title()}s ({len(type_shares)})**")
                
                for share in type_shares[:5]:  # Show max 5 per type
                    _render_shared_entity_item(share)
                
                if len(type_shares) > 5:
                    st.caption(f"... and {len(type_shares) - 5} more")
                
                st.markdown("---")
        else:
            st.info("No entities shared with you yet.")
            
    except Exception as e:
        st.error(f"Failed to load shared entities: {str(e)}")


def _render_shared_entity_item(share: Dict):
    """Render a shared entity item."""
    entity_type = share["entity_type"]
    entity_id = share["entity_id"]
    permission = share["permission"]
    
    # Get entity name (simplified - in real app would fetch from respective API)
    entity_name = f"{entity_type}_{str(entity_id)[:8]}"
    
    # Permission color coding
    permission_colors = {
        "view": "🟢",
        "edit": "🟡", 
        "admin": "🔴"
    }
    
    permission_icon = permission_colors.get(permission, "⚪")
    
    with st.container():
        col_name, col_perm, col_btn = st.columns([3, 1, 1])
        
        with col_name:
            st.markdown(f"**{entity_name}**")
            st.caption(f"Shared by: User {str(share.get('shared_by_id', 'Unknown'))[:8]}")
        
        with col_perm:
            st.markdown(f"{permission_icon} {permission}")
        
        with col_btn:
            if st.button("Open", key=f"open_shared_{entity_id}"):
                _navigate_to_entity(entity_type, entity_id)


def _navigate_to_entity(entity_type: str, entity_id: str):
    """Navigate to the appropriate page for the entity type."""
    navigation_map = {
        "dataset": "pages/Datasets.py",
        "experiment": "pages/Experiments.py", 
        "compound": "pages/Chemistry.py",
        "signature": "pages/Signatures.py"
    }
    
    page = navigation_map.get(entity_type)
    if page:
        # Store entity ID in session state for the target page
        st.session_state[f"selected_{entity_type}_id"] = entity_id
        st.switch_page(page)
    else:
        st.error(f"Unknown entity type: {entity_type}")


def _create_team(name: str, description: str):
    """Create a new team."""
    try:
        data = {
            "name": name,
            "description": description if description else None
        }
        
        response = _api_post("/api/v1/teams/teams", data)
        
        if response:
            st.success(f"Team '{name}' created successfully!")
            st.rerun()  # Refresh the page
        else:
            st.error("Failed to create team")
            
    except Exception as e:
        st.error(f"Error creating team: {str(e)}")


def _create_project(team_id: str, name: str, description: str, is_public: bool):
    """Create a new project."""
    try:
        data = {
            "team_id": team_id,
            "name": name,
            "description": description if description else None,
            "is_public": is_public
        }
        
        response = _api_post("/api/v1/projects", data)
        
        if response:
            st.success(f"Project '{name}' created successfully!")
            st.rerun()  # Refresh the page
        else:
            st.error("Failed to create project")
            
    except Exception as e:
        st.error(f"Error creating project: {str(e)}")


def _api_get(endpoint: str) -> Optional[Dict]:
    """Make GET request to API."""
    try:
        api_url = os.getenv("API_URL", "http://localhost:8000")
        
        with httpx.Client() as client:
            response = client.get(
                f"{api_url}{endpoint}",
                headers={"X-User-Id": "mock-user-id"},  # TODO: Replace with real auth
                timeout=30.0
            )
            
            if response.status_code == 200:
                return response.json()
            else:
                st.error(f"API Error: {response.status_code}")
                return None
                
    except Exception as e:
        st.error(f"Network error: {str(e)}")
        return None


def _api_post(endpoint: str, data: Dict) -> Optional[Dict]:
    """Make POST request to API."""
    try:
        api_url = os.getenv("API_URL", "http://localhost:8000")
        
        with httpx.Client() as client:
            response = client.post(
                f"{api_url}{endpoint}",
                json=data,
                headers={"X-User-Id": "mock-user-id"},  # TODO: Replace with real auth
                timeout=30.0
            )
            
            if response.status_code in [200, 201]:
                return response.json()
            else:
                st.error(f"API Error: {response.status_code} - {response.text}")
                return None
                
    except Exception as e:
        st.error(f"Network error: {str(e)}")
        return None


if __name__ == "__main__":
    render_workspaces_page()
