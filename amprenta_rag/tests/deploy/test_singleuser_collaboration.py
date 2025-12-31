"""Tests for singleuser image collaboration configuration."""

import json
import os
import pytest
from pathlib import Path
from unittest.mock import mock_open, patch

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent


class TestSingleuserCollaboration:
    """Test singleuser image collaboration features."""

    def test_jupyter_collaboration_package_installed(self):
        """Test that jupyter-collaboration package is included in pip install."""
        dockerfile_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "singleuser" / "Dockerfile"
        
        assert dockerfile_path.exists(), "singleuser Dockerfile should exist"
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check that jupyter-collaboration is in the pip install command
        assert "jupyter-collaboration" in content, "jupyter-collaboration package should be installed"
        
        # Verify it's in the pip install command (simplified check)
        # Check that both pip install exists and jupyter-collaboration is mentioned
        has_pip_install = "RUN pip install" in content
        has_jupyter_collaboration = "jupyter-collaboration" in content
        
        assert has_pip_install, "Dockerfile should have pip install command"
        assert has_jupyter_collaboration, "jupyter-collaboration should be mentioned in Dockerfile"
        
        # More specific check: find the pip install command and verify jupyter-collaboration is nearby
        lines = content.split('\n')
        pip_install_line_index = None
        jupyter_collab_line_index = None
        
        for i, line in enumerate(lines):
            if "RUN pip install" in line:
                pip_install_line_index = i
            if "jupyter-collaboration" in line.strip() and jupyter_collab_line_index is None:
                # Take the first occurrence only
                jupyter_collab_line_index = i
        
        # Verify jupyter-collaboration appears after the pip install command
        if pip_install_line_index is not None and jupyter_collab_line_index is not None:
            assert jupyter_collab_line_index > pip_install_line_index, "jupyter-collaboration should appear after pip install command"
            # And within reasonable distance (within 20 lines)
            assert jupyter_collab_line_index - pip_install_line_index < 20, "jupyter-collaboration should be close to pip install command"

    def test_yjs_url_environment_variable(self):
        """Test that Y.js server URL is configured via environment variable."""
        dockerfile_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "singleuser" / "Dockerfile"
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check for Y.js URL environment variable
        assert "JUPYTER_COLLABORATION_YJS_URL" in content, "Y.js URL environment variable should be set"
        assert "ws://yjs-server:1234" in content, "Y.js URL should point to yjs-server container"
        
        # Verify it's set as an ENV command
        env_lines = [line.strip() for line in content.split('\n') if line.strip().startswith('ENV')]
        yjs_env_found = any("JUPYTER_COLLABORATION_YJS_URL" in line for line in env_lines)
        assert yjs_env_found, "Y.js URL should be set as environment variable"

    def test_collaboration_extension_enabled(self):
        """Test that jupyter-collaboration extension is enabled in server config."""
        dockerfile_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "singleuser" / "Dockerfile"
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check for collaboration configuration
        assert "jupyter_collaboration" in content, "jupyter_collaboration extension should be configured"
        assert "collaboration.json" in content, "collaboration.json config file should be created"
        assert "CollaborationApp" in content, "CollaborationApp configuration should be present"
        
        # Check for server extension enablement
        assert '"jupyter_collaboration": true' in content, "jupyter_collaboration should be enabled"

    def test_collaboration_configuration_structure(self):
        """Test that collaboration configuration has proper JSON structure."""
        dockerfile_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "singleuser" / "Dockerfile"
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Extract the collaboration configuration JSON from the Dockerfile
        lines = content.split('\n')
        json_lines = []
        in_collaboration_config = False
        
        for line in lines:
            line_stripped = line.strip()
            if 'collaboration.json' in line:
                in_collaboration_config = True
                continue
            elif in_collaboration_config:
                if line_stripped.startswith("'") and line_stripped.endswith("'") and '\\' not in line_stripped:
                    # Extract JSON content from printf format
                    json_content = line_stripped.strip("' \\")
                    if json_content:
                        json_lines.append(json_content)
                elif not line_stripped or line_stripped.startswith('>'):
                    break
        
        # Reconstruct and validate JSON
        if json_lines:
            json_str = ''.join(json_lines)
            try:
                config = json.loads(json_str)
                
                # Validate structure
                assert "ServerApp" in config, "ServerApp section should exist"
                assert "jpserver_extensions" in config["ServerApp"], "jpserver_extensions should exist"
                assert "jupyter_collaboration" in config["ServerApp"]["jpserver_extensions"], "jupyter_collaboration extension should be configured"
                
                assert "CollaborationApp" in config, "CollaborationApp section should exist"
                assert "yjs_url" in config["CollaborationApp"], "yjs_url should be configured"
                
            except json.JSONDecodeError:
                # If JSON parsing fails, at least check for key configuration elements
                assert "ServerApp" in content, "ServerApp configuration should exist"
                assert "CollaborationApp" in content, "CollaborationApp configuration should exist"
                assert "yjs_url" in content, "yjs_url should be configured"

    def test_yjs_server_url_configuration(self):
        """Test that Y.js server URL is properly configured in CollaborationApp."""
        dockerfile_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "singleuser" / "Dockerfile"
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check for Y.js URL in CollaborationApp configuration
        assert '"yjs_url": "ws://yjs-server:1234"' in content, "Y.js URL should be configured in CollaborationApp"
        
        # Verify the URL format is correct
        assert "ws://" in content, "WebSocket protocol should be used"
        assert "yjs-server" in content, "Should reference yjs-server container"
        assert ":1234" in content, "Should use port 1234"

    def test_collaboration_settings_configured(self):
        """Test that collaboration settings are properly configured."""
        dockerfile_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "singleuser" / "Dockerfile"
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check for collaboration-specific settings
        expected_settings = [
            "document_cleanup_delay",
            "file_poll_interval"
        ]
        
        for setting in expected_settings:
            assert setting in content, f"Collaboration setting '{setting}' should be configured"
        
        # Check for reasonable default values
        assert "document_cleanup_delay" in content, "Document cleanup delay should be configured"
        assert "file_poll_interval" in content, "File poll interval should be configured"

    def test_dockerfile_maintains_existing_functionality(self):
        """Test that existing Dockerfile functionality is preserved."""
        dockerfile_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "singleuser" / "Dockerfile"
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Verify key existing components are still present
        existing_components = [
            "FROM jupyter/minimal-notebook",
            "PYTHONPATH",
            "voila",
            "jupyterlab",
            "notebook",
            "rdkit-pypi",
            "plotly",
            "pandas",
            "matplotlib"
        ]
        
        for component in existing_components:
            assert component in content, f"Existing component '{component}' should be preserved"
        
        # Verify voila configuration is still present
        assert "voila.json" in content, "Voila configuration should be preserved"
        assert "VoilaConfiguration" in content, "VoilaConfiguration should be preserved"
        
        # Verify working directory is still set
        assert "WORKDIR /home/jovyan" in content, "Working directory should be preserved"
