"""Tests for JupyterHub collaboration configuration.

Validates that JupyterHub is properly configured to support real-time collaborative
editing via Y.js server integration.
"""

import os
import tempfile
import pytest
from unittest.mock import patch, MagicMock


class TestHubCollaborationConfig:
    """Test JupyterHub configuration for collaborative editing support."""

    @pytest.fixture
    def jupyterhub_config_path(self):
        """Fixture providing path to JupyterHub configuration file."""
        return "deploy/jupyterhub/jupyterhub_config.py"

    @pytest.fixture
    def config_content(self, jupyterhub_config_path):
        """Fixture providing JupyterHub configuration content."""
        with open(jupyterhub_config_path, 'r') as f:
            return f.read()

    def test_yjs_environment_variable_configured(self, config_content):
        """Test that Y.js server URL is configured in spawner environment."""
        # Verify JUPYTER_COLLABORATION_YJS_URL is in environment configuration
        assert "JUPYTER_COLLABORATION_YJS_URL" in config_content
        assert "yjs-server:1234" in config_content
        
        # Verify it uses environment variable with fallback
        assert 'os.environ.get("YJS_SERVER_URL", "ws://yjs-server:1234")' in config_content

    def test_docker_spawner_network_configured(self, config_content):
        """Test that DockerSpawner network is configured for Y.js connectivity."""
        # Verify network configuration allows container-to-container communication
        assert 'c.DockerSpawner.network_name = "jupyterhub-network"' in config_content
        
        # Verify spawner is using DockerSpawner
        assert 'c.JupyterHub.spawner_class = "dockerspawner.DockerSpawner"' in config_content

    def test_environment_variables_structure(self, config_content):
        """Test that environment variables are properly structured."""
        # Verify environment dictionary contains required variables
        assert "c.DockerSpawner.environment = {" in config_content
        assert '"API_URL"' in config_content
        assert '"JUPYTERHUB_SINGLEUSER_APP"' in config_content
        assert '"JUPYTER_COLLABORATION_YJS_URL"' in config_content

    @patch.dict(os.environ, {"YJS_SERVER_URL": "ws://custom-yjs:9999"})
    def test_custom_yjs_url_from_environment(self, jupyterhub_config_path):
        """Test that custom Y.js URL can be set via environment variable."""
        # Create a temporary config file to test environment variable substitution
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as temp_config:
            with open(jupyterhub_config_path, 'r') as original:
                content = original.read()
                # Remove the problematic import for testing
                content = content.replace("from token_authenticator import TokenAuthenticator", "# TokenAuthenticator import mocked")
                content = content.replace("c.JupyterHub.authenticator_class = TokenAuthenticator", "# c.JupyterHub.authenticator_class = TokenAuthenticator")
                temp_config.write(content)
            temp_config_path = temp_config.name

        try:
            # Execute the config file to simulate JupyterHub loading it
            config_globals = {"os": os, "c": MagicMock(), "sys": MagicMock()}
            with open(temp_config_path, 'r') as f:
                config_code = f.read()
            
            exec(config_code, config_globals)
            
            # Verify the custom Y.js URL was used
            spawner_env = config_globals["c"].DockerSpawner.environment
            assert spawner_env["JUPYTER_COLLABORATION_YJS_URL"] == "ws://custom-yjs:9999"
            
        finally:
            os.unlink(temp_config_path)

    def test_default_yjs_url_fallback(self, jupyterhub_config_path):
        """Test that default Y.js URL is used when environment variable not set."""
        # Create a temporary config file without the problematic import
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as temp_config:
            with open(jupyterhub_config_path, 'r') as original:
                content = original.read()
                # Remove the problematic import for testing
                content = content.replace("from token_authenticator import TokenAuthenticator", "# TokenAuthenticator import mocked")
                content = content.replace("c.JupyterHub.authenticator_class = TokenAuthenticator", "# c.JupyterHub.authenticator_class = TokenAuthenticator")
                temp_config.write(content)
            temp_config_path = temp_config.name
        
        try:
            # Ensure YJS_SERVER_URL is not in environment
            with patch.dict(os.environ, {}, clear=False):
                if "YJS_SERVER_URL" in os.environ:
                    del os.environ["YJS_SERVER_URL"]
                
                # Execute the config file
                config_globals = {"os": os, "c": MagicMock(), "sys": MagicMock()}
                with open(temp_config_path, 'r') as f:
                    config_code = f.read()
                
                exec(config_code, config_globals)
                
                # Verify the default Y.js URL was used
                spawner_env = config_globals["c"].DockerSpawner.environment
                assert spawner_env["JUPYTER_COLLABORATION_YJS_URL"] == "ws://yjs-server:1234"
                
        finally:
            os.unlink(temp_config_path)

    def test_collaboration_config_integration(self, config_content):
        """Test that collaboration configuration integrates properly with existing config."""
        # Verify existing configuration is preserved
        assert "API_URL" in config_content
        assert "JUPYTERHUB_SINGLEUSER_APP" in config_content
        
        # Verify new collaboration config doesn't break existing functionality
        assert "jupyter_server.serverapp.ServerApp" in config_content
        assert "host.docker.internal:8000" in config_content
        
        # Verify collaboration config is properly integrated
        assert "# Y.js WebSocket server URL for real-time collaborative editing" in config_content

    def test_spawner_configuration_completeness(self, config_content):
        """Test that spawner configuration is complete for collaboration support."""
        # Verify all required spawner settings for collaboration
        required_configs = [
            'c.JupyterHub.spawner_class = "dockerspawner.DockerSpawner"',
            'c.DockerSpawner.image = "amprenta-singleuser:latest"',
            'c.DockerSpawner.network_name = "jupyterhub-network"',
            'c.DockerSpawner.environment = {',
        ]
        
        for config_line in required_configs:
            assert config_line in config_content, f"Missing required config: {config_line}"

    def test_network_connectivity_requirements(self, config_content):
        """Test that network configuration meets collaboration requirements."""
        # Verify network name allows container-to-container communication
        assert '"jupyterhub-network"' in config_content
        
        # Verify the network configuration is applied to DockerSpawner
        assert "c.DockerSpawner.network_name" in config_content
        
        # Ensure the Y.js server URL uses the correct protocol and port
        yjs_url_pattern = "ws://yjs-server:1234"
        assert yjs_url_pattern in config_content
