"""Tests for Y.js WebSocket server deployment configuration."""

import json
import os
import pytest
from pathlib import Path
from unittest.mock import mock_open, patch

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent


class TestYjsServerDeployment:
    """Test Y.js WebSocket server deployment files."""

    def test_package_json_structure(self):
        """Test that package.json has correct structure and dependencies."""
        package_json_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "yjs-server" / "package.json"
        
        assert package_json_path.exists(), "package.json should exist"
        
        with open(package_json_path, 'r') as f:
            package_data = json.load(f)
        
        # Verify basic structure
        assert package_data["name"] == "yjs-websocket-server"
        assert "version" in package_data
        assert package_data["main"] == "server.js"
        
        # Verify required dependencies
        dependencies = package_data.get("dependencies", {})
        required_deps = ["yjs", "y-websocket", "ws", "pg", "express", "dotenv"]
        
        for dep in required_deps:
            assert dep in dependencies, f"Missing dependency: {dep}"
        
        # Verify Node.js version requirement
        assert "engines" in package_data
        assert "node" in package_data["engines"]
        assert "18" in package_data["engines"]["node"]

    def test_server_js_exists_and_executable(self):
        """Test that server.js exists and is properly configured."""
        server_js_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "yjs-server" / "server.js"
        
        assert server_js_path.exists(), "server.js should exist"
        
        with open(server_js_path, 'r') as f:
            content = f.read()
        
        # Check for essential components
        assert "#!/usr/bin/env node" in content
        assert "require('y-websocket')" in content or "y-websocket" in content
        assert "require('yjs')" in content or "yjs" in content
        assert "require('pg')" in content or "pg" in content
        assert "WebSocket" in content
        assert "/health" in content
        assert "PostgreSQL" in content or "postgres" in content
        
        # Check for proper server setup
        assert "PORT" in content
        assert "HOST" in content
        assert "setupWSConnection" in content or "WebSocket" in content

    def test_dockerfile_configuration(self):
        """Test that Dockerfile has correct configuration."""
        dockerfile_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "yjs-server" / "Dockerfile"
        
        assert dockerfile_path.exists(), "Dockerfile should exist"
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check base image
        assert "FROM node:" in content
        assert "alpine" in content.lower()
        
        # Check for security best practices
        assert "USER" in content  # Non-root user
        assert "HEALTHCHECK" in content
        assert "EXPOSE 1234" in content
        
        # Check for required commands
        assert "npm ci" in content or "npm install" in content
        assert "COPY" in content
        assert "WORKDIR" in content
        
        # Check environment setup
        assert "NODE_ENV" in content
        assert "CMD" in content or "ENTRYPOINT" in content

    def test_docker_compose_integration(self):
        """Test that docker-compose.yml includes yjs-server service."""
        compose_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "docker-compose.yml"
        
        assert compose_path.exists(), "docker-compose.yml should exist"
        
        with open(compose_path, 'r') as f:
            content = f.read()
        
        # Check for yjs-server service
        assert "yjs-server:" in content
        assert "amprenta-yjs-server" in content
        assert "1234:1234" in content
        
        # Check for proper network configuration
        assert "jupyterhub-network" in content
        
        # Check for PostgreSQL environment variables
        assert "POSTGRES_HOST" in content
        assert "POSTGRES_DB" in content
        assert "POSTGRES_USER" in content
        
        # Check for health check configuration
        assert "healthcheck:" in content
        assert "/health" in content
        
        # Check for restart policy
        assert "restart:" in content
        assert "unless-stopped" in content or "always" in content

    def test_yjs_server_environment_variables(self):
        """Test that server.js handles environment variables correctly."""
        server_js_path = PROJECT_ROOT / "deploy" / "jupyterhub" / "yjs-server" / "server.js"
        
        with open(server_js_path, 'r') as f:
            content = f.read()
        
        # Check for environment variable usage
        env_vars = [
            "process.env.PORT",
            "process.env.HOST", 
            "process.env.POSTGRES_HOST",
            "process.env.POSTGRES_PORT",
            "process.env.POSTGRES_DB",
            "process.env.POSTGRES_USER",
            "process.env.POSTGRES_PASSWORD"
        ]
        
        for env_var in env_vars:
            assert env_var in content, f"Missing environment variable: {env_var}"
        
        # Check for default values
        assert "|| 1234" in content  # Default port
        assert "|| '0.0.0.0'" in content  # Default host
        assert "|| 'postgres'" in content  # Default postgres host
        assert "|| 'amprenta'" in content  # Default database name
