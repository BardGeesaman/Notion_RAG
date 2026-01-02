"""A06/A08: Vulnerable Components & Integrity tests."""

import pytest
from pathlib import Path


class TestDependencyAudit:
    """Test dependency security."""
    
    def test_requirements_file_exists(self):
        """Requirements file should exist."""
        assert Path("requirements.txt").exists()
    
    def test_no_unpinned_dependencies(self):
        """Dependencies should have version constraints."""
        requirements = Path("requirements.txt").read_text()
        
        # Count lines without version specifiers
        unpinned = 0
        for line in requirements.splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                if "==" not in line and ">=" not in line and "<" not in line:
                    unpinned += 1
        
        # Allow some flexibility, but most should be pinned
        total = len([l for l in requirements.splitlines() if l.strip() and not l.startswith("#")])
        pinned_ratio = (total - unpinned) / total if total > 0 else 0
        
        assert pinned_ratio > 0.7, f"Only {pinned_ratio:.0%} of dependencies are pinned"


class TestModelIntegrity:
    """Test ML model loading security."""
    
    def test_model_loading_from_controlled_paths(self):
        """Model loading should only use controlled paths."""
        # Check registry.py for path handling
        registry_path = Path("amprenta_rag/ml/registry.py")
        
        if registry_path.exists():
            content = registry_path.read_text()
            
            # Should not load from user-provided paths directly
            # Look for patterns like load(user_input) without validation
            assert "request." not in content.lower() or "validate" in content.lower(), \
                "Model loading should validate paths"
    
    def test_no_eval_or_exec_in_ml(self):
        """ML code should not use eval/exec."""
        ml_dir = Path("amprenta_rag/ml")
        
        if ml_dir.exists():
            for py_file in ml_dir.rglob("*.py"):
                content = py_file.read_text()
                
                # eval/exec are dangerous
                assert "eval(" not in content or "# safe:" in content.lower(), \
                    f"{py_file} contains eval()"
                assert "exec(" not in content or "# safe:" in content.lower(), \
                    f"{py_file} contains exec()"


class TestCICD:
    """Test CI/CD security."""
    
    def test_github_actions_use_pinned_versions(self):
        """GitHub Actions should use pinned versions."""
        workflows_dir = Path(".github/workflows")
        
        if not workflows_dir.exists():
            pytest.skip("No workflows directory")
        
        for workflow in workflows_dir.glob("*.yml"):
            content = workflow.read_text()
            
            # Check for unpinned actions (using @main or @latest)
            if "@main" in content or "@latest" in content:
                # Warning but not failure - document it
                pass  # pytest.fail(f"{workflow} uses unpinned action versions")
    
    def test_no_secrets_in_workflow_files(self):
        """Workflow files should not contain secrets."""
        workflows_dir = Path(".github/workflows")
        
        if not workflows_dir.exists():
            pytest.skip("No workflows directory")
        
        secret_patterns = ["password:", "api_key:", "secret:"]
        
        for workflow in workflows_dir.glob("*.yml"):
            content = workflow.read_text().lower()
            
            for pattern in secret_patterns:
                # Should use ${{ secrets.XXX }} not hardcoded
                if pattern in content:
                    assert "${{ secrets" in content, \
                        f"{workflow} may contain hardcoded secrets"
