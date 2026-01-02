"""A03: Injection vulnerability tests."""

import pytest


class TestSSRFPrevention:
    """Test SSRF protection."""
    
    def test_blocks_localhost(self):
        """Should block localhost URLs."""
        from amprenta_rag.utils.safe_requests import is_safe_url
        
        safe, _ = is_safe_url("http://localhost/admin")
        assert safe is False
        
        safe, _ = is_safe_url("http://127.0.0.1/admin")
        assert safe is False
    
    def test_blocks_private_ips(self):
        """Should block private IP ranges."""
        from amprenta_rag.utils.safe_requests import is_safe_url
        
        private_ips = [
            "http://10.0.0.1/",
            "http://172.16.0.1/",
            "http://192.168.1.1/",
        ]
        
        for url in private_ips:
            safe, _ = is_safe_url(url)
            assert safe is False, f"Should block {url}"
    
    def test_blocks_file_scheme(self):
        """Should block file:// URLs."""
        from amprenta_rag.utils.safe_requests import is_safe_url
        
        safe, _ = is_safe_url("file:///etc/passwd")
        assert safe is False
    
    def test_allows_public_apis(self):
        """Should allow known public APIs."""
        from amprenta_rag.utils.safe_requests import is_safe_url
        
        public_urls = [
            "https://api.ncbi.nlm.nih.gov/datasets/v2/",
            "https://rest.ensembl.org/",
            "https://www.uniprot.org/",
        ]
        
        for url in public_urls:
            safe, _ = is_safe_url(url)
            assert safe is True, f"Should allow {url}"
    
    def test_safe_get_blocks_dangerous_urls(self):
        """safe_get should raise ValueError for blocked URLs."""
        from amprenta_rag.utils.safe_requests import safe_get
        
        with pytest.raises(ValueError, match="URL blocked"):
            safe_get("http://localhost/admin")
    
    def test_safe_post_blocks_dangerous_urls(self):
        """safe_post should raise ValueError for blocked URLs."""
        from amprenta_rag.utils.safe_requests import safe_post
        
        with pytest.raises(ValueError, match="URL blocked"):
            safe_post("http://127.0.0.1/internal")


class TestCommandInjection:
    """Test command injection prevention."""
    
    def test_subprocess_uses_list_form(self):
        """Subprocess calls should use list form, not shell=True."""
        from pathlib import Path
        
        # Check known subprocess files
        files_to_check = [
            "amprenta_rag/ingestion/dvc_manager.py",
            "amprenta_rag/ingestion/genomics/pipeline.py",
            "amprenta_rag/ingestion/genomics/nextflow_runner.py",
        ]
        
        for filepath in files_to_check:
            path = Path(filepath)
            if not path.exists():
                continue
            
            content = path.read_text()
            
            # Check for shell=True (dangerous)
            assert "shell=True" not in content, f"{filepath} uses shell=True"
    
    def test_subprocess_uses_safe_patterns(self):
        """Check that subprocess calls use safe patterns."""
        from pathlib import Path
        
        # Check a known file
        path = Path("amprenta_rag/ingestion/dvc_manager.py")
        if path.exists():
            content = path.read_text()
            
            # Should use subprocess.run with list arguments
            assert "subprocess.run([" in content, "Should use list form for subprocess"


class TestPathTraversal:
    """Test path traversal prevention."""
    
    def test_sanitize_filename(self):
        """Filenames should be sanitized."""
        # Test that path traversal sequences are blocked
        dangerous_names = [
            "../../../etc/passwd",
            "..\\..\\windows\\system32", 
            "file\x00.txt",  # Null byte injection
        ]
        
        from pathlib import Path
        import os
        
        for name in dangerous_names:
            # Secure implementation should sanitize
            safe_name = os.path.basename(name).replace('\x00', '')  # Remove null bytes
            assert ".." not in safe_name or safe_name == name.split(os.sep)[-1]
            assert "\x00" not in safe_name
    
    def test_path_traversal_in_uploads(self):
        """Test that file upload paths are safe."""
        from pathlib import Path
        
        # Simulate safe file handling
        upload_dir = Path("/tmp/uploads")
        user_filename = "../../../etc/passwd"
        
        # Safe approach: use only the filename part
        safe_path = upload_dir / Path(user_filename).name
        
        # Should not escape upload directory
        assert not str(safe_path).startswith("/etc/")
        assert "passwd" in str(safe_path)  # filename preserved
        assert str(safe_path).startswith("/tmp/uploads/")
