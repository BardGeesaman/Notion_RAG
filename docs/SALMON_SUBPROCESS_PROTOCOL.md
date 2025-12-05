# Salmon Subprocess Protocol

## Critical Protocol Rules

This document describes the **STRICT Subprocess Protocol** for using Salmon in our genomics pipeline.

### 1. NO PYTHON IMPORTS

**⚠️ CRITICAL WARNING:**

- **DO NOT** write `import salmon`
- **DO NOT** attempt to install it via `pip install salmon`
- **Reason:** The package on PyPI named `salmon` is a mail server utility, NOT the bioinformatics quantifier. Installing it will break the environment.

**Assumption:** The user has already installed the Salmon binary via Conda (`conda install -c bioconda salmon`) or Homebrew. It is available in the system PATH.

### 2. THE SUBPROCESS PROTOCOL

Instead of importing a library, we execute Salmon as a **command-line subprocess** from within Python:

- Use `subprocess.run(..., check=True)` to execute the shell command
- Always define command arguments as a list of strings to avoid shell injection issues
- Check installation before running commands

### 3. IMPLEMENTATION PATTERN

All Salmon operations follow this exact structure:

```python
import subprocess

def check_salmon_installed() -> bool:
    """Verifies that the Salmon binary is accessible."""
    try:
        result = subprocess.run(
            ["salmon", "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        version_output = result.stdout.strip()
        logger.info("✅ Found Salmon: %s", version_output)
        return True
    except FileNotFoundError:
        logger.error(
            "❌ Error: 'salmon' command not found. "
            "Please install it via Conda: 'conda install -c bioconda salmon'"
        )
        return False

def run_salmon_command(...):
    """Run Salmon via subprocess."""
    if not check_salmon_installed():
        return None
    
    cmd = [
        "salmon",
        "quant",  # or "index", etc.
        "-i", index_path,
        "-l", "A",  # Auto-detect library type
        "-r", fastq_file,
        "-o", output_dir,
        "--validateMappings",  # Recommended flag
        "--quiet",  # Reduce console spam
    ]
    
    subprocess.run(cmd, check=True, capture_output=True, text=True)
```

### 4. RECOMMENDED FLAGS

When running Salmon quantification:

- `--validateMappings`: Recommended flag for accuracy
- `--quiet`: Reduce console spam
- `-l A`: Auto-detect library type (or specify: `U`=unstranded, `SR`=single-end reverse, etc.)

### 5. FILE LOCATIONS

Our implementation:

- **Pipeline Module:** `amprenta_rag/ingestion/genomics/pipeline.py`
  - `check_salmon_installed()`: Installation verification
  - `quantify_with_salmon()`: Quantification function
  - `build_salmon_index()`: Index building (in setup script)

- **Setup Script:** `scripts/setup_salmon_index.py`
  - Downloads reference transcriptome
  - Builds Salmon index

### 6. ERROR HANDLING

Always handle:

1. **FileNotFoundError**: Salmon not installed
2. **subprocess.CalledProcessError**: Command failed (check return code, stdout, stderr)
3. **File existence checks**: Verify input files and output files exist

### 7. VERIFICATION

After installation, verify:

```bash
which salmon
salmon --version
```

Expected output:
```
/opt/anaconda3/bin/salmon
salmon 1.10.3
```

---

**Last Updated:** Based on Gemini agent recommendations (Subprocess Protocol)

