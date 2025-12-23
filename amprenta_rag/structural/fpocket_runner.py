"""fpocket execution helper (docker-first with local fallback)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def _docker_available() -> bool:
    return shutil.which("docker") is not None


def _fpocket_available() -> bool:
    return shutil.which("fpocket") is not None


def run_fpocket(pdb_path: str, output_dir: str) -> bool:
    """Run fpocket on a PDB file.

    Tries:
    1) Docker: `docker run --rm -v {dir}:/data fpocket:latest fpocket -f /data/input.pdb`
    2) Local binary: `fpocket -f {pdb_path}`

    Notes:
    - When using docker, we copy input into output_dir as input.pdb so results land there.
    - When using local fpocket, it writes next to input file (fpocket output folder). We run in cwd=output_dir.
    """
    in_path = Path(pdb_path)
    if not in_path.exists():
        raise FileNotFoundError(str(in_path))
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Docker path
    if _docker_available():
        input_copy = out_dir / "input.pdb"
        input_copy.write_bytes(in_path.read_bytes())
        cmd = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{str(out_dir)}:/data",
            "fpocket:latest",
            "fpocket",
            "-f",
            "/data/input.pdb",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode == 0:
            return True

    # Local fallback
    if _fpocket_available():
        proc = subprocess.run(["fpocket", "-f", str(in_path)], cwd=str(out_dir), capture_output=True, text=True, check=False)
        return proc.returncode == 0

    return False


__all__ = ["run_fpocket"]


