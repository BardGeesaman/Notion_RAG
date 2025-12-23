"""PLIP execution helper (docker-first with local fallback)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def _docker_available() -> bool:
    return shutil.which("docker") is not None


def _plip_available() -> bool:
    return shutil.which("plip") is not None


def run_plip(complex_pdb: str, output_dir: str) -> bool:
    """Run PLIP on a protein-ligand complex PDB, producing XML output.

    Docker-first:
      docker run --rm -v {output_dir}:/data pharmai/plip -f /data/complex.pdb -x

    Local fallback:
      plip -f complex.pdb -x
    """
    in_path = Path(complex_pdb)
    if not in_path.exists():
        raise FileNotFoundError(str(in_path))

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Docker path: copy input into output_dir so container can read it.
    if _docker_available():
        (out_dir / "complex.pdb").write_bytes(in_path.read_bytes())
        cmd = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{str(out_dir)}:/data",
            "pharmai/plip",
            "-f",
            "/data/complex.pdb",
            "-x",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return proc.returncode == 0

    if _plip_available():
        # Run in output_dir with a local copy for consistent output placement
        local_in = out_dir / "complex.pdb"
        local_in.write_bytes(in_path.read_bytes())
        proc = subprocess.run(
            ["plip", "-f", str(local_in), "-x"],
            cwd=str(out_dir),
            capture_output=True,
            text=True,
            check=False,
        )
        return proc.returncode == 0

    return False


__all__ = ["run_plip"]


