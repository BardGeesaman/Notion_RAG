"""AutoDock Vina execution helper (docker-first with local fallback)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Tuple


def _docker_available() -> bool:
    return shutil.which("docker") is not None


def _vina_available() -> bool:
    return shutil.which("vina") is not None


def run_vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    output_dir: str,
    *,
    exhaustiveness: int = 8,
    num_modes: int = 1,
) -> bool:
    """Run AutoDock Vina docking.

    Docker-first:
      docker run --rm -v {output_dir}:/data ccsb/autodock-vina vina ...

    Local fallback:
      vina ...

    Writes:
    - {output_dir}/out.pdbqt
    - {output_dir}/vina.log
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rec_path = Path(receptor_pdbqt)
    lig_path = Path(ligand_pdbqt)
    if not rec_path.exists():
        raise FileNotFoundError(str(rec_path))
    if not lig_path.exists():
        raise FileNotFoundError(str(lig_path))

    out_pdbqt = out_dir / "out.pdbqt"
    out_log = out_dir / "vina.log"

    cx, cy, cz = center
    sx, sy, sz = size

    # Docker path (copy inputs so container reads from /data)
    if _docker_available():
        (out_dir / "receptor.pdbqt").write_bytes(rec_path.read_bytes())
        (out_dir / "ligand.pdbqt").write_bytes(lig_path.read_bytes())
        cmd = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{str(out_dir)}:/data",
            "ccsb/autodock-vina",
            "vina",
            "--receptor",
            "/data/receptor.pdbqt",
            "--ligand",
            "/data/ligand.pdbqt",
            "--center_x",
            str(cx),
            "--center_y",
            str(cy),
            "--center_z",
            str(cz),
            "--size_x",
            str(sx),
            "--size_y",
            str(sy),
            "--size_z",
            str(sz),
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(num_modes),
            "--out",
            "/data/out.pdbqt",
            "--log",
            "/data/vina.log",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode == 0 and out_pdbqt.exists():
            return True

    # Local fallback
    if _vina_available():
        cmd = [
            "vina",
            "--receptor",
            str(rec_path),
            "--ligand",
            str(lig_path),
            "--center_x",
            str(cx),
            "--center_y",
            str(cy),
            "--center_z",
            str(cz),
            "--size_x",
            str(sx),
            "--size_y",
            str(sy),
            "--size_z",
            str(sz),
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(num_modes),
            "--out",
            str(out_pdbqt),
            "--log",
            str(out_log),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=str(out_dir))
        return proc.returncode == 0 and out_pdbqt.exists()

    return False


__all__ = ["run_vina"]


