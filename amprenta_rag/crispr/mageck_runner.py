"""MAGeCK execution helper (docker-first with local fallback)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def _docker_available() -> bool:
    return shutil.which("docker") is not None


def _mageck_available() -> bool:
    return shutil.which("mageck") is not None


def run_mageck_test(count_file: str, control: str, treatment: str, output_prefix: str) -> str:
    """Run `mageck test` and return the path to the gene summary output.

    Docker-first:
      docker run --rm -v {out_dir}:/data biocontainers/mageck mageck test ...

    Local fallback:
      mageck test ...

    MAGeCK typically writes:
    - {output_prefix}.gene_summary.txt
    - {output_prefix}.sgrna_summary.txt
    - {output_prefix}.log
    """
    count_path = Path(count_file)
    if not count_path.exists():
        raise FileNotFoundError(str(count_path))

    out_prefix = Path(output_prefix)
    out_dir = out_prefix.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Ensure count file is accessible inside container by copying into out_dir.
    count_copy = out_dir / "counts.txt"
    count_copy.write_bytes(count_path.read_bytes())

    gene_summary = out_dir / f"{out_prefix.name}.gene_summary.txt"

    def _build_args(count_in_container: str) -> list[str]:
        return [
            "mageck",
            "test",
            "-k",
            count_in_container,
            "-t",
            treatment,
            "-c",
            control,
            "-n",
            f"/data/{out_prefix.name}" if count_in_container.startswith("/data/") else str(out_prefix),
        ]

    # Docker path
    if _docker_available():
        cmd = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{str(out_dir)}:/data",
            "biocontainers/mageck",
        ] + _build_args("/data/counts.txt")
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode == 0 and gene_summary.exists():
            return str(gene_summary)

    # Local fallback
    if _mageck_available():
        cmd = _build_args(str(count_copy))
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=str(out_dir))
        if proc.returncode == 0 and gene_summary.exists():
            return str(gene_summary)

    raise RuntimeError("MAGeCK execution failed (docker and local mageck unavailable or returned error).")


__all__ = ["run_mageck_test"]


