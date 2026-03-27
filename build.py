"""Poetry build script — compiles libssw_aligner.so via CMake.

Invoked automatically by `poetry install` and `poetry build`.
The compiled library is placed inside the ssw_aligner/ package directory
so it is bundled into wheels and found by the ctypes loader at runtime.
"""
import multiprocessing
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.resolve()
BUILD_DIR = ROOT / "build"
DEST = ROOT / "ssw_aligner" / "libssw_aligner.so"


def build(setup_kwargs: dict) -> None:  # noqa: ARG001 — Poetry passes setup_kwargs
    cmake = shutil.which("cmake")
    if cmake is None:
        print("ERROR: cmake not found on PATH. Install CMake >= 3.14.", file=sys.stderr)
        sys.exit(1)

    jobs = str(multiprocessing.cpu_count())

    configure_cmd = [
        cmake,
        "-B", str(BUILD_DIR),
        "-S", str(ROOT),
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    build_cmd = [cmake, "--build", str(BUILD_DIR), "-j", jobs]

    print(f"[build.py] Configuring: {' '.join(configure_cmd)}")
    subprocess.run(configure_cmd, check=True)

    print(f"[build.py] Building:    {' '.join(build_cmd)}")
    subprocess.run(build_cmd, check=True)

    built_so = BUILD_DIR / "libssw_aligner.so"
    if not built_so.exists():
        print(f"ERROR: expected {built_so} after build.", file=sys.stderr)
        sys.exit(1)

    print(f"[build.py] Copying {built_so} -> {DEST}")
    shutil.copy2(built_so, DEST)
