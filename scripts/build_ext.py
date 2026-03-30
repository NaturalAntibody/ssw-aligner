"""Poetry build script — compiles libssw_aligner.so via CMake.

Invoked automatically by ``poetry install``, ``poetry build``, and any
PEP 517 frontend (``python -m build``, ``pip install``).

poetry-core runs this script as ``python scripts/build_ext.py`` (NOT as a
function call), so the actual work happens at the module level.
The compiled library is placed inside the ssw_aligner/ package directory
so it is bundled into wheels and found by the ctypes loader at runtime.
"""
import multiprocessing
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
BUILD_DIR = ROOT / "build"
DEST = ROOT / "ssw_aligner" / "libssw_aligner.so"


def _compile() -> None:
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

    print(f"[build_ext] Configuring: {' '.join(configure_cmd)}")
    subprocess.run(configure_cmd, check=True)

    print(f"[build_ext] Building:    {' '.join(build_cmd)}")
    subprocess.run(build_cmd, check=True)

    built_so = BUILD_DIR / "libssw_aligner.so"
    if not built_so.exists():
        print(f"ERROR: expected {built_so} after build.", file=sys.stderr)
        sys.exit(1)

    print(f"[build_ext] Copying {built_so} -> {DEST}")
    shutil.copy2(built_so, DEST)


# poetry-core runs the script directly (not as a module import),
# so call _compile() unconditionally at the top level.
_compile()
