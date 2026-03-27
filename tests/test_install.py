"""Integration tests — install ssw-aligner from wheel and sdist in clean venvs.

These tests build a wheel and an sdist (if not already present in dist/),
create a fresh virtual environment, install from the artifact, and exercise
the public API on example sequences.

Requirements to run:
    pytest (in current env)
    cmake (on PATH, for sdist build)
    python3 -m venv (stdlib)

Invocation:
    pytest tests/test_install.py -v --timeout=300
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DIST = ROOT / "dist"
SMOKE_SCRIPT = Path(__file__).resolve().parent / "smoke_test.py"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_dist_if_needed():
    """Ensure dist/ contains a current wheel and sdist."""
    if not DIST.exists() or not any(DIST.glob("*.whl")) or not any(DIST.glob("*.tar.gz")):
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "poetry-core>=2.0.0", "build"],
            check=True,
            capture_output=True,
        )
        subprocess.run(
            [sys.executable, "-m", "build", str(ROOT)],
            check=True,
            capture_output=True,
        )


def _latest_artifact(pattern: str) -> Path:
    """Return the newest file in dist/ matching *pattern*."""
    files = sorted(DIST.glob(pattern), key=os.path.getmtime, reverse=True)
    if not files:
        pytest.skip(f"No {pattern} found in {DIST}")
    return files[0]


def _run_in_clean_venv(artifact: Path) -> dict:
    """Create a temp venv, install *artifact* + numpy, run the smoke test."""
    import tempfile

    with tempfile.TemporaryDirectory(prefix="ssw_test_") as tmpdir:
        venv_dir = Path(tmpdir) / "venv"
        subprocess.run(
            [sys.executable, "-m", "venv", str(venv_dir)],
            check=True,
            capture_output=True,
        )
        pip = str(venv_dir / "bin" / "pip")
        python = str(venv_dir / "bin" / "python")

        # Upgrade pip so it handles modern metadata
        subprocess.run(
            [pip, "install", "--upgrade", "pip"],
            check=True,
            capture_output=True,
        )

        # Install the artifact (wheel or sdist)
        subprocess.run(
            [pip, "install", str(artifact)],
            check=True,
            capture_output=True,
        )

        # Copy and run the smoke-test script
        script = Path(tmpdir) / "smoke_test.py"
        shutil.copy2(SMOKE_SCRIPT, script)

        result = subprocess.run(
            [python, str(script)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            pytest.fail(
                f"Smoke test failed (exit {result.returncode}).\n"
                f"stdout: {result.stdout}\nstderr: {result.stderr}"
            )
        return json.loads(result.stdout.strip())


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestWheelInstall:
    """Install from a pre-built wheel and verify the API works."""

    @pytest.fixture(autouse=True, scope="class")
    def _ensure_dist(self):
        _build_dist_if_needed()

    @pytest.fixture(scope="class")
    def smoke_results(self):
        whl = _latest_artifact("*.whl")
        return _run_in_clean_venv(whl)

    def test_nucleotide_score_positive(self, smoke_results):
        assert smoke_results["nt_score"] > 0

    def test_nucleotide_has_cigar(self, smoke_results):
        assert smoke_results["nt_cigar"]

    def test_nucleotide_has_aligned_sequences(self, smoke_results):
        assert smoke_results["nt_has_aligned"] is True

    def test_protein_score_positive(self, smoke_results):
        assert smoke_results["prot_score"] > 0

    def test_profile_scores_match_seq_seq(self, smoke_results):
        assert smoke_results["profile_scores_match"] is True

    def test_mmseqs_aa_freqs_available(self, smoke_results):
        assert smoke_results["freqs_count"] == 20

    def test_version(self, smoke_results):
        assert smoke_results["version"] == "0.2.0"


class TestSdistInstall:
    """Install from an sdist (triggers CMake build) and verify the API works."""

    @pytest.fixture(autouse=True, scope="class")
    def _ensure_dist(self):
        _build_dist_if_needed()

    @pytest.fixture(scope="class")
    def smoke_results(self):
        sdist = _latest_artifact("*.tar.gz")
        return _run_in_clean_venv(sdist)

    def test_nucleotide_score_positive(self, smoke_results):
        assert smoke_results["nt_score"] > 0

    def test_nucleotide_has_cigar(self, smoke_results):
        assert smoke_results["nt_cigar"]

    def test_nucleotide_has_aligned_sequences(self, smoke_results):
        assert smoke_results["nt_has_aligned"] is True

    def test_protein_score_positive(self, smoke_results):
        assert smoke_results["prot_score"] > 0

    def test_profile_scores_match_seq_seq(self, smoke_results):
        assert smoke_results["profile_scores_match"] is True

    def test_mmseqs_aa_freqs_available(self, smoke_results):
        assert smoke_results["freqs_count"] == 20

    def test_version(self, smoke_results):
        assert smoke_results["version"] == "0.2.0"
