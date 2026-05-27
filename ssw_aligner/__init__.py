# ssw-aligner
# Extracted StripedSmithWaterman (SSW) C implementation and Python/Cython wrapper.
# Drop-in replacement for skbio.alignment.StripedSmithWaterman.

from ssw_aligner._blosum_62 import BLOSUM_62
from ssw_aligner.aligners import ProteinAligner, NucleotideAligner
from ssw_aligner.metrics import (
	GumbellParams,
	calculate_seq_identity,
	compute_bit_score,
	compute_evalue,
	compute_raw_score_aa,
	compute_raw_score_from_bit_score,
)

_extension_import_error = None

try:
	from ssw_aligner._ssw_wrapper import AlignmentStructure, StripedSmithWaterman
except ImportError as exc:
	_extension_import_error = exc

__all__ = [
	"AlignmentStructure",
	"BLOSUM_62",
	"GumbellParams",
	"NucleotideAligner",
	"ProteinAligner",
	"StripedSmithWaterman",
	"calculate_seq_identity",
	"compute_bit_score",
	"compute_evalue",
	"compute_raw_score_aa",
	"compute_raw_score_from_bit_score",
]
__version__ = "0.1.0"


def __getattr__(name: str):
	if name in {"AlignmentStructure", "StripedSmithWaterman"} and _extension_import_error is not None:
		raise ImportError(
			"ssw_aligner compiled extension could not be imported. "
			"Install a wheel or rebuild the extension for the active Python environment."
		) from _extension_import_error
	raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
