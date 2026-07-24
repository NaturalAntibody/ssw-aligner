# ssw-aligner
# Extracted StripedSmithWaterman (SSW) C implementation and Python/Cython wrapper.
# Drop-in replacement for skbio.alignment.StripedSmithWaterman.

from typing import Any

from ssw_aligner._blosum_62 import BLOSUM_62
from ssw_aligner.aligners import ProteinAligner, NucleotideAligner

_extension_import_error: BaseException | None = None

try:
	from ssw_aligner._ssw_wrapper import AlignmentStructure, StripedSmithWaterman
except ImportError as exc:
	_extension_import_error = exc

__all__ = [
	"AlignmentStructure",
	"BLOSUM_62",
	"NucleotideAligner",
	"ProteinAligner",
	"StripedSmithWaterman",
]
__version__ = "0.1.0"


def __getattr__(name: str) -> Any:
	if name in {"AlignmentStructure", "StripedSmithWaterman"} and _extension_import_error is not None:
		raise ImportError(
			"ssw_aligner compiled extension could not be imported. "
			"Install a wheel or rebuild the extension for the active Python environment."
		) from _extension_import_error
	raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
