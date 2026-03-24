# ssw-aligner
# Extracted StripedSmithWaterman (SSW) C implementation and Python/Cython wrapper.
# Drop-in replacement for skbio.alignment.StripedSmithWaterman.

from ssw_aligner._ssw_wrapper import AlignmentStructure, StripedSmithWaterman

__all__ = ["AlignmentStructure", "StripedSmithWaterman"]
__version__ = "0.1.0"
