# ssw-aligner
# Extracted StripedSmithWaterman (SSW) C++ implementation and Python wrapper.
# Drop-in replacement for skbio.alignment.StripedSmithWaterman.

from ssw_aligner._wrapper import (
    AlignmentStructure,
    GumbelParams,
    SmithWatermanProfileAligner,
    StripedSmithWaterman,
    compute_gumbel_params,
    _MMSEQS_AA_FREQS as MMSEQS_AA_FREQS,
)

__all__ = [
    "AlignmentStructure",
    "GumbelParams",
    "MMSEQS_AA_FREQS",
    "SmithWatermanProfileAligner",
    "StripedSmithWaterman",
    "compute_gumbel_params",
]
__version__ = "0.2.0"
