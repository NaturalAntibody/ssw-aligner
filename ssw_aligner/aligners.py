from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, Optional

from ssw_aligner._blosum_62 import BLOSUM_62

if TYPE_CHECKING:
    from ssw_aligner._ssw_wrapper import AlignmentStructure

ALIGNER_PARAMS: dict[str, Any] = {
    "match_score": 1,
    "mismatch_score": -1,
    "gap_open_penalty": 4,
    "gap_extend_penalty": 1,
}

AA_ALIGNER_PARAMS: dict[str, Any] = {
    "gap_open_penalty": 11,
    "gap_extend_penalty": 1,
    "protein": True,
    "substitution_matrix": BLOSUM_62,
}

_EXTENSION_IMPORT_MESSAGE = (
    "ssw_aligner compiled extension could not be imported. "
    "Install a wheel or rebuild the extension for the active Python environment."
)

_extension_import_error: Optional[BaseException]
try:
    from ssw_aligner._ssw_wrapper import StripedSmithWaterman as _StripedSmithWaterman
except ImportError as exc:
    _extension_import_error = exc
    _StripedSmithWaterman = None
else:
    _extension_import_error = None


class _AlignerWrapper:
    DEFAULTS: Mapping[str, object] = {}

    def __init__(self, query_sequence: str, **kwargs: Any) -> None:
        if _extension_import_error is not None or _StripedSmithWaterman is None:
            raise ImportError(_EXTENSION_IMPORT_MESSAGE) from _extension_import_error

        aligner_kwargs = dict(self.DEFAULTS)
        aligner_kwargs.update(kwargs)
        self._aligner = _StripedSmithWaterman(query_sequence, **aligner_kwargs)

    def __call__(self, target_sequence: str) -> "AlignmentStructure":
        return self._aligner(target_sequence)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._aligner, name)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._aligner!r})"


class NucleotideAligner(_AlignerWrapper):
    DEFAULTS = ALIGNER_PARAMS


class ProteinAligner(_AlignerWrapper):
    DEFAULTS = AA_ALIGNER_PARAMS
