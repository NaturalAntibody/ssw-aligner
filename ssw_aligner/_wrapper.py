"""
Pure-Python wrapper around the MMseqs2-based SSW shared library.

Provides :class:`StripedSmithWaterman` and :class:`AlignmentStructure` that are
API-compatible with the ``skbio.alignment`` equivalents (and with the old
Cython-based ``ssw_aligner``).

The alignment engine lives in ``libssw_aligner.so`` which exposes four
``extern "C"`` entry points:

    ssw_init, ssw_align, ssw_free_cigar, ssw_destroy
"""

from __future__ import annotations

import ctypes
import ctypes.util
import os
import pathlib
from typing import Dict, Optional

import numpy as np

# ---------------------------------------------------------------------------
# Locate and load the shared library
# ---------------------------------------------------------------------------

def _find_library() -> ctypes.CDLL:
    """Locate libssw_aligner.so from known build paths."""
    # Common locations, relative to this file or a build/ directory.
    here = pathlib.Path(__file__).resolve().parent
    candidates = [
        here / "libssw_aligner.so",
        here.parent / "build" / "libssw_aligner.so",
        here.parent / "build" / "release" / "libssw_aligner.so",
    ]
    # Also honour an env-var override.
    env = os.environ.get("SSW_ALIGNER_LIB")
    if env:
        candidates.insert(0, pathlib.Path(env))

    for p in candidates:
        if p.exists():
            return ctypes.CDLL(str(p))

    # Last resort: let the linker search LD_LIBRARY_PATH / system
    name = ctypes.util.find_library("ssw_aligner")
    if name:
        return ctypes.CDLL(name)

    raise OSError(
        "Cannot find libssw_aligner.so. Set SSW_ALIGNER_LIB or place the "
        "library next to this file or in <project>/build/."
    )


_lib = _find_library()

# ---------------------------------------------------------------------------
# ctypes result struct
# ---------------------------------------------------------------------------

class _ssw_result(ctypes.Structure):
    _fields_ = [
        ("score1", ctypes.c_uint32),
        ("score2", ctypes.c_uint32),
        ("ref_begin1", ctypes.c_int32),
        ("ref_end1", ctypes.c_int32),
        ("read_begin1", ctypes.c_int32),
        ("read_end1", ctypes.c_int32),
        ("ref_end2", ctypes.c_int32),
        ("cigar", ctypes.POINTER(ctypes.c_uint32)),
        ("cigarLen", ctypes.c_int32),
    ]

# ---------------------------------------------------------------------------
# Declare C function signatures
# ---------------------------------------------------------------------------

_lib.ssw_init.argtypes = [
    ctypes.POINTER(ctypes.c_int8),  # read_num
    ctypes.c_int32,                  # readLen
    ctypes.POINTER(ctypes.c_int8),  # mat
    ctypes.c_int32,                  # n
    ctypes.c_int8,                   # score_size
]
_lib.ssw_init.restype = ctypes.c_void_p

_lib.ssw_init_profile.argtypes = [
    ctypes.POINTER(ctypes.c_int8),  # read_num
    ctypes.c_int32,                  # readLen
    ctypes.POINTER(ctypes.c_int8),  # pssm
    ctypes.c_int32,                  # n
]
_lib.ssw_init_profile.restype = ctypes.c_void_p

_lib.ssw_align.argtypes = [
    ctypes.c_void_p,                 # handle
    ctypes.POINTER(ctypes.c_int8),  # ref_num
    ctypes.c_int32,                  # refLen
    ctypes.c_uint8,                  # gap_open
    ctypes.c_uint8,                  # gap_extend
    ctypes.c_uint8,                  # flag
    ctypes.c_uint16,                 # filters
    ctypes.c_int32,                  # filterd
    ctypes.c_int32,                  # maskLen
]
_lib.ssw_align.restype = _ssw_result

_lib.ssw_free_cigar.argtypes = [ctypes.POINTER(ctypes.c_uint32)]
_lib.ssw_free_cigar.restype = None

_lib.ssw_destroy.argtypes = [ctypes.c_void_p]
_lib.ssw_destroy.restype = None

# ---------------------------------------------------------------------------
# Gumbel parameter computation (ALP library)
# ---------------------------------------------------------------------------

class _gumbel_params(ctypes.Structure):
    """Mirrors the C ``gumbel_params`` struct."""
    _fields_ = [
        ("lambda_", ctypes.c_double),
        ("K", ctypes.c_double),
        ("a_I", ctypes.c_double),
        ("b_I", ctypes.c_double),
        ("a_J", ctypes.c_double),
        ("b_J", ctypes.c_double),
        ("alpha_I", ctypes.c_double),
        ("beta_I", ctypes.c_double),
        ("alpha_J", ctypes.c_double),
        ("beta_J", ctypes.c_double),
        ("sigma", ctypes.c_double),
        ("tau", ctypes.c_double),
        ("valid", ctypes.c_int),
    ]

_lib.compute_gumbel_params.argtypes = [
    ctypes.POINTER(ctypes.c_int8),   # mat
    ctypes.c_int32,                   # n
    ctypes.POINTER(ctypes.c_double),  # bg_freqs (nullable)
    ctypes.c_int32,                   # gap_open
    ctypes.c_int32,                   # gap_extend
    ctypes.c_double,                  # max_seconds
]
_lib.compute_gumbel_params.restype = _gumbel_params

# ---------------------------------------------------------------------------
# Sequence encoding tables
# ---------------------------------------------------------------------------

# 20 standard amino acids only (ARNDCQEGHILKMFPSTWYV, indices 0-19).
# Non-standard symbols: B→D(3), Z→E(6), others→A(0).
_AA_TABLE = np.array([
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  3,  4,  3,  6, 13,  7,  8,  9,  0, 11, 10, 12,  2,  0,
    14,  5,  1, 15, 16,  0, 19, 17,  0, 18,  6,  0,  0,  0,  0,  0,
     0,  0,  3,  4,  3,  6, 13,  7,  8,  9,  0, 11, 10, 12,  2,  0,
    14,  5,  1, 15, 16,  0, 19, 17,  0, 18,  6,  0,  0,  0,  0,  0,
], dtype=np.int8)

_NT_TABLE = np.array([
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
], dtype=np.int8)

_MID_TABLE = ["M", "I", "D"]

# Sequence orderings used for building flat matrices
_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
_NT_ORDER = "ACGTN"

# Robinson-Robinson amino acid background frequencies from MMseqs2's blosum62.out.
# These are extracted from the "# Background (precomputed optional):" line in
# data/blosum62.out.  Using these with use_mmseqs_aa_freqs=True reproduces
# Gumbel parameters matching MMseqs2's precomputed table.
_MMSEQS_AA_FREQS: Dict[str, float] = {
    "A": 0.07422, "R": 0.05161, "N": 0.04465, "D": 0.05363,
    "C": 0.02469, "Q": 0.03426, "E": 0.05431, "G": 0.07415,
    "H": 0.02621, "I": 0.06792, "L": 0.09891, "K": 0.05815,
    "M": 0.02499, "F": 0.04742, "P": 0.03854, "S": 0.05723,
    "T": 0.05089, "W": 0.01303, "Y": 0.03228, "V": 0.07292,
}


# ---------------------------------------------------------------------------
# AlignmentStructure
# ---------------------------------------------------------------------------

class AlignmentStructure:
    """Wraps the result of an alignment so that it is accessible to Python.

    Supports both dict-style access (``result['score1']``) **and** attribute
    access (``result.score1``).
    """

    __slots__ = (
        "_score1", "_score2",
        "_ref_begin1", "_ref_end1",
        "_read_begin1", "_read_end1", "_ref_end2",
        "_cigar_ops",
        "_read_sequence", "_reference_sequence",
        "_index_starts_at",
        "_cigar_string",
    )

    def __init__(
        self,
        result: _ssw_result,
        read_sequence: str,
        reference_sequence: str,
        index_starts_at: int,
    ):
        self._score1 = int(result.score1)
        self._score2 = int(result.score2)
        self._ref_begin1 = int(result.ref_begin1)
        self._ref_end1 = int(result.ref_end1)
        self._read_begin1 = int(result.read_begin1)
        self._read_end1 = int(result.read_end1)
        self._ref_end2 = int(result.ref_end2)

        # Decode and copy cigar so the C memory can be freed immediately.
        cigar_ops: list[tuple[int, str]] = []
        cigar_parts: list[str] = []
        if result.cigar and result.cigarLen > 0:
            for i in range(result.cigarLen):
                c = result.cigar[i]
                length = c >> 4
                op = _MID_TABLE[c & 0xF]
                cigar_ops.append((length, op))
                cigar_parts.append(f"{length}{op}")
        self._cigar_ops = cigar_ops
        self._cigar_string = "".join(cigar_parts)

        self._read_sequence = read_sequence
        self._reference_sequence = reference_sequence
        self._index_starts_at = index_starts_at

    # -- dict-style access --------------------------------------------------

    def __getitem__(self, key: str):
        return getattr(self, key)

    # -- repr / str ---------------------------------------------------------

    def __repr__(self) -> str:
        keys = [
            "optimal_alignment_score", "suboptimal_alignment_score",
            "query_begin", "query_end",
            "target_begin", "target_end_optimal", "target_end_suboptimal",
            "cigar", "query_sequence", "target_sequence",
        ]
        return "{\n%s\n}" % ",\n".join(
            "    {!r}: {!r}".format(k, self[k]) for k in keys
        )

    def __str__(self) -> str:
        score = "Score: %d" % self.optimal_alignment_score
        if self.query_sequence and self.cigar:
            target = self.aligned_target_sequence
            query = self.aligned_query_sequence
            align_len = len(query) if query else 0
            if align_len > 13:
                target = target[:10] + "..."
                query = query[:10] + "..."
            length = "Length: %d" % align_len
            return "\n".join([query, target, score, length])
        return score

    # -- score properties ---------------------------------------------------

    @property
    def optimal_alignment_score(self) -> int:
        return self._score1

    @property
    def suboptimal_alignment_score(self) -> int:
        return self._score2

    # -- position properties (with index offset) ----------------------------

    @property
    def target_begin(self) -> int:
        if self._ref_begin1 >= 0:
            return self._ref_begin1 + self._index_starts_at
        return -1

    @property
    def target_end_optimal(self) -> int:
        return self._ref_end1 + self._index_starts_at

    @property
    def target_end_suboptimal(self) -> int:
        return self._ref_end2 + self._index_starts_at

    @property
    def query_begin(self) -> int:
        if self._read_begin1 >= 0:
            return self._read_begin1 + self._index_starts_at
        return -1

    @property
    def query_end(self) -> int:
        return self._read_end1 + self._index_starts_at

    # -- cigar --------------------------------------------------------------

    @property
    def cigar(self) -> str:
        return self._cigar_string

    # -- sequences ----------------------------------------------------------

    @property
    def query_sequence(self) -> str:
        return self._read_sequence

    @property
    def target_sequence(self) -> str:
        return self._reference_sequence

    # -- aligned sequences --------------------------------------------------

    @property
    def aligned_query_sequence(self) -> Optional[str]:
        if self._read_sequence:
            return self._get_aligned_sequence(
                self._read_sequence, self._cigar_ops,
                self._read_begin1, self._read_end1, "D",
            )
        return None

    @property
    def aligned_target_sequence(self) -> Optional[str]:
        if self._reference_sequence:
            return self._get_aligned_sequence(
                self._reference_sequence, self._cigar_ops,
                self._ref_begin1, self._ref_end1, "I",
            )
        return None

    # -- zero_based helpers -------------------------------------------------

    def set_zero_based(self, is_zero_based: bool) -> None:
        self._index_starts_at = 0 if is_zero_based else 1

    def is_zero_based(self) -> bool:
        return self._index_starts_at == 0

    # -- internal -----------------------------------------------------------

    @staticmethod
    def _get_aligned_sequence(
        sequence: str,
        cigar_ops: list[tuple[int, str]],
        begin: int,
        end: int,
        gap_type: str,
    ) -> str:
        # Always use 0-based positions for slicing.
        seq = sequence[begin : end + 1]
        parts: list[str] = []
        idx = 0
        for length, op in cigar_ops:
            if op == gap_type:
                parts.append("-" * length)
            else:
                parts.append(seq[idx : idx + length])
                idx += length
        # Remainder beyond the cigar
        parts.append(seq[idx : end - begin + 1])
        return "".join(parts)


# ---------------------------------------------------------------------------
# StripedSmithWaterman
# ---------------------------------------------------------------------------

class StripedSmithWaterman:
    """Perform striped Smith-Waterman alignment.

    Instantiate with a query sequence, then call the resulting object on
    one or more target sequences.

    Parameters
    ----------
    query_sequence : str
    gap_open_penalty : int
    gap_extend_penalty : int
    score_size : int
        Ignored (kept for API compatibility).
    mask_length : int
    mask_auto : bool
    score_only : bool
    score_filter : int or None
    distance_filter : int or None
    override_skip_babp : bool
    protein : bool
    match_score : int
    mismatch_score : int
    substitution_matrix : dict[str, dict[str, int]] or None
    suppress_sequences : bool
    zero_index : bool
    """

    def __init__(
        self,
        query_sequence: str,
        *,
        gap_open_penalty: int = 5,
        gap_extend_penalty: int = 2,
        score_size: int = 2,
        mask_length: int = 15,
        mask_auto: bool = True,
        score_only: bool = False,
        score_filter: Optional[int] = None,
        distance_filter: Optional[int] = None,
        override_skip_babp: bool = False,
        protein: bool = False,
        match_score: int = 2,
        mismatch_score: int = -3,
        substitution_matrix: Optional[Dict[str, Dict[str, int]]] = None,
        suppress_sequences: bool = False,
        zero_index: bool = True,
    ):
        self._read_sequence = query_sequence

        if gap_open_penalty <= 0:
            raise ValueError("`gap_open_penalty` must be > 0")
        if gap_extend_penalty <= 0:
            raise ValueError("`gap_extend_penalty` must be > 0")

        self._gap_open = gap_open_penalty
        self._gap_extend = gap_extend_penalty
        self._distance_filter = distance_filter or 0
        self._score_filter = score_filter or 0
        self._suppress_sequences = suppress_sequences
        self._is_protein = protein
        self._index_starts_at = 0 if zero_index else 1

        # Compute bit-flag (same logic as original Cython wrapper)
        self._bit_flag = self._compute_bit_flag(override_skip_babp, score_only)

        # Mask length
        if mask_auto:
            self._mask_length = max(len(query_sequence) // 2, mask_length)
        else:
            self._mask_length = mask_length

        # Build substitution matrix (flat int8 array, row-major)
        if substitution_matrix is None:
            if protein:
                raise ValueError(
                    "Must provide a substitution matrix for protein sequences"
                )
            matrix = self._build_match_matrix(match_score, mismatch_score)
        else:
            matrix = self._convert_dict2d_to_matrix(substitution_matrix)

        # Encode query
        read_seq = self._encode_sequence(query_sequence)
        m_width = len(_AA_ORDER) if self._is_protein else len(_NT_ORDER)

        # Keep numpy arrays alive for ctypes
        self._matrix_buf = matrix
        self._read_buf = read_seq

        # Create the C handle
        self._handle = _lib.ssw_init(
            read_seq.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)),
            ctypes.c_int32(len(query_sequence)),
            matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)),
            ctypes.c_int32(m_width),
            ctypes.c_int8(score_size),
        )

    def __del__(self):
        handle = getattr(self, "_handle", None)
        if handle:
            _lib.ssw_destroy(handle)
            self._handle = None

    def __call__(self, target_sequence: str) -> AlignmentStructure:
        """Align *target_sequence* to the stored query."""
        ref_seq = self._encode_sequence(target_sequence)

        result = _lib.ssw_align(
            self._handle,
            ref_seq.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)),
            ctypes.c_int32(len(target_sequence)),
            ctypes.c_uint8(self._gap_open),
            ctypes.c_uint8(self._gap_extend),
            ctypes.c_uint8(self._bit_flag),
            ctypes.c_uint16(self._score_filter),
            ctypes.c_int32(self._distance_filter),
            ctypes.c_int32(self._mask_length),
        )

        read_seq = "" if self._suppress_sequences else self._read_sequence
        ref_str = "" if self._suppress_sequences else target_sequence

        alignment = AlignmentStructure(
            result, read_seq, ref_str, self._index_starts_at,
        )

        # Free the C-allocated cigar array immediately (data already copied)
        if result.cigar:
            _lib.ssw_free_cigar(result.cigar)

        return alignment

    # -- helpers ------------------------------------------------------------

    @staticmethod
    def _compute_bit_flag(override_skip_babp: bool, score_only: bool) -> int:
        if score_only:
            return 0
        bit_flag = 0
        if override_skip_babp:
            bit_flag |= 0x08
        # Note: we always produce full results (start pos + cigar) when not
        # score_only. The old filter modes (bit 1/2) are not used here.
        if bit_flag == 0 or bit_flag == 8:
            bit_flag |= 0x01
        return bit_flag

    def _encode_sequence(self, sequence: str) -> np.ndarray:
        table = _AA_TABLE if self._is_protein else _NT_TABLE
        arr = np.empty(len(sequence), dtype=np.int8)
        for i, ch in enumerate(sequence):
            arr[i] = table[ord(ch)]
        return arr

    def _build_match_matrix(
        self, match_score: int, mismatch_score: int
    ) -> np.ndarray:
        order = _NT_ORDER
        dict2d: Dict[str, Dict[str, int]] = {}
        for r in order:
            dict2d[r] = {}
            for c in order:
                if r == "N" or c == "N":
                    dict2d[r][c] = 0
                else:
                    dict2d[r][c] = match_score if r == c else mismatch_score
        return self._convert_dict2d_to_matrix(dict2d)

    def _convert_dict2d_to_matrix(
        self, dict2d: Dict[str, Dict[str, int]]
    ) -> np.ndarray:
        order = _AA_ORDER if self._is_protein else _NT_ORDER
        n = len(order)
        flat = np.empty(n * n, dtype=np.int8)
        i = 0
        for r in order:
            for c in order:
                flat[i] = dict2d[r][c]
                i += 1
        return flat


# ---------------------------------------------------------------------------
# SmithWatermanProfileAligner (PSSM-based profile alignment)
# ---------------------------------------------------------------------------

class SmithWatermanProfileAligner:
    """Align targets against a query using a position-specific scoring matrix.

    Unlike :class:`StripedSmithWaterman` which uses a single substitution
    matrix for all positions, this class accepts a per-position scoring
    matrix (PSSM) that assigns different scores at each query position.

    Parameters
    ----------
    query_sequence : str
        The query amino-acid sequence.  Used for consensus/traceback.
    pssm : numpy.ndarray
        Position-specific scoring matrix of shape ``(20, query_length)``
        with dtype ``int8``.  Rows correspond to amino acids in the order
        ``ARNDCQEGHILKMFPSTWYV`` (indices 0–19); columns correspond to
        query positions.  Element ``pssm[aa_idx, pos]`` is the score for
        amino acid ``aa_idx`` at query position ``pos``.
    gap_open_penalty : int
        Gap opening penalty (positive, default 11).
    gap_extend_penalty : int
        Gap extension penalty (positive, default 1).
    score_only : bool
        If True, skip traceback (faster but no CIGAR / aligned sequences).
    mask_length : int
        Minimum distance between primary and sub-optimal alignment ends.
    mask_auto : bool
        Automatically set mask_length to ``max(query_length // 2, 15)``.
    zero_index : bool
        If True (default), positions are 0-based.

    Examples
    --------
    >>> import numpy as np
    >>> from ssw_aligner import SmithWatermanProfileAligner
    >>> query = "EVQLVES"
    >>> pssm = np.zeros((20, len(query)), dtype=np.int8)
    >>> # Fill PSSM with position-specific scores ...
    >>> aligner = SmithWatermanProfileAligner(query, pssm)
    >>> result = aligner("EVQLVES")
    >>> result.optimal_alignment_score
    """

    def __init__(
        self,
        query_sequence: str,
        pssm: np.ndarray,
        *,
        gap_open_penalty: int = 11,
        gap_extend_penalty: int = 1,
        score_only: bool = False,
        mask_length: int = 15,
        mask_auto: bool = True,
        zero_index: bool = True,
    ):
        self._query_sequence = query_sequence
        query_len = len(query_sequence)
        n = len(_AA_ORDER)  # always 20

        # Validate PSSM shape
        pssm = np.asarray(pssm, dtype=np.int8)
        if pssm.shape != (n, query_len):
            raise ValueError(
                f"pssm must have shape ({n}, {query_len}), "
                f"got {pssm.shape}"
            )

        # Make sure it's C-contiguous
        if not pssm.flags["C_CONTIGUOUS"]:
            pssm = np.ascontiguousarray(pssm)
        self._pssm = pssm

        self._gap_open = gap_open_penalty
        self._gap_extend = gap_extend_penalty
        self._score_only = score_only
        self._index_starts_at = 0 if zero_index else 1

        if mask_auto:
            self._mask_length = max(query_len // 2, 15)
        else:
            self._mask_length = mask_length

        # Encode query to numeric
        read_num = np.array(
            [_AA_TABLE[ord(ch)] for ch in query_sequence], dtype=np.int8
        )

        self._handle = _lib.ssw_init_profile(
            read_num.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)),
            ctypes.c_int32(query_len),
            pssm.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)),
            ctypes.c_int32(n),
        )

        # prevent GC of arrays while handle is alive
        self._read_num = read_num

    def __del__(self):
        if hasattr(self, "_handle") and self._handle:
            _lib.ssw_destroy(self._handle)
            self._handle = None

    def __call__(self, target_sequence: str) -> AlignmentStructure:
        """Align *target_sequence* against the stored query profile.

        Parameters
        ----------
        target_sequence : str
            Target (reference) amino-acid sequence.

        Returns
        -------
        AlignmentStructure
        """
        ref_num = np.array(
            [_AA_TABLE[ord(ch)] for ch in target_sequence], dtype=np.int8
        )

        flag = 0 if self._score_only else 2

        raw = _lib.ssw_align(
            self._handle,
            ref_num.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)),
            ctypes.c_int32(len(target_sequence)),
            ctypes.c_uint8(self._gap_open),
            ctypes.c_uint8(self._gap_extend),
            ctypes.c_uint8(flag),
            ctypes.c_uint16(0),
            ctypes.c_int32(0),
            ctypes.c_int32(self._mask_length),
        )

        result = AlignmentStructure(
            raw,
            self._query_sequence,
            target_sequence,
            self._index_starts_at,
        )

        # Free C-allocated cigar
        if raw.cigar:
            _lib.ssw_free_cigar(raw.cigar)

        return result


# ---------------------------------------------------------------------------
# Gumbel statistical parameters
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Precomputed Gumbel parameters (from MMseqs2's EvalueComputation.h)
# ---------------------------------------------------------------------------
# Each entry maps (matrix_name, gap_open, gap_extend) → 12-tuple of Gumbel
# parameters in order: (lambda, K, a_I, b_I, a_J, b_J, alpha_I, beta_I,
# alpha_J, beta_J, sigma, tau).  These were estimated by the ALP Monte Carlo
# method with seed=42 and then hard-coded in MMseqs2 for the most common
# matrix/gap combinations.

_PRECOMPUTED_PARAMS: Dict[tuple, tuple] = {
    # nucleotide.out, gap_open=7, gap_extend=1
    ("nucleotide", 7, 1): (
        1.0960171987681839,    0.33538787507026158,
        2.0290734315292083,   -0.46514786408422282,
        2.0290734315292083,   -0.46514786408422282,
        5.0543294182155085,   15.130999712620039,
        5.0543294182155085,   15.130999712620039,
        5.0543962679167036,   15.129930117400917,
    ),
    # nucleotide.out, gap_open=5, gap_extend=2
    ("nucleotide", 5, 2): (
        0.62092274139392822363, 0.35177597988201619872,
        0.74528059208662511548, -0.71027220445456995535,
        0.74528059208662511548, -0.71027220445456995535,
        1.0135243407674570104, -2.5226486486783059604,
        1.0135243407674570104, -2.5226486486783059604,
        1.0031949332622873694, -2.3780369436059309862,
    ),
    # blosum62.out, gap_open=11, gap_extend=1
    ("blosum62", 11, 1): (
        0.27359865037097330642, 0.044620920658722244834,
        1.5938724404943873658, -19.959867650284412122,
        1.5938724404943873658, -19.959867650284412122,
        30.455610143099914211, -622.28684628915891608,
        30.455610143099914211, -622.28684628915891608,
        29.602444874818868215, -601.81087985041381216,
    ),
}


def _lookup_precomputed(
    matrix_name: Optional[str], gap_open: int, gap_extend: int
) -> Optional["GumbelParams"]:
    """Return precomputed GumbelParams if the configuration is known."""
    if matrix_name is None:
        return None
    key = (matrix_name.lower(), gap_open, gap_extend)
    vals = _PRECOMPUTED_PARAMS.get(key)
    if vals is None:
        return None
    return GumbelParams(*vals)


class GumbelParams:
    """Gumbel distribution parameters for Smith-Waterman E-value computation.

    These parameters describe the extreme-value distribution of local
    alignment scores and are either looked up from a precomputed table
    (for well-known matrix/gap combinations) or computed via the ALP
    (Ascending Ladder Points) Monte Carlo simulation.

    Attributes
    ----------
    lambda_ : float
        Gumbel distribution scale parameter.
    K : float
        Gumbel distribution prefactor.
    a_I, b_I : float
        Length adjustment (slope, intercept) for sequence I.
    a_J, b_J : float
        Length adjustment (slope, intercept) for sequence J.
    alpha_I, beta_I : float
        Variance (slope, intercept) for sequence I.
    alpha_J, beta_J : float
        Variance (slope, intercept) for sequence J.
    sigma, tau : float
        Aggregate variance (slope, intercept).

    Notes
    -----
    E-value is computed as:

    .. math::

        E = K \\cdot e^{-\\lambda \\cdot S} \\cdot m' \\cdot n'

    where *S* is the raw alignment score and *m'*, *n'* are the
    finite-size-corrected sequence lengths.

    Bit-score is:

    .. math::

        \\text{bit\\_score} = \\frac{\\lambda \\cdot S - \\ln K}{\\ln 2}
    """

    __slots__ = (
        "lambda_", "K",
        "a_I", "b_I", "a_J", "b_J",
        "alpha_I", "beta_I", "alpha_J", "beta_J",
        "sigma", "tau",
    )

    def __init__(
        self,
        lambda_: float, K: float,
        a_I: float, b_I: float,
        a_J: float, b_J: float,
        alpha_I: float, beta_I: float,
        alpha_J: float, beta_J: float,
        sigma: float, tau: float,
    ):
        self.lambda_ = lambda_
        self.K = K
        self.a_I = a_I
        self.b_I = b_I
        self.a_J = a_J
        self.b_J = b_J
        self.alpha_I = alpha_I
        self.beta_I = beta_I
        self.alpha_J = alpha_J
        self.beta_J = beta_J
        self.sigma = sigma
        self.tau = tau

    def bit_score(self, raw_score: float) -> float:
        """Convert a raw alignment score to a bit-score."""
        import math
        return (self.lambda_ * raw_score - math.log(self.K)) / math.log(2.0)

    def evalue(
        self, raw_score: float, query_length: int, db_residues: int
    ) -> float:
        """Compute an approximate E-value (without finite-size correction).

        Parameters
        ----------
        raw_score : float
            Raw alignment score from Smith-Waterman.
        query_length : int
            Length of the query sequence.
        db_residues : int
            Total number of residues in the database (or target length
            for a single pairwise comparison).

        Returns
        -------
        float
            E-value.
        """
        import math
        m = max(query_length - self.a_I * raw_score - self.b_I, 1.0)
        n = max(db_residues - self.a_J * raw_score - self.b_J, 1.0)
        return self.K * math.exp(-self.lambda_ * raw_score) * m * n

    def __repr__(self) -> str:
        return (
            f"GumbelParams(lambda_={self.lambda_:.6g}, K={self.K:.6g}, "
            f"sigma={self.sigma:.6g}, tau={self.tau:.6g})"
        )


def compute_gumbel_params(
    substitution_matrix: Dict[str, Dict[str, int]],
    gap_open: int,
    gap_extend: int,
    *,
    protein: bool = True,
    matrix_name: Optional[str] = None,
    bg_freqs: Optional[Dict[str, float]] = None,
    use_mmseqs_aa_freqs: bool = False,
    recalculate_gumbel_params: bool = False,
    max_seconds: float = 60.0,
) -> GumbelParams:
    """Compute Gumbel statistical parameters for Smith-Waterman alignment.

    For well-known matrix/gap combinations (BLOSUM62 11/1, nucleotide 7/1
    and 5/2) **precomputed** parameters matching MMseqs2 are returned
    instantly.  Pass ``recalculate_gumbel_params=True`` to run the ALP
    Monte Carlo simulation instead.

    Parameters
    ----------
    substitution_matrix : dict of dict of int
        2D dictionary mapping ``{row_char: {col_char: score}}``.  
        For proteins, keys should cover the 20 standard amino acids
        (``ARNDCQEGHILKMFPSTWYV``); for nucleotides, ``ACGTN``.
    gap_open : int
        Gap opening penalty (positive integer, e.g. 11 for BLOSUM62).
    gap_extend : int
        Gap extension penalty (positive integer, e.g. 1 for BLOSUM62).
    protein : bool, default True
        Whether the matrix is for protein (20-letter) or nucleotide
        (5-letter) alphabets.
    matrix_name : str, optional
        A short identifier for the substitution matrix (e.g.
        ``"blosum62"`` or ``"nucleotide"``).  When provided and
        ``recalculate_gumbel_params`` is ``False``, the function first checks a table
        of precomputed Gumbel parameters for the combination
        ``(matrix_name, gap_open, gap_extend)``.

        If ``None``, the function attempts to auto-detect the name by
        inspecting the matrix scores (BLOSUM62 if protein, nucleotide
        if not protein).
    bg_freqs : dict of {str: float}, optional
        Background frequencies for each symbol. If ``None``, uniform
        frequencies (1/n) are used.  Ignored when *use_mmseqs_aa_freqs*
        is ``True``.
    use_mmseqs_aa_freqs : bool, default False
        When ``True`` (protein mode only), use the Robinson-Robinson
        amino acid background frequencies from MMseqs2's built-in
        BLOSUM62.
    recalculate_gumbel_params : bool, default False
        When ``True``, always run the ALP Monte Carlo simulation even
        if precomputed parameters are available.  Useful for
        non-standard background frequencies or custom matrices.
    max_seconds : float, default 60.0
        Maximum wall-clock time for the Monte Carlo simulation.

    Returns
    -------
    GumbelParams
        The estimated (or precomputed) Gumbel parameters.

    Raises
    ------
    RuntimeError
        If ALP fails to converge within the time limit.

    Examples
    --------
    >>> import blosum
    >>> mat = blosum.BLOSUM(62)
    >>> # Instant — uses precomputed parameters
    >>> params = compute_gumbel_params(mat, gap_open=11, gap_extend=1,
    ...                               matrix_name="blosum62")
    >>> params.lambda_  # 0.2736
    >>> # Force ALP recalculation
    >>> params2 = compute_gumbel_params(mat, gap_open=11, gap_extend=1,
    ...                                recalculate_gumbel_params=True,
    ...                                use_mmseqs_aa_freqs=True)
    """
    if use_mmseqs_aa_freqs:
        if not protein:
            raise ValueError(
                "use_mmseqs_aa_freqs=True is only valid for protein matrices"
            )
        bg_freqs = _MMSEQS_AA_FREQS

    # --- auto-detect matrix name if not provided --------------------------
    if matrix_name is None:
        matrix_name = _detect_matrix_name(substitution_matrix, protein)

    # --- try precomputed table first (unless forced to recalculate) -------
    if not recalculate_gumbel_params:
        precomputed = _lookup_precomputed(matrix_name, gap_open, gap_extend)
        if precomputed is not None:
            return precomputed

    # --- fall back to ALP Monte Carlo -------------------------------------
    order = _AA_ORDER if protein else _NT_ORDER
    n = len(order)

    # Build flat int8 matrix from dict2d
    flat = np.empty(n * n, dtype=np.int8)
    idx = 0
    for r in order:
        for c in order:
            flat[idx] = substitution_matrix[r][c]
            idx += 1

    # Build background frequency array if provided
    bg_ptr = None
    bg_arr = None
    if bg_freqs is not None:
        bg_arr = np.array([bg_freqs.get(ch, 0.0) for ch in order],
                          dtype=np.float64)
        bg_ptr = bg_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    result = _lib.compute_gumbel_params(
        flat.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)),
        ctypes.c_int32(n),
        bg_ptr,
        ctypes.c_int32(gap_open),
        ctypes.c_int32(gap_extend),
        ctypes.c_double(max_seconds),
    )

    if not result.valid:
        raise RuntimeError(
            "ALP did not converge for the given substitution matrix and gap "
            "penalties. Try increasing max_seconds or check your inputs."
        )

    return GumbelParams(
        lambda_=result.lambda_,
        K=result.K,
        a_I=result.a_I,
        b_I=result.b_I,
        a_J=result.a_J,
        b_J=result.b_J,
        alpha_I=result.alpha_I,
        beta_I=result.beta_I,
        alpha_J=result.alpha_J,
        beta_J=result.beta_J,
        sigma=result.sigma,
        tau=result.tau,
    )


def _detect_matrix_name(
    substitution_matrix: Dict[str, Dict[str, int]],
    protein: bool,
) -> Optional[str]:
    """Best-effort auto-detection of the substitution matrix name."""
    if not protein:
        # There is only one nucleotide matrix in the precomputed table.
        return "nucleotide"
    # Quick BLOSUM62 fingerprint: check a few signature scores.
    try:
        if (
            substitution_matrix["A"]["A"] == 4
            and substitution_matrix["W"]["W"] == 11
            and substitution_matrix["A"]["R"] == -1
            and substitution_matrix["D"]["E"] == 2
        ):
            return "blosum62"
    except KeyError:
        pass
    return None
