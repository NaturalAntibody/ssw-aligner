from enum import Enum


class EqualScoreSide(str, Enum):
    """
    Alignment end(s) extended over equal-score residues (StripedSmithWaterman's extend_equal_score), see
    docs/equal_score_extension.md. The string values ("start", "end", "both") are accepted as well.
    """

    START = "start"
    END = "end"
    BOTH = "both"
