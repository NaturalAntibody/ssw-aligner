from collections.abc import Iterator
from itertools import groupby


def get_cigar_op_groups(cigar: str) -> Iterator[tuple[int, str]]:
    cigar_groups = (list(grouper) for _, grouper in groupby(cigar, lambda character: character.isdigit()))
    for group_size_grouper, op_grouper in zip(cigar_groups, cigar_groups):
        yield int("".join(group_size_grouper)), op_grouper[0]


def unfold_cigar(cigar: str) -> str:
    result = []
    for group_size, operation in get_cigar_op_groups(cigar):
        for _ in range(group_size):
            result.append(operation)
    return "".join(result)


def fold_cigar(alignment_str: str) -> str:
    return "".join(f"{sum(1 for _ in op_group)}{op}" for op, op_group in groupby(alignment_str))