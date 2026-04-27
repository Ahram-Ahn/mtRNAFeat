"""Project a global dot-bracket structure onto a window or prefix, dropping
pairs whose partner falls outside. Also sanitize_dangling, used when the
DMS structure is truncated to a 5' prefix."""
from __future__ import annotations

from mtrnafeat.core.structure import pair_table


def sanitize_dangling(structure: str) -> str:
    """Replace any unmatched bracket with '.' (left- and right-dangling)."""
    chars = list(structure)
    stack: list[int] = []
    for i, ch in enumerate(chars):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                stack.pop()
            else:
                chars[i] = "."
    for i in stack:
        chars[i] = "."
    return "".join(chars)


def project_structure_to_window(structure: str, start_0based: int, end_exclusive: int) -> str:
    """Keep only pairs whose two endpoints both fall in [start, end).

    Returns a new dot-bracket string of length (end - start), reindexed to the
    local window. Pairs that span out of the window are replaced by '.'.
    """
    pt = pair_table(structure)
    n = end_exclusive - start_0based
    local = ["."] * n
    for ig in range(start_0based, end_exclusive):
        jg = pt.get(ig)
        if jg is None:
            continue
        if start_0based <= jg < end_exclusive:
            i_local = ig - start_0based
            j_local = jg - start_0based
            if i_local < j_local:
                local[i_local] = "("
                local[j_local] = ")"
    return "".join(local)


def truncate_prefix(structure: str, length: int) -> str:
    """Truncate structure to the first `length` chars and sanitize dangling."""
    return sanitize_dangling(structure[:length])
