"""Dot-bracket ↔ pair-table conversions and structural-element decomposition.

The element parser is lifted from legacy/10.feature_parsing.py and refined.
"""
from __future__ import annotations

from dataclasses import dataclass, field


def pair_table(structure: str) -> dict[int, int]:
    """Return a dict mapping every paired index → its partner. 0-based."""
    stack: list[int] = []
    pairs: dict[int, int] = {}
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise ValueError("Unbalanced ')' in structure")
            j = stack.pop()
            pairs[i] = j
            pairs[j] = i
    if stack:
        raise ValueError(f"{len(stack)} unclosed '(' in structure")
    return pairs


def extract_pairs(structure: str) -> list[tuple[int, int]]:
    """Return [(i, j) ...] with i<j, sorted by i. 0-based."""
    pt = pair_table(structure)
    return sorted({(min(i, j), max(i, j)) for i, j in pt.items()}, key=lambda p: p[0])


def filter_max_bp_span(structure: str, max_bp_span: int | None) -> tuple[str, int]:
    """Replace pairs with span > max_bp_span by '.'. Returns (new_struct, n_removed)."""
    if max_bp_span is None:
        return structure, 0
    chars = list(structure)
    removed = 0
    for i, j in extract_pairs(structure):
        if (j - i) > int(max_bp_span):
            chars[i] = "."
            chars[j] = "."
            removed += 1
    return "".join(chars), removed


def paired_fraction(structure: str) -> float:
    n = len(structure)
    if n == 0:
        return 0.0
    return sum(1 for ch in structure if ch in "()") / n


@dataclass
class ElementCensus:
    macro_helix: list[int] = field(default_factory=list)  # bp counts
    hairpin: list[int] = field(default_factory=list)      # loop sizes
    bulge: list[int] = field(default_factory=list)
    internal_loop: list[int] = field(default_factory=list)
    multi_branch: list[int] = field(default_factory=list)
    spans: list[int] = field(default_factory=list)         # j-i for every pair


def parse_element_sizes(structure: str, max_loop_artifact_size: int = 50) -> ElementCensus:
    """Decompose dot-bracket into structural elements with the legacy logic."""
    pairs = extract_pairs(structure)
    if not pairs:
        return ElementCensus()
    pair_set = set(pairs)
    spans = [j - i for i, j in pairs]

    # Group pairs into stems (i+1, j-1) consecutive runs
    stems: list[list[tuple[int, int]]] = []
    visited: set[tuple[int, int]] = set()
    for p in pairs:
        if p in visited:
            continue
        stem = [p]
        visited.add(p)
        ci, cj = p
        while (ci + 1, cj - 1) in pair_set:
            ci += 1
            cj -= 1
            stem.append((ci, cj))
            visited.add((ci, cj))
        stems.append(stem)
    stems.sort(key=lambda s: s[0][0])

    children_of: dict[int, list[int]] = {i: [] for i in range(len(stems))}
    for i in range(len(stems)):
        in_i, in_j = stems[i][-1]
        candidates = []
        for j in range(len(stems)):
            if i == j:
                continue
            o_i, o_j = stems[j][0]
            if in_i < o_i and o_j < in_j:
                candidates.append((j, o_i, o_j))
        candidates.sort(key=lambda x: x[1])
        max_j = -1
        for idx, c_i, c_j in candidates:
            if c_i > max_j:
                children_of[i].append(idx)
                max_j = c_j

    macro_helices: list[list[int]] = []
    hairpins: list[int] = []
    bulges: list[int] = []
    internals: list[int] = []
    multiloops: list[int] = []
    stem_to_macro: dict[int, int] = {}

    for i in range(len(stems)):
        if i in stem_to_macro:
            continue
        current_macro = [i]
        stem_to_macro[i] = len(macro_helices)
        cur = i
        while len(children_of[cur]) == 1:
            child = children_of[cur][0]
            in_i, in_j = stems[cur][-1]
            o_i, o_j = stems[child][0]
            gap_left = o_i - in_i - 1
            gap_right = in_j - o_j - 1
            gap_total = gap_left + gap_right

            if 0 < gap_total <= max_loop_artifact_size:
                if gap_left == 0 or gap_right == 0:
                    bulges.append(gap_total)
                else:
                    internals.append(gap_total)

            if gap_total <= 4:
                current_macro.append(child)
                stem_to_macro[child] = len(macro_helices)
                cur = child
            else:
                break

        if len(children_of[cur]) == 0:
            in_i, in_j = stems[cur][-1]
            hp = in_j - in_i - 1
            if hp <= max_loop_artifact_size:
                hairpins.append(hp)
        elif len(children_of[cur]) > 1:
            in_i, in_j = stems[cur][-1]
            child_span = sum(stems[c][0][1] - stems[c][0][0] + 1 for c in children_of[cur])
            mb = (in_j - in_i - 1) - child_span
            if mb <= max_loop_artifact_size:
                multiloops.append(mb)

        macro_helices.append(current_macro)

    macro_sizes = [sum(len(stems[idx]) for idx in m) for m in macro_helices]
    return ElementCensus(
        macro_helix=macro_sizes,
        hairpin=hairpins,
        bulge=bulges,
        internal_loop=internals,
        multi_branch=multiloops,
        spans=spans,
    )
