"""Altschul-Erikson dinucleotide-preserving shuffle.

This is the same algorithm as `uShuffle` (Jiang et al. 2008,
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192).

Used by mtrnafeat.analysis.significance to build a per-window null distribution
that respects the dinucleotide composition of the sequence.
"""
from __future__ import annotations

import numpy as np

# Altschul-Erikson Euler-tour shuffle of dinucleotides.

def dinuc_shuffle(seq: str, rng: np.random.Generator) -> str:
    """Return one dinucleotide-preserving shuffle of `seq`.

    Implements the Altschul-Erikson algorithm: build a multigraph of
    dinucleotide transitions, generate a random Eulerian path that starts
    with seq[0] and ends with seq[-1], read off the resulting sequence.
    """
    seq = seq.upper()
    n = len(seq)
    if n < 3:
        return seq

    alphabet = sorted(set(seq))
    last_char = seq[-1]
    # Successor multisets keyed by current char.
    succ: dict[str, list[str]] = {a: [] for a in alphabet}
    for i in range(n - 1):
        succ[seq[i]].append(seq[i + 1])

    # Step 1: Build a random spanning tree rooted at last_char by reversing
    # the algorithm: for each node != last_char, pick a random outgoing edge
    # and mark it as a tree edge.
    while True:
        last_edge: dict[str, str] = {}
        ok = True
        for a in alphabet:
            if a == last_char:
                continue
            if not succ[a]:
                ok = False
                break
            last_edge[a] = succ[a][int(rng.integers(0, len(succ[a])))]
        if not ok:
            return seq
        # Verify reachability of last_char from every node via last_edge tree.
        reachable = {last_char}
        changed = True
        while changed:
            changed = False
            for a in alphabet:
                if a in reachable:
                    continue
                if a in last_edge and last_edge[a] in reachable:
                    reachable.add(a)
                    changed = True
        if all(a in reachable for a in alphabet):
            break

    # Step 2: For each node, shuffle the *non-tree* outgoing edges and append
    # the tree edge (if any) to the END.
    edges: dict[str, list[str]] = {a: [] for a in alphabet}
    for a in alphabet:
        bucket = list(succ[a])
        if a in last_edge:
            tree_target = last_edge[a]
            # Remove exactly one occurrence of tree_target from bucket.
            bucket.remove(tree_target)
            rng.shuffle(bucket)
            bucket.append(tree_target)
        else:
            rng.shuffle(bucket)
        edges[a] = bucket

    # Step 3: Walk from seq[0], consuming edges in order.
    out = [seq[0]]
    cur = seq[0]
    while edges[cur]:
        nxt = edges[cur].pop(0)
        out.append(nxt)
        cur = nxt
    return "".join(out)


def random_gc_sequence(length: int, gc_fraction: float, rng: np.random.Generator) -> str:
    """Generate a random RNA sequence at given GC fraction (G:C symmetric, A:U symmetric)."""
    p_g = p_c = gc_fraction / 2
    p_a = p_u = (1.0 - gc_fraction) / 2
    return random_sequence_with_freqs(length, {"A": p_a, "U": p_u, "G": p_g, "C": p_c}, rng)


def random_sequence_with_freqs(length: int, freqs: dict[str, float],
                                rng: np.random.Generator) -> str:
    """Generate a random RNA sequence sampled from arbitrary nucleotide
    frequencies. `freqs` is a dict over {A, U, G, C}; values are normalized
    to sum to 1.

    Use this instead of `random_gc_sequence` when you care about G/C asymmetry
    (e.g. human mt-mRNA: C ≫ G on H-strand transcripts).
    """
    chars = np.array(list("AUGC"))
    weights = np.array([float(freqs.get(b, 0.0)) for b in "AUGC"], dtype=float)
    s = weights.sum()
    if s <= 0:
        raise ValueError("nucleotide frequencies must be positive and sum > 0")
    weights = weights / s
    idx = rng.choice(4, size=length, p=weights)
    return "".join(chars[idx].tolist())
