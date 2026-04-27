"""The single source of randomness for the package."""
from __future__ import annotations

import random
import numpy as np


def make_rng(seed: int) -> np.random.Generator:
    """Return a new numpy Generator seeded from `seed`."""
    return np.random.default_rng(seed)


def seeded_python_random(seed: int) -> random.Random:
    """Return a stdlib Random instance seeded from `seed`. Use only when an
    API requires the stdlib `random` interface (e.g. random.choices)."""
    return random.Random(seed)
