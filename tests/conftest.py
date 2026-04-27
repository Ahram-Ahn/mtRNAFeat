"""Common pytest fixtures and skip markers."""
from __future__ import annotations

from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
TEST_DATA = PROJECT / "test_data"


@pytest.fixture(scope="session")
def mini_config_path() -> Path:
    return TEST_DATA / "mini.config.yaml"


@pytest.fixture(scope="session")
def mini_human_db() -> Path:
    return TEST_DATA / "mini_human.db"


@pytest.fixture(scope="session")
def mini_yeast_db() -> Path:
    return TEST_DATA / "mini_yeast.db"


def has_rna() -> bool:
    try:
        import RNA  # noqa: F401
        return True
    except ImportError:
        return False


needs_rna = pytest.mark.skipif(not has_rna(), reason="ViennaRNA not installed")
