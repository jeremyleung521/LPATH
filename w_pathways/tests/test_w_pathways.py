"""
Unit and regression test for the w_pathways package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import w_pathways


def test_w_pathways_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "w_pathways" in sys.modules
