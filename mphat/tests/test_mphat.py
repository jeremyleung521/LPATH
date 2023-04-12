"""
Unit and regression test for the mPHAT package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
import mphat

def test_mphat_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "mphat" in sys.modules
