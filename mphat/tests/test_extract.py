"""
Unit and regression test for the mPHAT package.
"""

# Import package, test suite, and other packages as needed
import sys

import numpy

import mphat.extract


def test_mphat_imported():
    """
    Sample test. This will always pass so long as import statements worked.

    """
    assert "mphat.extract" in sys.modules


def test_find_transitions():
    """
    Tests whether the ``mphat.extract`` ``find_transitions()`` function that calculates
    transitions from standard MD actually works as expected.

    """
    test_array = numpy.asarray([0, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 2, 0, 2, 1, 2, 1, 0])
    test_output = (numpy.array([0, 19, 21]), numpy.array([6, 8, 10, 11, 12, 13, 23, 25]), [[0, 6], [19, 23], [21, 23]])

    output_array = mphat.extract.find_transitions(test_array, 0, 1)

    assert all(numpy.array_equal(x, y) for (x, y) in zip(output_array, test_output))


def test_clean_self_to_self():
    """
    Tests whether the ``mphat.extract`` ``clean_self_to_self()` function that removes
    self to self transitions from standard MD actually works as expected.

    """
    test_array = [[0, 6], [19, 23], [21, 23]]
    test_output = [[0, 6], [21, 23]]

    output_array = mphat.extract.clean_self_to_self(test_array)

    assert numpy.array_equal(output_array, test_output)
