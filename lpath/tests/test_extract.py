"""
Unit and regression test for the lpath package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
import numpy
import lpath.extract


def test_lpath_extract_imported():
    """
    Sample test. This will always pass so long as import statements worked.

    """
    assert "lpath.extract" in sys.modules


def test_find_transitions():
    """
    Tests whether the ``lpath.extract`` ``find_transitions()`` function that calculates
    transitions from standard MD actually works as expected.

    """
    test_array = numpy.asarray([0, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 2, 0, 2, 1, 2, 1, 0])
    test_output = (numpy.array([0, 19, 21, 26]), numpy.array([6, 8, 10, 11, 12, 13, 23, 25]),
                   [[0, 6], [19, 23], [21, 23]])

    output_array = lpath.extract.find_transitions(test_array, 0, 1)

    assert all(numpy.array_equal(x, y) for (x, y) in zip(output_array, test_output))


def test_clean_self_to_self():
    """
    Tests whether the ``lpath.extract`` ``clean_self_to_self()` function that removes
    self to self transitions from standard MD actually works as expected.

    """
    test_array = [[0, 6], [19, 23], [21, 23]]
    test_output = [[0, 6], [21, 23]]

    output_array = lpath.extract.clean_self_to_self(test_array)

    assert numpy.array_equal(output_array, test_output)


def test_count_tmatrix_row():
    """
    Tests whether the ``lpath.extract`` ``count_tmatrix_row()` function that calculates
    weights for standard MD successful trajectories actually works as expected.

    """
    test_array = numpy.asarray([0, 3, 3, 1, 2, 2, 1, 2, 1, 2, 0, 0, 2, 3, 1, 1, 2, 2, 2, 2, 0, 2, 2, 0, 3, 1, 2, 1, 0])
    test_source_index = numpy.array(numpy.where(test_array == 0))[0]
    test_output = 0.6

    output_array = lpath.extract.count_tmatrix_row(test_source_index, test_array, 3, 0, 1)

    assert output_array == test_output


def test_assign_color_frame():
    """
    Tests whether it successfully assigns color (i.e. 0 or first time it arrives at source after target) to target_frame

    """
    test_array = numpy.asarray([0, 3, 3, 1, 2, 2, 1, 2, 1, 2, 0, 0, 2, 3, 1, 1, 2, 2, 2, 2, 0, 2, 2, 0, 3, 1, 2, 1, 0])
    test_source_index = numpy.array(numpy.where(test_array == 0))[0]
    test_sink_index = numpy.array(numpy.where(test_array == 2))[0]
    test_output = {4: 0, 5: 0, 7: 0, 9: 0, 12: 10, 16: 10, 17: 10, 18: 10, 19: 10, 21: 20, 22: 20, 26: 23}

    output_array = lpath.extract.assign_color_frame(test_source_index, test_sink_index)

    assert output_array == test_output


def test_raise_warnings():
    """
    Tests whether it successfully raise goes through warnings.

    """
    # Catching case where no successful trajectories
    with pytest.raises(lpath.io.EmptyOutputError):
        lpath.extract.raise_warnings([],False)

    lpath.extract.raise_warnings([[]], True)


def test_create_pickle_obj():
    """
    Test to make sure create_pickle_obj is working.

    """
    output = lpath.extract.create_pickle_obj([], [], [], features=None)

    assert output == []

