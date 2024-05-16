"""
Unit and regression test for the lpath package.
"""
# Import package, test suite, and other packages as needed
import sys
import pytest

import argparse
from argparse import ArgumentTypeError

import lpath.argparser
from lpath.argparser import InvalidArgumentError


def test_lpath_argparser_imported():
    """
    Sample test. This will always pass so long as import statements worked.

    """
    assert "lpath.argparser" in sys.modules


class TestCheckNonNeg:
    """
    Test to see if check_non_neg() is working properly.

    """

    def test_fail(self):
        with pytest.raises(InvalidArgumentError):
            lpath.argparser.check_non_neg(-1)

        with pytest.raises(ArgumentTypeError):
            lpath.argparser.check_non_neg('abc')

    def test_success(self):
        assert lpath.argparser.check_non_neg(3) == 3


class TestCheckNonNegFloat:
    """
    Test to see if check_non_neg_float() is working properly.

    """

    def test_fail(self):
        with pytest.raises(InvalidArgumentError):
            lpath.argparser.check_non_neg_float(-1)

        with pytest.raises(ArgumentTypeError):
            lpath.argparser.check_non_neg_float('abc')

    def test_success(self):
        return_val = lpath.argparser.check_non_neg_float(3)
        assert return_val == 3
        assert isinstance(return_val, float)


class TestCheckPositive:
    """
    Test to see if check_positive() is working properly.

    """

    def test_fail(self):
        with pytest.raises(InvalidArgumentError):
            lpath.argparser.check_positive(0)

        with pytest.raises(ArgumentTypeError):
            lpath.argparser.check_positive('abc')

    def test_success(self):
        assert lpath.argparser.check_positive(3) == 3


class TestCheckLessThree:
    """
    Test to see if check_less_three works properly

    """

    def test_fail(self):
        """
        Test to see if check_less_three() is working properly.

        """

        with pytest.raises(InvalidArgumentError):
            lpath.argparser.check_less_three(-1)

        with pytest.raises(ArgumentTypeError):
            lpath.argparser.check_positive('abc')

    def test_success(self):
        """
        Test to see if check_less_three() is working properly.

        """

        assert lpath.argparser.check_less_three(2) == 2


class TestParsers:
    """
    Classes of tests to see if all arguments are added correctly.

    """

    def test_create_parser(self, create_ref_parser):
        """
        Test to see if a parser is created correctly.

        """
        output = lpath.argparser.create_parser()

        assert isinstance(output, argparse.ArgumentParser)

    def test_dimensions(self, create_ref_parser):
        """
        See if returned objects are legit.

        """
        output, _ = create_ref_parser

        assert isinstance(output, argparse.ArgumentParser)
        # Need to change as number of options increase
        assert len(output._actions) == 67

    def test_common_arguments(self, create_ref_parser):
        """
        Test to see if common args are in the parser.

        """
        output, test_output = create_ref_parser
        test_input = ['-h', '-od', '-st', '-s', '--debug', '-we', '-W', '-A', '-r', '--tui']

        for option in test_input:
            assert option in test_output

    def test_discretize_arguments(self, create_ref_parser):
        """
        Test to see if discretize args are in the parser.

        """
        output, test_output = create_ref_parser
        test_input = ['-i', '-o', '-af', '-ar']

        for option in test_input:
            assert option in test_output

    def test_extract_arguments(self, create_ref_parser):
        """
        Test to see if extract args are in the parser.

        """
        output, test_output = create_ref_parser
        test_input = ['-ei', '-eo', '-ss', '-ts', '-pc', '-ef', '-fs', '-tb', '-el',
                      '-R', '-NR', '-t', '--first', '--last', '--hdf5', '-a', '-aa',
                      '-rw', '-oj', '-oe', '-se', '-ot']

        for option in test_input:
            assert option in test_output

    def test_match_arguments(self, create_ref_parser):
        """
        Test to see if match args are in the parser.

        """
        output, test_output = create_ref_parser
        test_input = ['-ip', '-op', '-co', '-me', '-ra', '-seq', '-str', '-mm', '-mr',
                      '-re', '-cc', '-dR', '-dN', '-nd', '-dF', '-dP', '-c', '-ex', '-fp']

        for option in test_input:
            assert option in test_output

    def test_plot_arguments(self, create_ref_parser):
        """
        Test to see if plot args are in the parser.

        """
        output, test_output = create_ref_parser
        test_input = ['-ipl', '-icl', '-pdF', '-pod', '-sty', '-mpl', '-col', '-pdt', '-pds',
                      '-pdh', '-nc', '-prl', '-pto']

        for option in test_input:
            assert option in test_output
