"""
Unit and regression test for the lpath package.
"""
import argparse
# Import package, test suite, and other packages as needed
import sys
import lpath.argparser
from lpath.argparser import InvalidArgumentError

from argparse import ArgumentTypeError


def test_lpath_argparser_imported():
    """
    Sample test. This will always pass so long as import statements worked.

    """
    assert "lpath.argparser" in sys.modules


def test_check_non_neg():
    """
    Test to see if check_non_neg() is working properly.

    """
    lpath.argparser.check_non_neg(3)

    try:
        lpath.argparser.check_non_neg(-1)
    except InvalidArgumentError:
        pass

    try:
        lpath.argparser.check_non_neg('abc')
    except ArgumentTypeError:
        pass


def test_check_positive():
    """
    Test to see if check_positive() is working properly.

    """
    lpath.argparser.check_positive(3)

    try:
        lpath.argparser.check_positive(0)
    except InvalidArgumentError:
        pass

    try:
        lpath.argparser.check_positive('abc')
    except ArgumentTypeError:
        pass


def test_check_less_three():
    """
    Test to see if check_less_three() is working properly.

    """
    lpath.argparser.check_less_three(2)

    try:
        lpath.argparser.check_less_three(-1)
    except InvalidArgumentError:
        pass

    try:
        lpath.argparser.check_positive('abc')
    except ArgumentTypeError:
        pass


def test_create_parser():
    """
    Test to see if a parser is created correctly.

    """
    output = lpath.argparser.create_parser()

    assert isinstance(output, argparse.ArgumentParser)


def test_add_common_args():
    """
    Test to see if arguments are added correctly correctly.

    """
    output = lpath.argparser.add_common_args()

    assert isinstance(output, argparse.ArgumentParser)
    assert len(output._actions) == 9

    test_output = [action for actions in output._actions for action in actions.option_strings]
    test_input = ['-h', '-od', '-st', '-s', '--debug', '-we', '-W', '-A', '-r']

    for option in test_input :
        assert option in test_output

