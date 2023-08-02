"""
Unit and regression test for the lpath package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
import lpath.match


def test_lpath_match_imported():
    """
    Sample test. This will always pass so long as import statements worked.

    """
    assert "lpath.match" in sys.modules


class Test_tostr():
    """
    Test to see if tostr() works properly

    """

    def test_tostr_fail():
        """
        Test to see if tostr() fails as expected.

        """
        test_string = lpath.match.tostr(None)

        assert test_string == None



    def test_tostr_success():
        """
        Test to see if tostr() fails as expected.

        """
        input_string = 'abc'
        test_output = 'abc'
        test_string = lpath.match.tostr(input_string)

        assert test_string == test_output


def test_match_remove_none():
    """
    Test to see if string comprehension remove_consec_repeats() works
    with n = 0.

    """
    input_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'

    test_string = lpath.match.remove_consec_repeats(input_string, 0)

    assert test_string == test_output


def test_match_remove_single():
    """
    Test to see if string comprehension remove_consec_repeats() works
    with n = 1.

    """
    input_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'ABCDABABABABABABABCABCABCABC'

    test_string = lpath.match.remove_consec_repeats(input_string, 1)

    assert test_string == test_output


def test_match_remove_pair():
    """
    Test to see if string comprehension remove_consec_repeats() works
    with n = 2.

    """
    test_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'ABCDABCABCABCABC'

    test_string = lpath.match.remove_consec_repeats(test_string, 2)

    assert test_string == test_output


def test_match_remove_three():
    """
    Test to see if string comprehension remove_consec_repeats() works
    with n = 2.

    """
    test_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'ABCDABC'

    test_string = lpath.match.remove_consec_repeats(test_string, 3)

    assert test_string == test_output


def test_match_remove_single_old():
    """
    Test to see if string comprehension remove_consec_states() works.

    """
    input_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'ABCDABABABABABABABCABCABCABC'

    test_string = lpath.match.remove_consec_states(input_string)

    assert test_string == test_output


def test_match_remove_pairs_old():
    """
    Test to see if string comprehension remove_consec_states()
    followed by remove_consec_pairs() works.

    """
    test_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'ABCDABCABCABCABC'

    test_string = lpath.match.remove_consec_states(test_string)
    test_string = lpath.match.remove_consec_pairs(test_string)

    assert test_string == test_output
