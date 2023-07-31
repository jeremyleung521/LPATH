"""
Unit and regression test for the lpath package.
"""

# Import package, test suite, and other packages as needed
import sys
import lpath.match


def test_lpath_match_imported():
    """
    Sample test. This will always pass so long as import statements worked.

    """
    assert "lpath.match" in sys.modules


def test_lpath_match_remove_none():
    """
    Test to see if string comprehension remove_consec_repeats() works
    with n = 0.

    """
    input_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'

    test_string = lpath.match.remove_consec_repeats(input_string, 0)

    assert test_string == test_output


def test_lpath_match_remove_single():
    """
    Test to see if string comprehension remove_consec_repeats() works
    with n = 1.

    """
    input_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'ABCDABABABABABABABCABCABCABC'

    test_string = lpath.match.remove_consec_repeats(input_string, 1)

    assert test_string == test_output


def test_lpath_match_remove_pair():
    """
    Test to see if string comprehension remove_consec_repeats() works
    with n = 2.

    """
    test_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
    test_output = 'ABCDABCABCABCABC'

    test_string = lpath.match.remove_consec_repeats(test_string, 2)

    assert test_string == test_output
