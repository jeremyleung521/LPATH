"""
Unit and regression test for the lpath package.
"""

# Import package, test suite, and other packages as needed
import sys
import lpath.match
import numpy


def test_lpath_match_imported():
    """
    Sample test. This will always pass so long as import statements worked.

    """
    assert "lpath.match" in sys.modules


class TestToStr:
    """
    Class to test to see if tostr() works properly

    """

    def test_tostr_fail(self):
        """
        Test to see if tostr() fails as expected.

        """
        test_string = lpath.match.tostr(None)

        assert test_string is None

    def test_tostr_success(self):
        """
        Test to see if tostr() fails as expected.

        """
        input_strings = ['abc', b'abc']
        test_output = 'abc'
        for input_string in input_strings:
            test_string = lpath.match.tostr(input_string)
            assert test_string == test_output


class TestCalcDist:
    """
    Class to test all the calc_dist_*() functions

    """
    def test_calc_dist(self, create_calc_dist_env):
        state_dict, pbar = create_calc_dist_env
        ref_score = 0.2727272727
        seq1 = numpy.asarray([0, 4, 4, 4, 1, 2, 3, 4, 4, 4])
        seq2 = numpy.asarray([0, 1, 2, 3, 5, 5, 5, 5, 5, 5])

        test_score = lpath.match.calc_dist(seq1, seq2, state_dict, pbar, 0)

        assert numpy.isclose(ref_score, test_score)


class TestReassignMethod:
    """
    Class to test all the determine_reassign() functionalities.

    """
    def test_identity(self):
        test_output = lpath.match.determine_reassign('reassign_identity')

        assert test_output == lpath.match.reassign_identity

    def test_statelabel(self):
        test_output = lpath.match.determine_reassign('reassign_statelabel')

        assert test_output == lpath.match.reassign_statelabel

    def test_custom(self):
        test_output = lpath.match.determine_reassign('reassign_custom')

        assert test_output == lpath.match.reassign_custom

    def test_segid(self):
        test_output = lpath.match.determine_reassign('reassign_segid')

        assert test_output == lpath.match.reassign_segid

    def test_userprovided(self):
        test_output = lpath.match.determine_reassign('lpath.match.reassign_custom')

        assert test_output == lpath.match.reassign_custom


class TestCondenseString:
    """
    Class to test to see if condense_string() works properly.

    """

    def test_match_remove_none(self):
        """
        Test to see if string comprehension condense_string() works
        with n = 0.

        """
        input_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
        test_output = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'

        test_string = lpath.match.condense_string(input_string, 0)

        assert test_string == test_output

    def test_match_remove_single(self):
        """
        Test to see if string comprehension condense_string() works
        with n = 1.

        """
        input_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
        test_output = 'ABCDABABABABABABABCABCABCABC'

        test_string = lpath.match.condense_string(input_string, 1)

        assert test_string == test_output

    def test_match_remove_pair(self):
        """
        Test to see if string comprehension condense_string() works
        with n = 2.

        """
        test_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
        test_output = 'ABCDABCABCABCABC'

        test_string = lpath.match.condense_string(test_string, 2)

        assert test_string == test_output

    def test_match_remove_three(self):
        """
        Test to see if string comprehension condense_string() works
        with n = 3.

        """
        test_string = 'AAABBBCCCDDDAABBAABBABABABABABCABCABCABC'
        test_output = 'ABCDABC'

        test_string = lpath.match.condense_string(test_string, 3)

        assert test_string == test_output

    def test_match_remove_three_t2(self):
        """
        Test to see if string comprehension condense_string() works
        with n = 3.

        """
        test_string = 'AAAAAAAABCABCABCBBBBBBBBCDCDAWERWER'
        test_output = 'ABCDAWER'

        test_string = lpath.match.condense_string(test_string, 3)

        assert test_string == test_output
