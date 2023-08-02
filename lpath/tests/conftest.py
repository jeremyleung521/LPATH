"""
Fixtures for other tests.
"""

import pytest
import lpath

@pytest.fixture
def create_ref_parser():
    output = lpath.argparser.add_all_args()
    test_output = [action for actions in output._actions for action in actions.option_strings]

    return output, test_output