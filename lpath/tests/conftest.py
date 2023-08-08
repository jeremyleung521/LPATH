"""
Fixtures for other tests.
"""

import pytest
import lpath
from lpath import argparser
from tqdm.auto import tqdm


@pytest.fixture
def create_ref_parser():
    parser = argparser.add_all_args()
    test_output = [action for actions in parser._actions for action in actions.option_strings]

    return parser, test_output


@pytest.fixture
def create_calc_dist_env():
    pbar = tqdm(range(1), leave=True, disable=True, delay=1)
    state_dict = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: '!'}

    return state_dict, pbar
