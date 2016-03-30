import numpy as np
from numpy import ndarray
from pprint import pprint
from unittest import TestCase

import cmfj


class TestCmfj(TestCase):

    def test_boltzmann(self):
        array = np.array([1,2,3
                          ])
        pprint(cmfj.stephanboltzmann(array))