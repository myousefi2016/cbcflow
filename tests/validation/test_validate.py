"""Scheme validation.

Work in progress.
"""

import unittest
import time

from headflow import NSProblem, NSScheme, NSSolver, all_schemes

class ValidationProblem(NSProblem):
    def __init__(self):
        # FIXME: Setup problem
        pass

class TestValidation(unittest.TestCase):
    @unittest.skip("Enable validation tests when problem is set up properly!")
    def test_validate_all_schemes(self):
        pass # TODO

