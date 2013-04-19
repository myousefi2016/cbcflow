"""
Tests of the interaction between the NSSsolver
class and the user overloaded subclasses.
These tests should act as a documentation of
the control and data flow between the classes.
"""

import unittest

from headflow import ParamDict, NSProblem, NSPostProcessor, NSScheme, NSSolver

class MockProblem(NSProblem):
    def __init__(self, params):
        NSProblem.__init__(self, params)

    @classmethod
    def default_user_params(cls):
        return ParamDict(p=1)

class MockPostProcessor(NSPostProcessor):
    def __init__(self, params):
        NSPostProcessor.__init__(self, params)

    @classmethod
    def default_user_params(cls):
        return ParamDict(pp=2)

class MockScheme(NSScheme):
    def __init__(self, params):
        NSScheme.__init__(self, params, segregated=False)

    @classmethod
    def default_user_params(cls):
        return ParamDict(s=3)

class TestNSSolver(unittest.TestCase):
    def test_base_class_data_flow(self):
        # TODO: Set up mock classes that validate the generic data flow
        #       between the problem, scheme and post processor classes

        problem = MockProblem({})
        self.assertEqual(1, problem.params.p)

        postproc = MockPostProcessor({})
        self.assertEqual(2, postproc.params.pp)

        scheme = MockScheme({})
        self.assertEqual(3, scheme.params.s)

        solver = NSSolver(problem, postproc, scheme, {})

	#solver.solve()

        self.assertEqual(1, 1)
