import sys
sys.path.insert(0,"../../site-packages")

import unittest

from headflow import ParamDict, NSProblem, NSPostProcessor, NSScheme, NSSolver

class MockProblem(NSProblem):
    def __init__(self, params):
        NSProblem.__init__(self, params)

    @classmethod
    def default_problem_params(cls):
        return ParamDict(p=1)

class MockPostProcessor(NSPostProcessor):
    def __init__(self, params):
        NSPostProcessor.__init__(self, params)

    @classmethod
    def default_postprocessor_params(cls):
        return ParamDict(pp=1)

class MockScheme(NSScheme):
    def __init__(self, params):
        NSScheme.__init__(self, params)

    @classmethod
    def default_scheme_params(cls):
        return ParamDict(s=1)

class TestNSSolver(unittest.TestCase):
    def test_base_class_data_flow(self):
        # TODO: Set up mock classes that validate the generic data flow between the problem, scheme and post processor classes

        problem = MockProblem({})
        self.assertEqual(1, problem.params.p)

        postproc = MockPostProcessor({})
        self.assertEqual(1, postproc.params.pp)

        scheme = MockScheme({})
        self.assertEqual(1, scheme.params.s)

        solver = NSSolver(problem, postproc, scheme, {})

	#solver.solve()

        self.assertEqual(1, 1)

