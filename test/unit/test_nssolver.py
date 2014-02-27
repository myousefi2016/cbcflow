#!/usr/bin/env py.test
"""
Tests of the interaction between the NSSsolver
class and the user overloaded subclasses.
These tests should act as a documentation of
the control and data flow between the classes.
"""

from cbcflow import ParamDict, NSProblem, NSPostProcessor, NSScheme, NSSolver

class MockProblem(NSProblem):
    def __init__(self, params):
        NSProblem.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(T=1.0)
        params.update(p=1)
        return params

class MockPostProcessor(NSPostProcessor):
    def __init__(self, params):
        NSPostProcessor.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSPostProcessor.default_params()
        params.update(pp=2)
        return params

class MockScheme(NSScheme):
    def __init__(self, params):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        #params.replace(u_degree=1, p_degree=1)
        params.update(s=3)
        return params

def test_base_class_data_flow():
    # TODO: Set up mock classes that validate the generic data flow
    #       between the problem, scheme and post processor classes

    problem = MockProblem({})
    assert 1 == problem.params.p

    postproc = MockPostProcessor({})
    assert 2 == postproc.params.pp

    scheme = MockScheme({})
    assert 3 == scheme.params.s

    solver = NSSolver(problem, postproc, scheme, {})

    #solver.solve()

    assert 1 == 1
