import unittest

from headflow import NSProblem, NSScheme, NSSolver

class MockProblem(NSProblem):
    def __init__(self):
        pass

class MockScheme(NSScheme):
    def __init__(self):
        pass

class TestNSSolver(unittest.TestCase):
    def test_one_equals_one(self):
        problem = MockProblem()
        scheme = MockScheme()
        params = {}
        solver = NSSolver(problem, scheme, params)
        self.assertEqual(1, 1)
