import unittest

from headflow import NSProblem, NSScheme, NSSolver

# FIXME: Set up a small collection of representative problems for regression testing
# FIXME: Make some reusable code for regression tests, just swapping problems and schemes

class Channel(NSProblem):
    def __init__(self):
        pass

class TestRegression(unittest.TestCase):
    @unittest.skip("Set up some test problems for regression testing!")
    def test_channel_pressure_drop(self):
        problem = ExampleProblem()
        scheme = SomeScheme()
        solver = NSSolver(problem, scheme)
        namespace = solver.solve()

        computed = some_functional_computation_from_namespace
        reference = read_from_reference_file

        self.assertAlmostEqual(computed, reference)

