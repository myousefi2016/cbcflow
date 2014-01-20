"""Performance testing.

Work in progress.
"""

import unittest
import time

from cbcflow import NSProblem, NSScheme, NSSolver, all_schemes

class PerformanceProblem(NSProblem):
    def __init__(self):
        # FIXME: Setup problem
        # This problem should be:
        # - large enough to measure scheme performance
        # - parameterizable on problem size to measure scalablility
        pass

class TestPerformance(unittest.TestCase):
    @unittest.skip("Enable performance test when performance problem is set up properly!")
    def test_performance_of_all_schemes(self):
        "For each scheme, test that basic performance has not reduced."
        for Scheme in all_schemes:
            problem = PerformanceProblem()
            scheme = Scheme()
            solver = NSSolver(problem, scheme)
            t = time.time()
            solver.solve()
            t = time.time() - t

            # FIXME: Read reference time oldt from file, store time t to other file
            oldt = None

            if oldt is None:
                print "Missing reference time for scheme %s" % scheme.__class__.__name__
            else:
                self.assertLess(t, oldt*1.1)

    def test_all_schemes_scale_linearly(self):
        "For each scheme, test that cost increases reasonably close to linear with problem size."
        pass # TODO
