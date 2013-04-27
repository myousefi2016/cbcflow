#!/usr/bin/env python
"""
Tests of timestepping utilities.
"""

import unittest
from init_test import init_test
init_test(__name__)

from headflow import NSProblem
from headflow.core.timesteps import compute_regular_timesteps

class MockProblem(NSProblem):
    def __init__(self, params):
        NSProblem.__init__(self, params)

class TestTimestepComputation(unittest.TestCase):
    def test_compute_timesteps(self):
        problem = MockProblem({'dt': 1e-2, 'T': 1.1, 'T0': 0.2})

        dt, ts = compute_regular_timesteps(problem)

        self.assertEqual(dt, 1e-2)
        self.assertEqual(ts[0], 0.2)
        self.assertAlmostEqual(ts[1]-ts[0], dt)
        self.assertEqual(len(ts), 91)
        self.assertAlmostEqual(ts[-1], ts[0]+dt*(len(ts)-1))

if __name__ == "__main__":
    unittest.main()
