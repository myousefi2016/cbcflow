"""
Tests of timestepping utilities.
"""

import unittest

from headflow import NSProblem
from headflow.core.timesteps import compute_regular_timesteps

class MockProblem(NSProblem):
    def __init__(self, params):
        NSProblem.__init__(self, params)

class TestTimestepComputation(unittest.TestCase):
    def test_compute_timesteps(self):
        problem = MockProblem({'dt': 1e-2, 'T': 0.9})

        dt, t0, ts = compute_regular_timesteps(problem)

        self.assertEqual(dt, 1e-2)
        self.assertEqual(t0, 0.0)
        self.assertAlmostEqual(ts[1]-ts[0], dt)
        self.assertAlmostEqual(ts[-1], dt*len(ts))
        self.assertEqual(len(ts), 90)
