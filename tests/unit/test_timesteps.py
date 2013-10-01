#!/usr/bin/env python
"""
Tests of timestepping utilities.
"""

import unittest

from headflow import NSProblem
from headflow.core.adaptivetimestepping import AdaptiveTimestepping
from headflow.core.timesteps import compute_regular_timesteps

class MockProblem(NSProblem):
    def __init__(self, params):
        NSProblem.__init__(self, params)

class TestTimestepComputation(unittest.TestCase):
    def test_compute_timesteps(self):
        problem = MockProblem({'dt': 1e-2, 'T': 1.1, 'T0': 0.2})

        dt, ts, start_timestep = compute_regular_timesteps(problem)

        self.assertEqual(dt, 1e-2)
        self.assertEqual(ts[0], 0.2)
        self.assertAlmostEqual(ts[1]-ts[0], dt)
        self.assertEqual(len(ts), 91)
        self.assertAlmostEqual(ts[-1], ts[0]+dt*(len(ts)-1))


    def test_adaptive_timestepper(self):
        # Configure timestepper
        p = AdaptiveTimestepping.default_params()
        p.T0 = 1.0
        p.T = 2.0
        p.dt = 0.5
        p.increase_on_positive_gradient = True
        ts = AdaptiveTimestepping(p)

        # --- Start timestep 0
        # Check initial time
        t = ts.time()
        timestep = ts.current_timestep()
        self.assertEqual(float(t), 1.0)
        self.assertEqual(timestep, 0)

        # Store values at end of accepted timestep
        curr_t = float(t)
        curr_dt = p.dt
        curr_timestep = timestep
        self.assertEqual(curr_t, 1.0)
        self.assertEqual(curr_dt, 0.5)
        self.assertEqual(curr_timestep, 0)

        # --- Start timestep 1
        # Check first timestep
        timestep = ts.next()
        self.assertEqual(float(t), p.T0 + p.dt)
        self.assertEqual(timestep, 1)

        # Don't trigger adaption (just below max)
        timestep_modified, dt = ts.adjusted_timestep(p.max_value-1)
        self.assertFalse(timestep_modified)
        self.assertEqual(dt, p.dt)

        # Store values at end of accepted timestep
        curr_t = float(t)
        curr_dt = dt
        curr_timestep = timestep
        self.assertEqual(curr_t, p.T0 + p.dt)
        self.assertEqual(curr_dt, p.dt)
        self.assertEqual(curr_timestep, 1)

        # --- Start timestep 2
        # Check second timestep
        timestep = ts.next()
        self.assertEqual(float(t), p.T0 + 2*p.dt)
        self.assertEqual(timestep, 2)

        # Trigger adaption from high value
        timestep_modified, dt = ts.adjusted_timestep(p.max_value+1)
        self.assertTrue(timestep_modified)
        self.assertEqual(dt, p.dt/p.decrease_factor)
        curr_dt = dt

        # Check that time has been reset
        self.assertEqual(float(t), p.T0 + p.dt)
        self.assertEqual(ts.current_timestep(), 1)

        # "Retry" second timestep
        timestep = ts.next()
        self.assertEqual(float(t), p.T0 + (1+1.0/p.decrease_factor)*p.dt)
        self.assertEqual(timestep, 2)

        # Don't trigger adaption (just above min)
        timestep_modified, dt = ts.adjusted_timestep(p.min_value+1)
        self.assertFalse(timestep_modified)
        self.assertEqual(dt, curr_dt)

        # Store values at end of accepted timestep
        curr_t = float(t)
        curr_dt = dt
        curr_timestep = timestep

        # --- Start timestep 3
        # Check third timestep
        timestep = ts.next()
        self.assertEqual(float(t), curr_t + curr_dt)
        self.assertEqual(timestep, 3)

        # Trigger adaption from low value
        timestep_modified, dt = ts.adjusted_timestep(p.min_value-2)
        self.assertTrue(timestep_modified)
        self.assertEqual(dt, curr_dt*p.increase_factor)

        # Check that time has NOT been reset
        self.assertEqual(float(t), curr_t + curr_dt)
        self.assertEqual(ts.current_timestep(), 3)

        # Store values at end of accepted timestep
        curr_t = float(t)
        curr_dt = dt
        curr_timestep = timestep

        # --- Start timestep 4
        # TODO: Test behaviour when staying below min_value but increasing
        # TODO: Test that dt does not increase more often than p.min_timesteps_between_increases
        # TODO: Other situations to test?
