#!/usr/bin/env py.test
"""
Tests of timestepping utilities.
"""

from cbcflow import NSProblem
from cbcflow.schemes.utils import compute_regular_timesteps
#from cbcflow.core.adaptivetimestepping import AdaptiveTimestepping

class MockProblem(NSProblem):
    def __init__(self, params):
        NSProblem.__init__(self, params)

def test_compute_timesteps():
    problem = MockProblem({'dt': 1e-2, 'T': 1.1, 'T0': 0.2})

    dt, ts, start_timestep = compute_regular_timesteps(problem)

    assert dt == 1e-2
    assert ts[0] == 0.2
    assert abs( (ts[1]-ts[0]) - (dt) ) < 1e-8
    #assert len(ts) == 91
    assert abs( (ts[-1]) - (ts[0]+dt*(len(ts)-1)) ) < 1e-8

def xtest_adaptive_timestepper():
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
    assert float(t) == 1.0
    assert timestep == 0

    # Store values at end of accepted timestep
    curr_t = float(t)
    curr_dt = p.dt
    curr_timestep = timestep
    assert curr_t == 1.0
    assert curr_dt == 0.5
    assert curr_timestep == 0

    # --- Start timestep 1
    # Check first timestep
    timestep = ts.next()
    assert float(t) == p.T0 + p.dt
    assert timestep == 1

    # Don't trigger adaption (just below max)
    timestep_modified, dt = ts.adjusted_timestep(p.max_value-1)
    assert not timestep_modified
    assert dt == p.dt

    # Store values at end of accepted timestep
    curr_t = float(t)
    curr_dt = dt
    curr_timestep = timestep
    assert curr_t == p.T0 + p.dt
    assert curr_dt == p.dt
    assert curr_timestep == 1

    # --- Start timestep 2
    # Check second timestep
    timestep = ts.next()
    assert float(t) == p.T0 + 2*p.dt
    assert timestep == 2

    # Trigger adaption from high value
    timestep_modified, dt = ts.adjusted_timestep(p.max_value+1)
    assert timestep_modified
    assert dt == p.dt/p.decrease_factor
    curr_dt = dt

    # Check that time has been reset
    assert float(t) == p.T0 + p.dt
    assert ts.current_timestep() == 1

    # "Retry" second timestep
    timestep = ts.next()
    assert float(t) == p.T0 + (1+1.0/p.decrease_factor)*p.dt
    assert timestep == 2

    # Don't trigger adaption (just above min)
    timestep_modified, dt = ts.adjusted_timestep(p.min_value+1)
    assert not timestep_modified
    assert dt == curr_dt

    # Store values at end of accepted timestep
    curr_t = float(t)
    curr_dt = dt
    curr_timestep = timestep

    # --- Start timestep 3
    # Check third timestep
    timestep = ts.next()
    assert float(t) == curr_t + curr_dt
    assert timestep == 3

    # Trigger adaption from low value
    timestep_modified, dt = ts.adjusted_timestep(p.min_value-2)
    assert timestep_modified
    assert dt == curr_dt*p.increase_factor

    # Check that time has NOT been reset
    assert float(t) == curr_t + curr_dt
    assert ts.current_timestep() == 3

    # Store values at end of accepted timestep
    curr_t = float(t)
    curr_dt = dt
    curr_timestep = timestep

    # --- Start timestep 4
    # TODO: Test behaviour when staying below min_value but increasing
    # TODO: Test that dt does not increase more often than p.min_timesteps_between_increases
    # TODO: Other situations to test?

