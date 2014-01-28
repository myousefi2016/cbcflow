from __future__ import division
__author__ = "Martin Alnaes <martinal@simula.no>"
__date__ = "2013-04-26"
__copyright__ = "Copyright (C) 2013-2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from numpy import linspace, zeros, append

def iround(x):
    return int(x+0.5)

def compute_regular_timesteps(problem):
    """Compute fixed timesteps for problem.

    Returns (dt, t0, timesteps), where timesteps does not include t0.
    """
    # Get the time range for the problem
    T0 = problem.params.T0
    T = problem.params.T

    # The timestep must be given in the problem
    dt = problem.params.dt

    # Compute regular timesteps, not including t0
    num_intervals = iround((T-T0)/dt)
    timesteps = linspace(T0, T, num_intervals+1)
    timesteps = append(zeros(problem.params.start_timestep), timesteps)
    
    return dt, timesteps, problem.params.start_timestep
