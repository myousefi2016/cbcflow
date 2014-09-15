# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.
from __future__ import division


from cbcflow.dol import as_vector, project, FacetNormal, DirichletBC, Constant, dot, Dn
from cbcpost.utils import cbc_warning
from numpy import linspace, zeros, append, arange

# --- Functions for timestepping

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
    timesteps = arange(T0, T+dt, dt)
    assert abs( (timesteps[1]-timesteps[0]) - dt) < 1e-8, "Timestep size does not match specified dt."
    assert timesteps[-1] >= T-dt/1e6, "Timestep range not including end time."
    
    if timesteps[-1] > T+dt/1e6:
        cbc_warning("End time for simulation does not match end time set for problem (T-T0 not a multiple of dt).")
    
    return dt, timesteps, problem.params.start_timestep

# --- Initial condition helper functions for schemes

def assign_ics_mixed(up0, spaces, ics):
    up = as_vector(list(ics[0]) + [ics[1]])
    up0.assign(project(up, spaces.W)) #, name="up0_init_projection"))

def assign_ics_split(u0, p0, spaces, ics):
    u = as_vector(list(ics[0]))
    p = ics[1]
    u0.assign(project(u, spaces.V)) #, name="u0_init_projection"))
    p0.assign(project(p, spaces.Q)) #, name="p0_init_projection"))

def assign_ics_segregated(u0, p0, spaces, ics):
    for d in spaces.dims:
        u0[d].assign(project(ics[0][d], spaces.U)) #, name="u0_%d_init_projection"%d))
    p0.assign(project(ics[1], spaces.Q)) #, name="p0_init_projection"))


# --- Boundary condition helper functions for schemes

def _domainargs(problem, D):
    "Helper function to pass domain args if necessary."
    if isinstance(D, int):
        return (problem.facet_domains, D)
    else:
        return (D,)

def make_velocity_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs
    bcu = [DirichletBC(spaces.Ubc[d], functions[d], *_domainargs(problem, region))
           for functions, region in bcu_raw
           for d in spaces.dims]
    return bcu

def make_segregated_velocity_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs
    bcu = [[DirichletBC(spaces.Ubc[d], functions[d], *_domainargs(problem, region))
            for d in spaces.dims]
           for functions, region in bcu_raw]
    return bcu

def make_pressure_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs
    bcp = [DirichletBC(spaces.Q, function, *_domainargs(problem, region))
           for function, region in bcp_raw]
    return bcp

def make_rhs_pressure_bcs(problem, spaces, bcs, v):
    bcu_raw, bcp_raw = bcs
    ds = problem.ds
    n = FacetNormal(problem.mesh)
    Lbc = -sum(dot(function*n, v)*ds(region) for (function, region) in bcp_raw)
    return Lbc

def make_penalty_pressure_bcs(problem, spaces, bcs, gamma, test, trial):
    bcu_raw, bcp_raw = bcs

    # Define trial and test functions
    p = trial
    q = test
    Q = spaces.Q
    ds = problem.ds

    # Define Nitche discretization constants
    gamma = Constant(gamma, name="gamma")
    hE = Q.mesh().ufl_cell().max_facet_edge_length

    # The Nietche terms to integrate
    a_dirichlet = (gamma/hE)*(p*q) - Dn(p)*q - p*Dn(q)
    L_dirichlet = (gamma/hE)*q - Dn(q)

    # Collect Nietche terms for each subboundary
    a, L = [], []
    for pbc, D in bcp_raw:
        a += [    a_dirichlet*ds(D)]
        L += [pbc*L_dirichlet*ds(D)]

    # Accumulate and return both sides
    a_pbc = sum(a)
    L_pbc = sum(L)
    return a_pbc, L_pbc
