from __future__ import division

__author__ = "Martin Alnaes <martinal@simula.no>"
__date__ = "2013-06-06"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ..dol import *

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