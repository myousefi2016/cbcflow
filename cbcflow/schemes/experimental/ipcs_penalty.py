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

__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-02-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.
# Modified by Martin Alnaes, 2013. (added penalty bcs)

from cbcflow.core.nsscheme import *
from cbcflow.core.utils import Timer, epsilon, sigma
from cbcflow.core.timesteps import compute_regular_timesteps
from cbcflow.core.schemeutils import (assign_ics_split,
                                make_velocity_bcs,
                                make_pressure_bcs,
                                make_penalty_pressure_bcs)
from cbcflow.core.spaces import NSSpacePoolSplit


class PenaltyIPCS(NSScheme):
    "Incremental pressure-correction scheme with penalty terms for boundary conditions."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,

            # Penalty pressure bcs is necessary for optimization
            use_penalty_pressure_bcs=True,
            penalty_gamma=100.0,
            )
        return params

    def solve(self, problem, update):
        # Get spatial parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

        # Function spaces
        spaces = NSSpacePoolSplit(mesh, self.params.u_degree, self.params.p_degree)
        V = spaces.V
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u0 = Function(V, name="u0")
        u1 = Function(V, name="u1")
        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Apply initial conditions and use it as initial guess
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_split(u0, p0, spaces, ics)
        u1.assign(u0)
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_velocity_bcs(problem, spaces, bcs)
        if self.params.use_penalty_pressure_bcs:
            bcp = []
            a_pbc, L_pbc = make_penalty_pressure_bcs(problem, spaces, bcs,
                                                     self.params.penalty_gamma, q, p)
        else:
            bcp = make_pressure_bcs(problem, spaces, bcs)
            a_pbc, L_pbc = 0, 0

        # Problem parameters
        nu = Constant(problem.params.mu/problem.params.rho, name="nu")
        k  = Constant(dt, name="dt")
        f  = as_vector(problem.body_force(spaces, t))

        # Tentative velocity step forms
        u_mean = 0.5 * (u + u0)
        u_diff = (u - u0)
        F_u_tent = ((1/k) * inner(v, u_diff) * dx()
                    + inner(v, grad(u0)*u0) * dx()
                    + inner(epsilon(v), sigma(u_mean, p0, nu)) * dx()
                    - nu * inner(grad(u_mean).T*n, v) * ds()
                    + inner(v, p0*n) * ds()
                    - inner(v, f) * dx())
        a_u_tent = lhs(F_u_tent)
        L_u_tent = rhs(F_u_tent)

        # Pressure correction forms
        a_p_corr = inner(grad(q), grad(p))*dx() + a_pbc
        L_p_corr = inner(grad(q), grad(p0))*dx() - (1/k)*q*div(u1)*dx() + L_pbc

        # Velocity correction forms
        a_u_corr = inner(v, u)*dx()
        L_u_corr = inner(v, u1)*dx() - k*inner(v, grad(p1-p0))*dx()

        # Assemble matrices
        A_u_tent = assemble(a_u_tent)
        A_p_corr = assemble(a_p_corr)
        A_u_corr = assemble(a_u_corr)

        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0:
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        update(u0, p0, float(t), start_timestep, spaces)

        # Profiling object
        timer = Timer(self.params.enable_timer)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Compute tentative velocity step
            b = assemble(L_u_tent)
            for bc in bcu: bc.apply(A_u_tent, b)
            timer.completed("u1 construct rhs")
            print "u1 rhs norm:", norm(b)

            iter = solve(A_u_tent, u1.vector(), b, *self.params.solver_u_tent)
            timer.completed("u1 solve (%s, %d, %d)"%(', '.join(self.params.solver_u_tent), b.size(), iter))
            print "u1 norm:", norm(u1.vector())

            # Pressure correction
            b = assemble(L_p_corr)
            for bc in bcp: bc.apply(A_p_corr, b)
            timer.completed("p construct rhs")
            print "pressure rhs norm:", norm(b)

            iter = solve(A_p_corr, p1.vector(), b, *solver_p_params)
            timer.completed("p solve (%s, %d, %d)"%(', '.join(solver_p_params), b.size(), iter))
            print "pressure norm:", norm(p1.vector())

            # Velocity correction
            b = assemble(L_u_corr)
            for bc in bcu: bc.apply(A_u_corr, b)
            timer.completed("u2 construct rhs")
            print "u2 rhs norm:", norm(b)

            solver_params = self.params.solver_u_corr
            iter = solve(A_u_corr, u1.vector(), b, *solver_params)
            timer.completed("u2 solve (%s, %d, %d)"%(', '.join(solver_params), b.size(),iter))
            print "u2 norm:", norm(u1.vector())

            # Rotate functions for next timestep
            u0.assign(u1)
            p0.assign(p1)

            # Update postprocessing
            update(u0, p0, float(t), timestep, spaces)

        # Make sure annotation gets that the timeloop is over
        finalize_time(t)

        # Return some quantities from the local namespace
        states = (u0, p0)
        scope = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return scope