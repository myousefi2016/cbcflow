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


from cbcflow.core.nsscheme import *
from cbcflow.utils.schemes import (compute_regular_timesteps,
                                      assign_ics_mixed,
                                      make_velocity_bcs,
                                      make_rhs_pressure_bcs)
from cbcflow.utils.core import NSSpacePoolMixed


class Stokes(NSScheme):
    "Coupled solver for the transient Stokes problem."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P2-P1 (Taylor-Hood)
            u_degree=2,
            p_degree=1,
            verbose=False,
            )
        return params

    def solve(self, problem, update, timer):
        # Spatial parameters
        mesh = problem.mesh
        n  = FacetNormal(mesh)
        dx = problem.dx
        ds = problem.ds

        # Time parameters
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[0])

        # Function spaces
        spaces = NSSpacePoolMixed(mesh, self.params.u_degree, self.params.p_degree)
        W = spaces.W

        # Test and trial functions
        u, p = TrialFunctions(W)
        v, q = TestFunctions(W)

        # Solution functions
        up0 = Function(W, name="up0") # Previous timestep
        up1 = Function(W, name="up1") # Current timestep
        u0, p0 = split(up0)
        u1, p1 = split(up1)

        # Get problem specific functions
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)
        ics = problem.initial_conditions(spaces, controls)
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)

        # Apply initial conditions and use it as initial guess
        assign_ics_mixed(up0, spaces, ics)
        up1.assign(up0)

        # Make scheme-specific representation of bcs
        bcu = make_velocity_bcs(problem, spaces, bcs)
        Lbc = make_rhs_pressure_bcs(problem, spaces, bcs, v)

        # Problem parameters
        nu = Constant(problem.params.mu/problem.params.rho, name="nu")
        k  = Constant(dt, name="dt")
        f  = as_vector(problem.body_force(spaces, t))
        c0 = Constant(0.0, name="zero")

        # Transient Stokes forms
        a = (
            (1.0/k)*dot(u, v)*dx()
            + nu*inner(grad(u), grad(v))*dx()
            - p*div(v)*dx() - q*div(u)*dx()
            + c0*p*q*ds(1) # Hack to induce facet markers into a for SystemAssembler...
            )
        L = (1.0/k)*dot(u0, v)*dx() + dot(f, v)*dx() + Lbc
        # Build residual form from a, L
        F = action(a, up1) - L

        # Create solver
        #lvproblem = LinearVariationalProblem(a, L, up1, bcu)
        #solver = LinearVariationalSolver(lvproblem)

        # Set solver and prec TODO: Are these the best? Does it work without block prec?
        #solver.parameters["linear_solver"] = "lu"
        #solver.parameters["preconditioner"] = "amg"
        #solver.parameters["symmetric"] = False #True

        # TODO: Call update() with initial conditions in all schemes?
        update(u0, p0, float(t), start_timestep, spaces)

        # Profiling object

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)

            # Solve for up1
            #solver.solve()
            solve(a == L, up1, bcu)

            # Rotate functions for next timestep
            up0.assign(up1)

            # Update postprocessing
            # TODO: Pass controls and observations here?
            update(u0, p0, float(t), timestep, spaces)

        # Make sure annotation gets that the timeloop is over
        finalize_time(t)

        # Return some quantities from the local namespace
        # TODO: The design of this is a bit ad-hoc...
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
