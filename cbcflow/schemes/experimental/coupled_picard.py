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

__author__ = "Martin Alnaes <martinal@simula.no>"
__date__ = "2013-05-22"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from cbcflow.core.nsscheme import *
from cbcflow.utils.schemes import compute_regular_timesteps,assign_ics_mixed, make_velocity_bcs, make_rhs_pressure_bcs
from cbcflow.utils.core import NSSpacePoolMixed


class CoupledPicard(NSScheme):
    "Incremental pressure-correction scheme with penalty terms for boundary conditions."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        fc = ParamDict(
            optimize=True,
            cpp_optimize=True,
            cpp_optimize_flags="-O2", #"-O3 -march=native -fno-math-errno",
            quadrature_degree="auto",
            )
        params.update(
            # Default to P2-P1 (Taylor-Hood)
            u_degree = 2,
            p_degree = 1,

            picard_tolerance=1e-13,
            picard_max_iterations=20,
            picard_error_on_nonconvergence=False,

            equation=1,

            reuse_lu_data=True,
            verbose=False,
            form_compiler_parameters=fc,
            )
        return params

    def solve(self, problem, update):
        # Spatial parameters
        mesh = problem.mesh
        n  = FacetNormal(mesh)
        dx = problem.dx
        ds = problem.ds

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

        # Function spaces
        spaces = NSSpacePoolMixed(mesh, self.params.u_degree, self.params.p_degree)
        W = spaces.W

        # Test and trial functions
        u, p = TrialFunctions(W)
        v, q = TestFunctions(W)

        # Solution functions
        up0 = Function(W, name="up0") # Previous timestep
        up1 = Function(W, name="up1") # Last iterate of fixed-point in current timestep
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

        # Picard linearization of Navier-Stokes, F = a*u - L = 0
        eqchoice = self.params.equation
        if eqchoice == 1:
            # Not scaled by k
            a = (
                dot((1.0/k)*u + (grad(u)*u1), v)*dx()
                + nu*inner(grad(u), grad(v))*dx()
                - p*div(v)*dx() - q*div(u)*dx()
                )
            L = dot((1.0/k)*u0 + f, v)*dx() + Lbc
        if eqchoice == 2:
            # Scaled by k (smaller residual, nonlinear solver hits absolute stopping criteria faster)
            a = (
                dot(u + k*(grad(u)*u1), v)*dx()
                + (k*nu)*inner(grad(u), grad(v))*dx()
                - (k*p)*div(v)*dx() - (k*div(u))*q*dx()
                )
            L = dot(u0 + k*f, v)*dx() + k*Lbc
        if eqchoice == 3:
            # Stokes
            a = (
                (1.0/k)*dot(u, v)*dx()
                + nu*inner(grad(u), grad(v))*dx()
                - p*div(v)*dx() - q*div(u)*dx()
                )
            L = (1.0/k)*dot(u0, v)*dx() + dot(f, v)*dx() + Lbc
        # Build residual form from a, L
        F = action(a, up1) - L

        # Create solver
        picard_problem = NonlinearVariationalProblem(F, up1, bcu, J=a,
                                                     form_compiler_parameters=self.params.form_compiler_parameters)
        solver = NonlinearVariationalSolver(picard_problem)

        # Turn on the noise
        if self.params.verbose:
            solver.parameters["lu_solver"]["report"] = True
            solver.parameters["lu_solver"]["verbose"] = True
            solver.parameters["newton_solver"]["report"] = True

        # Speed up solvers with reuse
        if self.params.reuse_lu_data:
            solver.parameters["lu_solver"]["reuse_factorization"] = True
            solver.parameters["lu_solver"]["same_nonzero_pattern"] = True
            solver.parameters["reset_jacobian"] = False

        # Define stricter stopping criteria for "newton" (actually fixed point) solver
        solver.parameters["newton_solver"]["absolute_tolerance"] = self.params.picard_tolerance
        solver.parameters["newton_solver"]["relative_tolerance"] = self.params.picard_tolerance
        solver.parameters["newton_solver"]["maximum_iterations"] = self.params.picard_max_iterations
        solver.parameters["newton_solver"]["error_on_nonconvergence"] = self.params.picard_error_on_nonconvergence

        # Call update() with initial conditions
        update(u0, p0, float(t), start_timestep, spaces)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)

            # Solve for up1
            solver.solve()

            # Rotate functions for next timestep
            up0.assign(up1)

            # Update postprocessing
            # TODO: Pass controls and observations here?
            update(u0, p0, float(t), timestep, spaces)

        # Make sure annotation gets that the timeloop is over
        finalize_time(t)

        # Return some quantities from the local namespace
        states = (u0, p0)
        namespace = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return namespace
