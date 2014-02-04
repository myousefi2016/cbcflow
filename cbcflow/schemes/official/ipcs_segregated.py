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

__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2012-02-15"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from cbcflow.core.nsscheme import *
from cbcflow.core.rhsgenerator import *
from cbcflow.core.utils import Timer, epsilon, sigma, is_periodic
from cbcflow.core.timesteps import compute_regular_timesteps
from cbcflow.core.schemeutils import (assign_ics_segregated,
                                make_segregated_velocity_bcs,
                                make_pressure_bcs,
                                make_penalty_pressure_bcs)
from cbcflow.core.spaces import NSSpacePoolSegregated


class SegregatedIPCS(NSScheme):
    "Segregated incremental pressure-correction scheme."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            )
        return params

    def solve(self, problem, update, restart=None):
        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        dims = range(mesh.topology().dim())

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

        # Define function spaces
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)
        U = spaces.U
        V = spaces.V
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(U)
        q = TestFunction(Q)
        u = TrialFunction(U)
        p = TrialFunction(Q)

        # Functions
        u0 = as_vector([Function(U, name="u0_%d"%d) for d in dims])
        u1 = as_vector([Function(U, name="u1_%d"%d) for d in dims])
        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Apply initial conditions and use it as initial guess
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_segregated(u0, p0, spaces, ics)
        for d in dims: u1[d].assign(u0[d])
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_segregated_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Remove boundary stress term is problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        # Tentative velocity
        F_u_tent = []
        r = Index()
        for d in dims:
            u_mean = 0.5 * (u + u0[d])
            u_diff = (u - u0[d])
            F_u_tent += [(1/k) * inner(v, u_diff) * dx()
                         + v * u0[d].dx(r)*u0[r] * dx()
                         + inner(grad(v), nu*grad(u_mean)) * dx()
                         - v.dx(d) * p0 * dx()
                         + v * p0 * n[d] * ds()
                         - v * f[d] * dx()]
        a_u_tent = [lhs(F) for F in F_u_tent]
        L_u_tent = [rhs(F) for F in F_u_tent]

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx()
        L_p_corr = inner(grad(q), grad(p0))*dx() - (1/k)*q*div(u1)*dx()

        # Velocity correction
        a_u_corr = [inner(v, u)*dx() for r in dims]
        L_u_corr = [v*u1[r]*dx() - k*inner(v, grad(p1-p0)[r])*dx() for r in dims]

        # Assemble matrices
        A_u_tent = [assemble(a) for a in a_u_tent]
        A_p_corr = assemble(a_p_corr)
        A_u_corr = [assemble(a) for a in a_u_corr]

        # Create solvers
        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0 or is_periodic(bcp):
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        solver_u_tent = LinearSolver(*self.params.solver_u_tent)
        solver_p_corr = LinearSolver(*solver_p_params)
        solver_u_corr = LinearSolver(*self.params.solver_u_corr)

        # Call update() with initial conditions
        update(u0, p0, float(t), start_timestep, spaces)

        timer = Timer(self.params.enable_timer)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Solve tentative velocity
            for d in dims:
                b = assemble(L_u_tent[d])
                for bc in bcu: bc[d].apply(A_u_tent[d], b)
                timer.completed("u_tent construct rhs")

                iter = solver_u_tent.solve(A_u_tent[d], u1[d].vector(), b)
                timer.completed("u_tent solve (%s, %d dofs, %d iter)" % (
                    ', '.join(self.params.solver_u_tent), b.size(), iter))

            # Pressure correction
            b = assemble(L_p_corr)
            if len(bcp) == 0 or is_periodic(bcp):
                normalize(b)
            else:
                for bc in bcp: bc.apply(A_p_corr, b)
            timer.completed("p_corr construct rhs")

            iter = solver_p_corr.solve(A_p_corr, p1.vector(), b)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            timer.completed("p_corr solve (%s, %d dofs, %d iter)" % (
                ', '.join(solver_p_params), b.size(), iter))

            # Velocity correction
            for d in dims:
                b = assemble(L_u_corr[d])
                for bc in bcu: bc[d].apply(A_u_corr[d], b)
                timer.completed("u_corr construct rhs")

                iter = solver_u_corr.solve(A_u_corr[d], u1[d].vector(), b)
                timer.completed("u_corr solve (%s, %d dofs, %d iter)" % (
                    ', '.join(self.params.solver_u_corr), b.size(), iter))

            # Rotate functions for next timestep
            for d in dims: u0[d].assign(u1[d])
            p0.assign(p1)

            # Update postprocessing
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