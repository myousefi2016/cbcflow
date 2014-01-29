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

from ..core.nsscheme import *
from ..core.utils import Timer, epsilon, sigma, is_periodic
from ..core.timesteps import compute_regular_timesteps
from ..core.schemeutils import (assign_ics_split,
                                make_velocity_bcs,
                                make_pressure_bcs,
                                make_penalty_pressure_bcs)
from ..core.spaces import NSSpacePoolSplit


class IPCS_Stabilized(NSScheme):
    "Incremental pressure-correction scheme stabilized with SUPG."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.replace(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            )
        params.update(
            theta=0.5,
            )
        return params

    def solve(self, problem, update):
        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[0])

        # Define function spaces
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

        # Get boundary conditions
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        # Stabilised test function
        h = CellSize(mesh)
        h_avg = (h('+') + h('-'))/2
        delta = 0.0*Constant(mesh.hmin())
        #delta = Constant(1e-4)
        w = v + delta*grad(v)*u0

        # Tentative velocity step
        theta = self.params.theta
        u_mean = theta*u + (1-theta)*u0
        u_diff = (u - u0)
        F_u_tent = ((1/k) * inner(w, u_diff) * dx()
                    + inner(w, grad(u_mean)*u0) * dx()
                    + inner(epsilon(w), sigma(u_mean, p0, nu)) * dx()
                    - beta * nu * inner(grad(u_mean).T*n, w) * ds()
                    + inner(w, p0*n) * ds()
                    - inner(w, f) * dx())

        a_u_tent = lhs(F_u_tent)
        L_u_tent = rhs(F_u_tent)

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx()
        L_p_corr = inner(grad(q), grad(p0))*dx() - (1/k)*q*div(u1)*dx()

        # Velocity correction
        F_corr = inner(v, u)*dx() - inner(v, u1)*dx() + k*inner(v, grad(p1-p0))*dx()
        #F_corr = inner(w, u)*dx() - inner(w, u1)*dx() + k*inner(w, grad(p1-p0))*dx()

        #a_u_corr = inner(v, u)*dx()
        #L_u_corr = inner(v, u1)*dx() - k*inner(v, grad(p1-p0))*dx()

        a_u_corr = lhs(F_corr)
        L_u_corr = rhs(F_corr)

        # Assemble matrices
        A_u_tent = assemble(a_u_tent)
        A_p_corr = assemble(a_p_corr)
        A_u_corr = assemble(a_u_corr)

        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0 or is_periodic(bcp):
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        # Call update() with initial conditions
        update(u0, p0, float(t), start_timestep, spaces)

        timer = Timer(self.params.enable_timer)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Compute tentative velocity step
            assemble(a_u_tent, tensor=A_u_tent)
            b = assemble(L_u_tent)
            for bc in bcu: bc.apply(A_u_tent, b)
            timer.completed("u1 construct rhs")

            iter = solve(A_u_tent, u1.vector(), b, *self.params.solver_u_tent)
            print "Iterations: ", iter
            timer.completed("u1 solve (%s, %d, %d)"%(', '.join(self.params.solver_u_tent), b.size(), iter))

            # Pressure correction
            b = assemble(L_p_corr)
            if len(bcp) == 0 or is_periodic(bcp):
                normalize(b)
            else:
                for bc in bcp: bc.apply(A_p_corr, b)
            timer.completed("p construct rhs")

            iter = solve(A_p_corr, p1.vector(), b, *solver_p_params)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            timer.completed("p solve (%s, %d, %d)"%(', '.join(solver_p_params), b.size(), iter))

            # Velocity correction
            assemble(a_u_corr, tensor=A_u_corr)
            b = assemble(L_u_corr)
            for bc in bcu: bc.apply(A_u_corr, b)
            timer.completed("u2 construct rhs")

            solver_params = self.params.solver_u_corr
            iter = solve(A_u_corr, u1.vector(), b, *solver_params)
            timer.completed("u2 solve (%s, %d, %d)"%(', '.join(solver_params), b.size(),iter))

            print "u-norm: ", norm(u1.vector())
            print "Pressure drop: ", max(p1.vector())-min(p1.vector())

            # Rotate functions for next timestep
            u0.assign(u1)
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
