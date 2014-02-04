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

__author__ = "Oyvind Evju <oyvinev@simula.no>"
__date__ = "2013-04-30"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from cbcflow.core.nsscheme import *
from cbcflow.core.rhsgenerator import *
from cbcflow.core.utils import Timer, epsilon, sigma, is_periodic
from cbcflow.core.timesteps import compute_regular_timesteps
from cbcflow.core.schemeutils import assign_ics_mixed, make_velocity_bcs, make_rhs_pressure_bcs
from cbcflow.core.spaces import NSSpacePoolMixed


class CoupledNonLinear(NSScheme):
    "Coupled scheme with fixed-point iterations on the convection term. NB: Direct solver!"

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P2-P1 (Taylor-Hood)
            u_degree = 2,
            p_degree = 1,

            theta=0.5, # 0.5: Crank-Nicholson, 1.0: Backward Euler, 0.0: Forward Euler

            fixed_point_tolerance=1e-6,
            max_fixed_point_iterations=500,
            )
        return params

    def solve(self, problem, update, restart=None):
        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

        # Define function spaces
        spaces = NSSpacePoolMixed(mesh, self.params.u_degree, self.params.p_degree)
        V = spaces.V
        Q = spaces.Q
        W = spaces.W

        # Test and trial functions
        v, q = TestFunctions(W)
        u, p = TrialFunctions(W)

        # Functions
        up0 = Function(W, name="up0") # Previous timestep
        up1 = Function(W, name="up1") # Current timestep
        upk = Function(W, name="upk") # Previous iterate
        u0, p0 = split(up0)
        u1, p1 = split(up1)
        uk, pk = split(upk)

        # Get problem specific functions
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Apply initial conditions and use it as initial guess
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_mixed(up0, spaces, ics)
        up1.assign(up0)
        upk.assign(up1)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_velocity_bcs(problem, spaces, bcs)
        Lbc = make_rhs_pressure_bcs(problem, spaces, bcs, v)

        # Remove boundary stress term if problem is periodic
        #beta = Constant(0) if is_periodic(bcp) else Constant(1)
        beta = 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))
        theta = self.params.theta

        # Variational forms
        a1 = ( (1/k)*dot(u,v)*dx()
               + theta*2*nu*inner(epsilon(u), epsilon(v))*dx()
               - theta*p*div(v)*dx()
               - theta*q*div(u)*dx() )
        #a1 -= beta*theta*nu*inner(grad(u).T*n, v) * ds() # Unsure about the boundary terms
        a2 = theta*dot(grad(u)*u1, v)*dx()

        L = ( (1/k)*dot(u0,v)*dx()
              - (1-theta)*dot(grad(u0)*u0, v)*dx()
              - (1-theta)*2*nu*inner(epsilon(u0), epsilon(v))*dx()
              + (1-theta)*p0*div(v)*dx()
              + (1-theta)*q*div(u0)*dx()
              + Lbc
              + dot(f,v)*dx() ) # TODO: Should apply theta rule to f as well in principle
        #L += beta*(1-theta)*nu*inner(grad(u).T*n, v) * ds() # Unsure about the boundary terms

        # Residual form
        F = action(a1+a2, up1) - L

        # Preassemble matrices
        A = Matrix()
        A1 = assemble(a1)
        A2 = assemble(a2)
        A.assign(A1+A2)

        b = assemble(L)

        # Call update() with initial conditions
        update(u0, p0, float(t), start_timestep, spaces)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)

            # Fixed point loop
            count = 0
            err = 1e16
            errors = []
            while err > self.params.fixed_point_tolerance:
                count += 1
                if count > self.params.max_fixed_point_iterations:
                    raise Exception("Fixed point iteration did not converge, errors were %s." % str(errors))

                # Remember last iterate for convergence check
                upk.assign(up1)

                # TODO: For resistance BCs, do we want to update bcs for each iterate?

                # Reassemble rhs
                assemble(L, tensor=b) # TODO: With Patricks Picard formulation, I think this should be -F?
                for bc in bcu: bc.apply(b)

                # Reassemble convection term
                assemble(a2, tensor=A2, reset_sparsity=False)

                # Compute matrix from linear and nonlinear terms
                A.assign(A1+A2)
                for bc in bcu: bc.apply(A)

                # TODO: Preconditioned iterative solver
                #solve(A, up1.vector(), b, "gmres", "ilu")
                solve(A, up1.vector(), b)

                # Compute relative change in up1 since last iterate upk
                err = norm(up1.vector()-upk.vector()) / norm(up1.vector())
                errors.append(err)

            # Update last timestep
            up0.assign(up1)

            print "Fixed point iteration converged in %d iterations. (err=%.4e)" %(count, err)
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