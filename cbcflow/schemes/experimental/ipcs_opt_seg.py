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
__date__ = "2011-11-11"
__copyright__ = "Copyright (C) 2011 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from cbcflow.core.nsscheme import *
from cbcflow.core.rhsgenerator import *
from cbcflow.core.utils import Timer, is_periodic, epsilon
from cbcflow.core.timesteps import compute_regular_timesteps
from cbcflow.core.schemeutils import (assign_ics_segregated,
                                make_segregated_velocity_bcs,
                                make_pressure_bcs,
                                make_penalty_pressure_bcs)
from cbcflow.core.spaces import NSSpacePoolSegregated

class SegregatedIPCS_Optimized(NSScheme):
    "Incremental pressure-correction scheme, optimized version."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
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
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)
        U = spaces.U
        V = spaces.V
        Q = spaces.Q
        dims = spaces.dims

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

        # Get initial conditions
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_segregated(u0, p0, spaces, ics)
        for d in dims: u1[d].assign(u0[d])
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_segregated_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        # Tentative velocity
        M  = assemble(inner(v, u) * dx())
        K1 = assemble((1/k) * inner(v, u) * dx())
        K2 = assemble(0.5 * inner(grad(v), nu*grad(u)) * dx())
        a_conv = -v * dot(grad(u),u0) * dx()
        Kconv = Matrix() # assembled from a_conv in the time loop

        A_u_tent = []
        rhs_u_tent = []
        for d in dims:
            A_u_tent.append(K1+K2) # Separate matrices, because they may have different BCs
            K3 = assemble(-v*p*n[d]*ds() + v.dx(d)*p*dx())

            rhs = RhsGenerator(U)
            rhs += K1, u0[d]
            rhs -= K2, u0[d]
            rhs += K3, p0
            rhs += M, f[d]
            rhs += Kconv, u0[d]

            rhs_u_tent.append(rhs)

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx())
        rhs_p_corr = RhsGenerator(Q)
        rhs_p_corr += A_p_corr, p0
        for r in dims:
            Ku = assemble(-(1/k)*q*u.dx(r)*dx())
            rhs_p_corr += Ku, u1[r]

        # Velocity correction
        A_u_corr = [M.copy() for r in dims]
        rhs_u_corr = []
        for r in dims:
            Kp = assemble(-k*inner(v, grad(p)[r])*dx())
            rhs = RhsGenerator(U)
            rhs += M, u1[r]
            rhs += Kp, p1
            rhs -= Kp, p0
            rhs_u_corr.append(rhs)

        # Apply BCs to matrices
        for bc in bcu:
            for d in dims:
                bc[d].apply(A_u_tent[d])
                bc[d].apply(A_u_corr[d])
        for bc in bcp:
            bc.apply(A_p_corr)

        # Create solvers
        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0 or is_periodic(bcp):
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        solver_u_tent = [LinearSolver(*self.params.solver_u_tent) for d in dims]
        solver_p_corr = LinearSolver(*solver_p_params)
        solver_u_corr = [LinearSolver(*self.params.solver_u_corr) for d in dims]

        for A,S in zip(A_u_tent, solver_u_tent) \
                + [(A_p_corr, solver_p_corr)] \
                + zip(A_u_corr, solver_u_corr):
            S.set_operator(A)
            if 'preconditioner' in S.parameters:
                S.parameters['preconditioner']['reuse'] = True

        # Call update() with initial conditions
        update(u0, p0, float(t), start_timestep, spaces)

        timer = Timer(self.params.enable_timer)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Assemble the u0-dependent convection matrix. It is important that
            # it is assembled into the same tensor, because it is stored in rhs.
            assemble(a_conv, tensor=Kconv, reset_sparsity=(Kconv.size(0)==0))
            timer.completed("u0 reassemble convection matrix")

            # Compute tentative velocity step
            for d, S, rhs, u1_comp in zip(dims, solver_u_tent, rhs_u_tent, u1):
                b = rhs()
                for bc in bcu: bc[d].apply(b)
                timer.completed("u0 construct rhs")

                iter = S.solve(u1_comp.vector(), b)
                timer.completed("u0 solve (%s, %d, %d)" % (
                    ', '.join(self.params.solver_u_tent), A.size(0), iter))

            # Pressure correction
            b = rhs_p_corr()
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            for bc in bcp: bc.apply(b)
            timer.completed("p1 construct rhs")

            iter = solver_p_corr.solve(p1.vector(), b)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            timer.completed("p1 solve (%s, %d, %d)" % (
                ', '.join(solver_p_params), A_p_corr.size(0), iter))

            # Velocity correction
            for d, S, rhs, u1_comp in zip(dims, solver_u_corr, rhs_u_corr, u1):
                b = rhs()
                for bc in bcu: bc[d].apply(b)
                timer.completed("u1 construct rhs")

                iter = S.solve(u1_comp.vector(), b)
                timer.completed("u1 solve (%s, %d, %d)" % (
                    ', '.join(self.params.solver_u_corr), A.size(0),iter))

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