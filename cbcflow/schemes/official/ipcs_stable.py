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
r"""
This scheme follows the same logic as in :class:`.IPCS`, but with a few notable exceptions.

A parameter :math:`\theta` is added to the diffusion and convection terms,
allowing for different evaluation of these, and the convection is handled semi-implicitly:

.. math::
    \frac{1}{\Delta t}\left( \tilde{u}^{n+1}-u^{n} \right)-
    \nabla\cdot\nu\nabla \tilde{u}^{n+\theta}+
    u^*\cdot\nabla \tilde{u}^{n+\theta}+\nabla p^{n}=f^{n+1},
    
where

.. math::
    u^* = \frac{3}{2}u^n - \frac{1}{2}u^{n-1}, \\
    \tilde{u}^{n+\theta} = \theta \tilde{u}^{n+1}+\left(1-\theta\right)u^n.

This convection term is unconditionally stable, and with :math:`\theta=0.5`,
this equation is second order in time and space [1]_.


In addition, the solution process is significantly faster by solving for each of the
velocity components separately, making for D number of smaller linear systems compared
to a large system D times the size.



.. [1] Simo, J. C., and F. Armero. *Unconditional stability and long-term behavior
    of transient algorithms for the incompressible Navier-Stokes and Euler equations.*
    Computer Methods in Applied Mechanics and Engineering 111.1 (1994): 111-154.

"""


from __future__ import division


from cbcflow.core.nsscheme import *
from cbcflow.utils.common import is_periodic, cbcflow_log
from cbcflow.utils.schemes import (RhsGenerator,
                                   compute_regular_timesteps,
                                   assign_ics_segregated,
                                   make_segregated_velocity_bcs,
                                   make_pressure_bcs,
                                   make_penalty_pressure_bcs)
from cbcflow.utils.core import NSSpacePoolSegregated


class IPCS_Stable(NSScheme):
    "Incremental pressure-correction scheme, fast and stable version."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            theta = 0.5,
            rebuild_prec_frequency = 1e16,
            u_tent_prec_structure = "same_nonzero_pattern",
            u_tent_solver_parameters = {},
            p_corr_solver_parameters = {},
            u_corr_solver_parameters = {},
            )
        return params

    def solve(self, problem, update, timer):

        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        dims = range(mesh.topology().dim())
        theta = self.params.theta

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

        # Define function spaces
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)
        V = spaces.V
        U = spaces.U
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(U)
        q = TestFunction(Q)
        u = TrialFunction(U)
        p = TrialFunction(Q)

        # Functions
        u0 = as_vector([Function(U, name="u0_%d"%d) for d in dims]) # u^{n-1}
        u1 = as_vector([Function(U, name="u1_%d"%d) for d in dims]) # u^n
        u2 = as_vector([Function(U, name="u2_%d"%d) for d in dims]) # u^{n+1}
        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")
        p2 = Function(Q, name="p2")

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Apply initial conditions and use it as initial guess
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_segregated(u0, p0, spaces, ics)
        for d in dims: u1[d].assign(u0[d])
        for d in dims: u2[d].assign(u1[d])
        p1.assign(p0)
        p2.assign(p1)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u1, p1, t, controls)
        bcu = make_segregated_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Remove boundary stress term is problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        rho = float(problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        # Tentative velocity step. Crank-Nicholson time-stepping is used for diffusion and convection.
        a1 = (1/k) * inner(v, u) * dx()
        a2 = inner(grad(v), nu*grad(u)) * dx()

        # Convection linearized as in Simo/Armero (1994)
        a_conv = v * sum((1.5*u1[r] - 0.5*u0[r]) * u.dx(r) for r in dims) * dx()
        if theta < 1.0:
            # Set a_conv to match rhs theta-weighting for RHSGenerator
            a_conv = Constant(1-theta)*a_conv
            Kconv_axpy_factor = theta/(1-theta)
        else:
            Kconv_axpy_factor = 1.0
        Kconv = Matrix() # assembled from a_conv in the time loop

        # Create the static part of the coefficient matrix for the tentative
        # velocity step. The convection is added in the time loop. We use a
        # single matrix for all dimensions, which means that the BCs must match
        # (they must apply to the full vector space at the vertex, not a
        # subspace.
        A_u_tent = assemble(a1+theta*a2)

        # Create matrices for generating the RHS
        B = assemble(a1-(1-theta)*a2)
        M = assemble(v*u*dx())

        # Define how to create the RHS for the tentative velocity. The RHS is
        # (of course) different for each dimension.
        rhs_u_tent = [None]*len(dims)
        for d in dims:
            C = assemble(-v*p*n[d]*ds() + v.dx(d)*p*dx())
            rhs_u_tent[d] = RhsGenerator(U)
            rhs_u_tent[d] += B, u1[d]
            rhs_u_tent[d] += C, p1
            rhs_u_tent[d] += M, f[d]
            if theta < 1.0:
                rhs_u_tent[d] -= Kconv, u1[d]

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx())
        rhs_p_corr = RhsGenerator(Q)
        rhs_p_corr += A_p_corr, p1
        Ku = [None]*len(dims)
        for d in dims:
            Ku[d] = assemble(-(1/k)*q*u.dx(d)*dx()) # TODO: Store forms in list, this is copied below
            rhs_p_corr += Ku[d], u2[d]

        # Velocity correction. Like for the tentative velocity, a single LHS is used.
        A_u_corr = M
        rhs_u_corr = [None]*len(dims)
        Kp = [None]*len(dims)
        for d in dims:
            Kp[d] = assemble(-k*inner(v, grad(p)[d])*dx())
            rhs_u_corr[d] = RhsGenerator(U)
            rhs_u_corr[d] += M, u2[d]
            rhs_u_corr[d] += Kp[d], p2
            rhs_u_corr[d] -= Kp[d], p1

        # Apply BCs to LHS
        for bc in bcu:
            bc[0].apply(A_u_tent)
            bc[0].apply(A_u_corr)

        for bc in bcp:
            bc.apply(A_p_corr)

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

        for A,S in [(A_u_tent, solver_u_tent), (A_p_corr, solver_p_corr), (A_u_corr, solver_u_corr)]:
            S.set_operator(A)
            if 'preconditioner' in S.parameters:
                S.parameters['preconditioner']['structure'] = 'same'
                
        solver_u_tent.parameters.update(self.params.u_tent_solver_parameters)
        solver_p_corr.parameters.update(self.params.p_corr_solver_parameters)
        solver_u_corr.parameters.update(self.params.u_corr_solver_parameters)

        timer.completed("problem initialization")

        # Call update() with initial conditions
        update(u1, p1, float(t), start_timestep, spaces)
        timer.completed("initial postprocessor update")

        # Time loop
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u1, p1, t, timestep, bcs, observations, controls)
            timer.completed("problem update")
            
            p1.vector()[:] *= 1.0/rho

            # Assemble the u-dependent convection matrix. It is important that
            # it is assembled into the same tensor, because the tensor is
            # also stored in rhs. (And it's faster).
            if Kconv.size(0) == 0:
                # First time, just assemble normally
                assemble(a_conv, tensor=Kconv, reset_sparsity=True)
            else:
                # Subtract the convection for previous time step before re-assembling Kconv
                A_u_tent.axpy(-Kconv_axpy_factor, Kconv, True)
                assemble(a_conv, tensor=Kconv, reset_sparsity=False)

            # Either zero BC rows in Kconv, or re-apply BCs to A_u_tent after
            # the axpy (it doesn't matter which)
            #for bc in bcu:
            #    bc[0].zero(Kconv)

            A_u_tent.axpy(Kconv_axpy_factor, Kconv, True)
            for bc in bcu:
                bc[0].apply(A_u_tent)
            timer.completed("u_tent assemble convection & construct lhs")

            # Check if preconditioner is to be rebuilt
            if timestep % self.params.rebuild_prec_frequency == 0 and 'preconditioner' in solver_u_tent.parameters:
                solver_u_tent.parameters['preconditioner']['structure'] = self.params.u_tent_prec_structure

            # Compute tentative velocity step
            for d in dims:
                b = rhs_u_tent[d]()

                for bc in bcu: bc[d].apply(b)
                timer.completed("u_tent construct rhs")
                iter = solver_u_tent.solve(u2[d].vector(), b)
                
                # Preconditioner is the same for all three components, so don't rebuild several times
                if 'preconditioner' in solver_u_tent.parameters:
                    solver_u_tent.parameters['preconditioner']['structure'] = "same"

                timer.completed("u_tent solve (%s, %d dofs)"%(', '.join(self.params.solver_u_tent), b.size()), {"iter": iter})

            # Pressure correction
            b = rhs_p_corr()
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            for bc in bcp:
                b *= rho
                
                bc.apply(b)
                b *= 1.0/rho
            timer.completed("p_corr construct rhs")

            iter = solver_p_corr.solve(p2.vector(), b)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p2.vector())
            timer.completed("p_corr solve (%s, %d dofs)"%(', '.join(solver_p_params), b.size()), {"iter": iter})

            # Velocity correction
            for d in dims:
                b = rhs_u_corr[d]()
                for bc in bcu: bc[d].apply(b)
                timer.completed("u_corr construct rhs")

                iter = solver_u_corr.solve(u2[d].vector(), b)
                timer.completed("u_corr solve (%s, %d dofs)"%(', '.join(self.params.solver_u_corr), b.size()),{"iter": iter})

             # Rotate functions for next timestep
            for d in dims: u0[d].assign(u1[d])
            for d in dims: u1[d].assign(u2[d])
            p1.assign(p2)

            p1.vector()[:] *= rho
            
            # Update postprocessing
            update(u1, p1, float(t), timestep, spaces)
            timer.completed("updated postprocessing (completed timestep)")
            timer.increment()

        # Return some quantities from the local namespace
        states = (u1, p1)
        namespace = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return namespace
