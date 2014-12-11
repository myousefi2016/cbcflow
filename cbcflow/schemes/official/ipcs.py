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
This incremental pressure correction scheme (IPCS) is an operator splitting scheme that
follows the idea of Goda [1]_.
This scheme preserves the exact same stability properties
as Navier-Stokes and hence does not introduce additional dissipation in the flow.

The idea is to replace the unknown pressure with an approximation. This is chosen as
the pressure solution from the previous solution.

The time discretization is done using backward Euler, the diffusion term is handled with Crank-Nicholson, and the convection is handled explicitly, making the
equations completely linear. Thus, we have a discretized version of the Navier-Stokes equations as

.. math:: \frac{1}{\Delta t}\left( u^{n+1}-u^{n} \right)-\nabla\cdot\nu\nabla u^{n+\frac{1}{2}}+u^n\cdot\nabla u^{n}+\frac{1}{\rho}\nabla p^{n+1}=f^{n+1}, \\
    \nabla \cdot u^{n+1} = 0,

where :math:`u^{n+\frac{1}{2}} = \frac{1}{2}u^{n+1}+\frac{1}{2}u^n.`

For the operator splitting, we use the pressure solution from the previous timestep as an estimation, giving an equation for a tentative velocity, :math:`\tilde{u}^{n+1}`:

.. math:: \frac{1}{\Delta t}\left( \tilde{u}^{n+1}-u^{n} \right)-\nabla\cdot\nu\nabla \tilde{u}^{n+\frac{1}{2}}+u^n\cdot\nabla u^{n}+\frac{1}{\rho}\nabla p^{n}=f^{n+1}.

This tenative velocity is not divergence free, and thus we define a velocity correction :math:`u^c=u^{n+1}-\tilde{u}^{n+1}`. Substracting the second equation from the first, we see that

.. math::
    \frac{1}{\Delta t}u^c-\frac{1}{2}\nabla\cdot\nu\nabla u^c+\frac{1}{\rho}\nabla\left( p^{n+1} - p^n\right)=0, \\
    \nabla \cdot u^c = -\nabla \cdot \tilde{u}^{n+1}.

The operator splitting is a first order approximation, :math:`O(\Delta t)`, so we can, without reducing the order of the approximation simplify the above to

.. math::
    \frac{1}{\Delta t}u^c+\frac{1}{\rho}\nabla\left( p^{n+1} - p^n\right)=0, \\
    \nabla \cdot u^c = -\nabla \cdot \tilde{u}^{n+1},

which is reducible to a Poisson problem:

.. math::
   \Delta p^{n+1} = \Delta p^n+\frac{\rho}{\Delta t}\nabla \cdot \tilde{u}^{n+1}.

The corrected velocity is then easily calculated from

.. math::
    u^{n+1} = \tilde{u}^{n+1}-\frac{\Delta t}{\rho}\nabla\left(p^{n+1}-p^n\right)

The scheme can be summarized in the following steps:
    #. Replace the pressure with a known approximation and solve for a tenative velocity :math:`\tilde{u}^{n+1}`.

    #. Solve a Poisson equation for the pressure, :math:`p^{n+1}`

    #. Use the corrected pressure to find the velocity correction and calculate :math:`u^{n+1}`

    #. Update t, and repeat.

.. [1] Goda, Katuhiko. *A multistep technique with implicit difference schemes for calculating two-or three-dimensional cavity flows.* Journal of Computational Physics 30.1 (1979): 76-95.

"""

from __future__ import division

from ufl import grad, Identity

from cbcflow.core.nsscheme import *

from cbcflow.schemes.utils import (
    # Time
    compute_regular_timesteps,
    # Spaces
    NSSpacePoolSplit,
    # ICs
    #assign_ics_split, # Modified implementation below
    # BCs
    make_velocity_bcs,
    make_pressure_bcs,
    make_penalty_pressure_bcs,
    )

def epsilon(u):
    "Symmetric gradient."
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, p, mu):
    "Stress tensor."
    return 2.0*mu*epsilon(u) - p*Identity(len(u))


def assign_ics_split(u, p, spaces, ics, annotate):
    """Assign initial conditions from ics to u0, p0.

    u0 is a vector valued function in spaces.V and p0 is a scalar function in spaces.Q,
    while ics = (icu, icp); icu = (icu0, icu1, ...).
    """
    #project(u, spaces.V, function=u, annotate=annotate) # TODO: Can do this in fenics dev
    #project(p, spaces.Q, function=p, annotate=annotate) # TODO: Can do this in fenics dev
    icu = project(as_vector(list(ics[0])), spaces.V, name="icu_projection")
    icp = project(ics[1], spaces.Q, name="icp_projection")
    u.assign(icu, annotate=annotate) # Force annotation
    p.assign(icp, annotate=annotate) # Force annotation


class IPCS(NSScheme):
    "Incremental pressure-correction scheme."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()

        # Set defalt Nietche parameters
        nietche = ParamDict(
            enable=True,
            formulation=1,
            stabilize=True,
            gamma=100.0,
            )

        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            #theta = 0.5,

            # Boundary condition method
            nietche=nietche,

            # Annotation on/off
            annotate=True,
            )
        return params

    def solve(self, problem, timer):
        # Get problem parameters
        mesh = problem.mesh
        n  = FacetNormal(mesh)
        dx = problem.dx
        ds = problem.ds

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        #t = Time(t0=timesteps[start_timestep])
        t = Constant(timesteps[start_timestep], name="TIME")

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
        u0 = Function(V, name="u0") # Previous timestep
        u1 = Function(V, name="u1") # Tentative / current timestep (could use one for tentative and one for current)
        p0 = Function(Q, name="p0") # Previous timestep
        p1 = Function(Q, name="p1") # Tentative / current timestep

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)
        cost_functionals = problem.cost_functionals(spaces, t, observations, controls)

        # TODO: Handle scaling of initial condition pressure?
        # Apply initial conditions
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_split(u1, p1, spaces, ics, annotate=self.params.annotate)

        # Get other problem specific functions
        bcs = problem.boundary_conditions(spaces, u1, p1, t, controls)
        body_force = problem.body_force(spaces, t)

        # Problem parameters
        nu = Constant(problem.params.mu/problem.params.rho, name="nu")
        rho = float(problem.params.rho)
        k  = Constant(dt, name="dt")
        f  = as_vector(body_force)

        # ======================================== FIXME Nietzche above FIXME ===================================

        # FIXME: Update BCs to allow Nietzche + strong
        # Make scheme-specific representation of bcs
        #bcu = make_velocity_bcs(problem, spaces, bcs)
        #bcp = make_pressure_bcs(problem, spaces, bcs)

        # Make scheme-specific representation of bcs
        a_bc, L_bc = 0, 0
        a_bc2, L_bc2 = 0, 0
        if self.params.nietche.enable:
            # Extract boundary conditions from bcs list (somewhat cumbersome now)
            bcu, bcp = bcs
            bcu_nietche = [bc[:2] for bc in bcu if bc[2] == "nietche"]
            bcu_strong = [bc[:2] for bc in bcu if bc[2] == "strong"]
            bcp_natural = [bc[:2] for bc in bcp]
            assert len(bcu) == len(bcu_strong) + len(bcu_nietche)

            # Define Nitche discretization constants
            gamma = Constant(self.params.nietche.gamma, name="gamma")
            hE = MinFacetEdgeLength(mesh)

            # Collect Nietche terms for each subboundary
            for g, region in bcu_nietche:
                g = as_vector(g)
                dsr = ds(region)

                # FIXME: Are these terms right? Maybe we need separate BC terms for u_tent and u_corr?
                # Add natural BC term
                a_bc += dot(p*n - nu*Dn(u), v)*dsr

                # Add Nietche terms
                # (all these formulations seem to work well _with_ stabilization but fail miserably without)
                s = self.params.nietche.formulation
                if s == 0:
                    # Don't touch B/B.T or A at all, stabilization terms only
                    pass
                elif s == 1:
                    # Symmetrize natural boundary term in B/B.T
                    a_bc += dot(u, q*n - nu*Dn(v))*dsr
                    L_bc += dot(g, q*n - nu*Dn(v))*dsr
                elif s == 2:
                    # Anti-symmetrize natural boundary term in B/B.T
                    a_bc += -dot(u, q*n - nu*Dn(v))*dsr
                    L_bc += -dot(g, q*n - nu*Dn(v))*dsr
                elif s == 3:
                    # Don't touch B/B.T but modify A (from Magnes experimental test)
                    a_bc += -dot(u, -nu*Dn(v))*dsr
                    L_bc += -dot(g, -nu*Dn(v))*dsr

                # Add Nietche stabilization terms
                if self.params.nietche.stabilize:
                    # Adding to separate a_bc2/L_bc2 to scale these terms by dt the same way as u_t,
                    # as doing so would weaken the terms vs the u_t terms for dt -> 0.
                    # We need to split a_bc/a_bc2 because the natural boundary terms above
                    # should be scaled as a natural part of the Navier-Stokes equation,
                    # and the symmetrization terms should be treated the same way.
                    a_bc2 += (gamma/hE)*dot(u, v)*dsr
                    L_bc2 += (gamma/hE)*dot(g, v)*dsr

        else:
            # Extract boundary conditions from bcs list (somewhat cumbersome now)
            bcu, bcp = bcs
            bcu_strong = [bc[:2] for bc in bcu]
            bcp_natural = [bc[:2] for bc in bcp]

        # ======================================== FIXME Nietzche above FIXME ===================================

        # Create Dirichlet boundary terms for velocity where it should be applied strongly
        # Note that this implies v = 0 and thus (p*n - nu*Dn(u)) . v = 0 on the same regions.
        bcu = [DirichletBC(V.sub(i), function, problem.facet_domains, region)
               for functions, region in bcu_strong
               for i, function in enumerate(functions)]

        # Add pressure boundary terms (should not overlap with velocity bc regions)
        # Note that this implies (p*n - nu*Dn(u)) = function*n, i.e. Dn(u) = 0.
        #L_bc += sum(-dot(function*n, v)*ds(region) for (function, region) in bcp_natural) # From coupled
        # Create DBCs for pressure from bcp_natural
        bcp = [DirichletBC(spaces.Q, function, problem.facet_domains, D)
               for function, region in bcp_natural]


        # FIXME: Recreate tent forms from F_tent, these are probably not correct. Also move boundary terms to bc code above.
        # Tentative velocity step
        a_u_tent = (
              (1/k) * inner(v, u) * dx
                    + 0.5 * inner(epsilon(v), sigma(u, p0, nu)) * dx
                    - 0.5 * nu * inner(grad(u).T*n, v) * ds
            )
        L_u_tent = (
              (1/k) * inner(v, u0) * dx
                    - inner(v, grad(u0)*u0) * dx
                    + 0.5 * inner(epsilon(v), sigma(u0, p0, nu)) * dx
                    + 0.5 * nu * inner(grad(u0).T*n, v) * ds
                    - inner(v, p0*n) * ds
                    + inner(v, f) * dx
            )

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx
        L_p_corr = inner(grad(q), grad(p0))*dx - (1/k)*q*div(u1)*dx

        # Velocity correction
        a_u_corr = inner(v, u)*dx
        L_u_corr = inner(v, u1)*dx - k*inner(v, grad(p1-p0))*dx

        # Add bc terms to both velocity equations FIXME: Figure out correct scaling with k or 1/k of these terms
        a_u_tent += a_bc + a_bc2
        L_u_tent += L_bc + L_bc2
        a_u_corr += a_bc + a_bc2
        L_u_corr += L_bc + L_bc2

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

        # Yield initial data for postprocessing
        state = (u1, p1)
        timestep = start_timestep
        data = ParamDict(timestep=timestep, t=float(t),
                         spaces=spaces,
                         u=u1, p=p1, # TODO: Remove u,p, only need NSSolver and other schemes to use state instead
                         state=state,
                         boundary_conditions=bcs,
                         controls=controls,
                         observations=observations, cost_functionals=cost_functionals)
        yield data

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            t0 = float(t)
            t.assign(timesteps[timestep])

            # Set up0 (u0, p0) to be u from previous timestep
            u0.assign(u1)
            p0.assign(p1)

            # Advance boundary conditions and body force from time t0 to time t
            problem.advance(t0, t, timestep, spaces, state,
                            bcs, body_force, controls)

            # Scale to solver pressure
            #p0.vector()[:] *= 1.0/rho # TODO: Handle otherwise

            # Compute tentative velocity step
            b = assemble(L_u_tent)
            for bc in bcu:
                bc.apply(A_u_tent, b)
            iter = solve(A_u_tent, u1.vector(), b, *self.params.solver_u_tent)

            # Compute pressure correction
            b = assemble(L_p_corr)
            if len(bcp) == 0:
                normalize(b)
            else:
                # Scale to physical pressure
                #b *= rho # TODO: Handle otherwise
                for bc in bcp:
                    bc.apply(A_p_corr, b)
                # ... and back to solver pressure
                #b *= 1.0/rho # TODO: Handle otherwise
            iter = solve(A_p_corr, p1.vector(), b, *solver_p_params)
            if len(bcp) == 0:
                normalize(p1.vector())

            # Compute velocity correction
            b = assemble(L_u_corr)
            for bc in bcu:
                bc.apply(A_u_corr, b)
            iter = solve(A_u_corr, u1.vector(), b, *self.params.solver_u_corr)

            # Scale to physical pressure
            #p0.vector()[:] *= rho # TODO: Handle otherwise

            # Callback to problem and yield data for postprocessing
            # (all variables here consistently time t)
            data = ParamDict(timestep=timestep, t=float(t),
                             spaces=spaces,
                             u=u1, p=p1, # TODO: Remove u,p, only need NSSolver and other schemes to use state instead
                             state=state,
                             boundary_conditions=bcs,
                             controls=controls,
                             observations=observations, cost_functionals=cost_functionals)
            yield data

        # Make sure annotation gets that the timeloop is over
        #finalize_time(t)
