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

from cbcflow.schemes.utils import (
    # Time
    compute_regular_timesteps,
    # Spaces
    NSSpacePoolMixed,
    # ICs
    assign_ics_mixed,
    # BCs
    # ... implemented inline below for now
    )

def assign_ics_mixed(up, spaces, ics):
    """Assign initial conditions from ics to up.

    up is a mixed function in spaces.W = spaces.V * spaces.Q,
    while ics = (icu, icp); icu = (icu0, icu1, ...).
    """
    icup = as_vector(list(ics[0]) + [ics[1]])
    #project(icup, spaces.W, function=up) # TODO: Can do this in fenics dev
    up.assign(project(icup, spaces.W, name="icup_projection"))

class CoupledScheme(NSScheme):
    "Coupled scheme using a fixed point (Picard) nonlinear solver."

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

        # Set default nonlinear solver params, making all dolfin solver params available
        nonlinear_solver = ParamDict(NonlinearVariationalSolver.default_parameters().to_dict())
        nonlinear_solver.newton_solver.absolute_tolerance = 1e-13
        nonlinear_solver.newton_solver.relative_tolerance = 1e-13
        nonlinear_solver.newton_solver.maximum_iterations = 10
        nonlinear_solver.newton_solver.error_on_nonconvergence = False
        nonlinear_solver.newton_solver.report = True
        nonlinear_solver.newton_solver.linear_solver = "lu"
        nonlinear_solver.reset_jacobian = True

        # Set default form compiler parameters
        fc = ParamDict(
            optimize=True,
            cpp_optimize=True,
            cpp_optimize_flags="-O2",
            #cpp_optimize_flags="-O3 -march=native -fno-math-errno",
            quadrature_degree="auto",
            )

        params.update(
            # Default to P2-P1 (Taylor-Hood)
            u_degree=2,
            p_degree=1,
            theta=1.0,

            # Choice of equation formulation
            scale_by_dt=True,
            enable_convection=True, # False for Stokes

            # Boundary condition method
            nietche=nietche,

            # Nonlinear solver parameters
            picard_newton_fraction=1.0, # 0.0 = Picard, 1.0 = Newton, (0.0,1.0) = mix
            nonlinear_solver=nonlinear_solver,

            # Form compiler parameters for optimizing run time
            form_compiler_parameters=fc,
            )
        return params

    def solve(self, problem, timer):
        # Spatial parameters
        mesh = problem.mesh
        n  = FacetNormal(mesh)
        dx = problem.dx
        ds = problem.ds

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)

        #t = Time(t0=timesteps[start_timestep])
        t = Constant(timesteps[start_timestep], name="TIME")

        # Function spaces
        spaces = NSSpacePoolMixed(mesh, self.params.u_degree, self.params.p_degree)
        W = spaces.W

        # Test and trial functions
        up = TrialFunction(W)
        vq = TestFunction(W)
        u, p = split(up)
        v, q = split(vq)

        # Solution functions
        up0 = Function(W, name="up0") # Previous timestep
        up1 = Function(W, name="up1") # Last iterate of fixed-point in current timestep
        u0, p0 = split(up0)
        u1, p1 = split(up1)

        # Get problem specific functions
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)
        cost_functionals = problem.cost_functionals(spaces, t, observations, controls)
        ics = problem.initial_conditions(spaces, controls)
        bcs = problem.boundary_conditions(spaces, u1, p1, t, controls)
        body_force = problem.body_force(spaces, t)

        # Problem parameters
        nu = Constant(problem.params.mu/problem.params.rho, name="nu")
        k  = Constant(dt, name="dt")
        f  = as_vector(body_force)

        if self.params.scale_by_dt:
            # Scaled by k (smaller residual, nonlinear solver hits absolute stopping criteria faster)
            kinv = 1
            kval = k
        else:
            # Not scaled by k (keep u_t = (u1-u0)/dt, larger residual, nonlinear solver may not converge properly)
            kinv = 1.0 / k
            kval = 1

        # Apply initial conditions and use it as initial guess
        assign_ics_mixed(up1, spaces, ics)

        # Make scheme-specific representation of bcs
        abc, Lbc = 0, 0
        abc2, Lbc2 = 0, 0
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

                # Add natural BC term
                abc += dot(p*n - nu*Dn(u), v)*dsr

                # Add Nietche terms
                # (all these formulations seem to work well _with_ stabilization but fail miserably without)
                s = self.params.nietche.formulation
                if s == 0:
                    # Don't touch B/B.T or A at all, stabilization terms only
                    pass
                elif s == 1:
                    # Symmetrize natural boundary term in B/B.T
                    abc += dot(u, q*n - nu*Dn(v))*dsr
                    Lbc += dot(g, q*n - nu*Dn(v))*dsr
                elif s == 2:
                    # Anti-symmetrize natural boundary term in B/B.T
                    abc += -dot(u, q*n - nu*Dn(v))*dsr
                    Lbc += -dot(g, q*n - nu*Dn(v))*dsr
                elif s == 3:
                    # Don't touch B/B.T but modify A (from Magnes experimental test)
                    abc += -dot(u, -nu*Dn(v))*dsr
                    Lbc += -dot(g, -nu*Dn(v))*dsr

                # Add Nietche stabilization terms
                if self.params.nietche.stabilize:
                    # Adding to separate abc2/Lbc2 to scale these terms by dt the same way as u_t,
                    # as doing so would weaken the terms vs the u_t terms for dt -> 0.
                    # We need to split abc/abc2 because the natural boundary terms above
                    # should be scaled as a natural part of the Navier-Stokes equation,
                    # and the symmetrization terms should be treated the same way.
                    abc2 += (gamma/hE)*dot(u, v)*dsr
                    Lbc2 += (gamma/hE)*dot(g, v)*dsr

        else:
            # Extract boundary conditions from bcs list (somewhat cumbersome now)
            bcu, bcp = bcs
            bcu_strong = [bc[:2] for bc in bcu]
            bcp_natural = [bc[:2] for bc in bcp]

        # Create Dirichlet boundary terms for velocity where it should be applied strongly
        # Note that this implies v = 0 and thus (p*n - nu*Dn(u)) . v = 0 on the same regions.
        bc_strong = [DirichletBC(spaces.Ubc[i], function, problem.facet_domains, region)
                     for functions, region in bcu_strong
                     for i, function in enumerate(functions)]

        # Add pressure boundary terms (should not overlap with velocity bc regions)
        # Note that this implies (p*n - nu*Dn(u)) = function*n, i.e. Dn(u) = 0.
        Lbc += sum(-dot(function*n, v)*ds(region) for (function, region) in bcp_natural)

        # Set up Picard linearization of (Navier-)Stokes, F = a*u - L = 0
        # Stokes terms
        a = kinv*dot(u, v)*dx + kval*(nu*inner(grad(u), grad(v)) - p*div(v) - q*div(u))*dx
        L = kinv*dot(u0, v)*dx + kval*dot(f, v)*dx

        # Add convection
        if self.params.enable_convection:
            # Note the convection linearization u1 . nabla u . v, which becomes u1 . nabla u1 . v in F below
            a += kval*dot(grad(u)*u1, v)*dx

        # BC terms scaled by dt or 1/dt
        a += kval*abc + kinv*abc2
        L += kval*Lbc + kinv*Lbc2

        # Build residual form from a, L
        F = action(a, up1) - L

        # Jacobian for Newton method
        J = derivative(F, up1, up)

        # Select between Picard, Newton or a mix with alpha in [0,1]
        alpha = float(self.params.picard_newton_fraction)
        if alpha == 0.0:
            # Use Picard linearization only
            J_approx = a
        elif alpha == 1.0:
            # Use full Newton linearization
            J_approx = J
        else:
            # Use a mix (avoid recompilation with Constant)
            alpha = Constant(alpha, name="picard_newton_fraction")
            J_approx = alpha*J + (1.0-alpha)*a

        # Create solver
        nonlinear_problem = NonlinearVariationalProblem(F, up1, bc_strong, J=J_approx,
                                                        form_compiler_parameters=self.params.form_compiler_parameters)
        solver = NonlinearVariationalSolver(nonlinear_problem)
        solver.parameters.update(self.params.nonlinear_solver)

        # TODO: Could use this to yield u0,p0 in proper spaces below:
        #u_assigner = FunctionAssigner(spaces.V, spaces.W.sub(0))
        #p_assigner = FunctionAssigner(spaces.Q, spaces.W.sub(1))
        #uv = Function(spaces.V)
        #pv = Function(spaces.Q)
        #u_assigner.assign(uv, up0.sub(0))
        #p_assigner.assign(pv, up0.sub(1))

        # Callback to problem and yield data for postprocessing
        # (all variables here are consistently at time t)
        state = (u1, p1)
        timestep = start_timestep
        data = ParamDict(timestep=timestep, t=float(t),
                         spaces=spaces,
                         u=u1, p=p1, bcs=bcs, # TODO: Remove
                         state=state,
                         boundary_conditions=bcs,
                         controls=controls,
                         observations=observations, cost_functionals=cost_functionals)
        #problem.post_update(**data)
        yield data

        # Loop over fixed timesteps
        ##adj_start_timestep(float(t))
        for timestep in xrange(start_timestep+1, len(timesteps)):
            t0 = float(t)
            #assign_time(t, timesteps[timestep])
            t.assign(timesteps[timestep])

            # Set up0 (u0, p0) to be u from previous timestep
            up0.assign(up1)

            # Advance boundary conditions and body force from time t0 to time t
            problem.advance(t0, t, timestep, spaces, state,
                            bcs, body_force, controls)

            # Solve for up1 (u1, p1 becomes u, p at time t)
            solver.solve()
            #solve(F == 0, up1, bc_strong, J=J_approx,
            #      form_compiler_parameters=self.params.form_compiler_parameters,
            #      solver_parameters=self.params.nonlinear_solver.to_dict())

            # Callback to problem and yield data for postprocessing
            # (all variables here consistently time t)
            data = ParamDict(timestep=timestep, t=float(t),
                             spaces=spaces,
                             u=u1, p=p1, bcs=bcs, # TODO: Remove
                             state=state,
                             boundary_conditions=bcs,
                             controls=controls,
                             observations=observations, cost_functionals=cost_functionals)
            #problem.post_update(**data)
            yield data

            ##adj_inc_timestep(float(t), finished=timestep == len(timesteps)-1)

        # Make sure annotation gets that the timeloop is over
        #finalize_time(t)
