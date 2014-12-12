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

from argparse import Namespace

from cbcflow.core.nsscheme import *

from cbcflow.schemes.utils import (compute_regular_timesteps,
#                                   assign_ics_segregated,
#                                   make_segregated_velocity_bcs,
#                                   make_pressure_bcs,
                                   NSSpacePoolSegregated)

def assign_ics_mixed(up, spaces, ics, annotate):
    """Assign initial conditions from ics to up.

    up is a mixed function in spaces.W = spaces.V * spaces.Q,
    while ics = (icu, icp); icu = (icu0, icu1, ...).
    """
    #project(as_vector(list(ics[0]) + [ics[1]]), spaces.W, function=up, annotate=annotate) # TODO: Can pass function in fenics-dev
    icup = project(as_vector(list(ics[0]) + [ics[1]]), spaces.W, name="icup_projection")
    up.assign(icup, annotate=annotate) # Force annotation

def SegregatedFunction(U, dim, name):
    return as_vector([Function(U, name="{0}[comp={1}]".format(name, i))
                      for i in range(dim)])

def SegregatedTimeFunction(U, dim, name, timesteps):
    return { ts: SegregatedFunction(U, dim, "{0}[ts={1}]".format(name, ts))
             for ts in timesteps }

def TimeFunction(Q, name, timesteps):
    return { ts: Function(Q, name="{0}[ts={1}]".format(name, ts))
             for ts in timesteps }

class IPCS_Stable_DA(NSScheme):
    "."

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

            # Annotation on/off
            annotate=True,

            # Solver params
            u_tent_prec_structure = "same_nonzero_pattern",
            u_tent_solver_parameters = {},
            p_corr_solver_parameters = {},
            u_corr_solver_parameters = {},

            # Form compiler parameters for optimizing run time
            form_compiler_parameters=fc,
            )
        return params

    def solve(self, problem, timer):
        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        dim = mesh.topology().dim()
        dims = range(dim)
        theta = self.params.theta

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Constant(timesteps[start_timestep], name="TIME")

        # Define function spaces
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)
        U = spaces.U
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(U)
        q = TestFunction(Q)
        ut = TrialFunction(U)
        pt = TrialFunction(Q)

        # Input functions with history
        f = SegregatedTimeFunction(U, dim, "f", (0,+1))
        g = SegregatedTimeFunction(U, dim, "g", (0,+1))

        # Solution functions with or without history
        u = SegregatedTimeFunction(U, dim, "u", (-1,0,+1))
        utent = SegregatedFunction(U, dim, "utent")
        p = TimeFunction(Q, "p", (0, 1))

        # Extrapolate from solutions at t[-1], t[0] to t=0.5*(t[0]+t[1])
        u_ext = as_vector([1.5*u[0][i] + 0.5*u[-1][i] for i in dims])
        #p_ext = 1.5*p[0] + 0.5*p[-1]
        p_ext = p[0]

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)
        cost_functionals = problem.cost_functionals(spaces, t, observations, controls)

        # Apply initial conditions and use it as initial guess
        ics = problem.initial_conditions(spaces, controls)
        # Set u^+1 to be u at t0
        assign_ics_segregated(u[1], p[1], spaces, ics) # FIXME
        # Set u^0 to be u at t0-dt
        for i in dims:
            u[0][i].assign(u[1][i])
        p0.assign(p1)
        # Let u^-1 be undefined, overwritten first thing in time loop.

        # Problem coefficients
        rho = Constant(problem.params.rho, name="rho")
        nu = Constant(float(problem.params.mu) / float(problem.params.rho), name="nu")
        k = Constant(dt, name="dt")
        body_force = problem.body_force(spaces, t)
        f  = as_vector(body_force)

        # Setting alpha = 0 gives mass matrix in velocity correction (dropping stiffness term)
        #alpha = 1
        alpha = 0

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u[1], p[1], t, controls)
        bcu = make_segregated_velocity_bcs(problem, spaces, bcs) # FIXME
        bcp = make_pressure_bcs(problem, spaces, bcs) # FIXME
        # FIXME: Setup Nietzche form terms and strong bcs here
        # Setup strong bcs
        bcs_u_tent = []
        bcs_p = []
        bcs_u_corr = []

        # Tentative velocity equation
        # Convection linearized as in Simo/Armero (1994) (see u_ext)
        a_u_tent = [(
            (1.0/k) * inner(ut, v)*dx
            + theta * inner(dot(grad(ut), u_ext), v)*dx
            + theta*nu * inner(grad(ut), grad(v))*dx # (iii)
            - theta*nu * inner(dot(grad(ut), n), v)*ds # (iii)
            ) for i in dims]
        b_u_tent = [(
            (1.0/k)*inner(u[0][i], v)*dx
            - (1-theta) * inner(dot(grad(u[0][i]), u_ext), v)*dx
            - (1-theta)*nu * inner(grad(u[0][i]), grad(v))*dx # (ii)
            + (1.0/rho) * p_ext*v.dx(i)*dx # (i)
            + inner(f[0][i], v)*dx
            # Boundary terms from partial integration:
            - (1.0/rho) * inner(p_ext*n[i], v)*ds # (i)
            - (1-theta)*nu * inner(dot(grad(u[0][i]), n), v)*ds # (ii)
            ) for i in dims]

        # Pressure equation
        a_p = (
            (1.0/rho) * inner(grad(pt), grad(q))*dx # (i)
            # Boundary terms from partial integration:
            - (1.0/rho) * inner(dot(grad(pt), n), q)*ds # (i)
            )
        b_p = (
            (1.0/rho) * inner(grad(p_ext), grad(q))*dx # (ii)
            #+ alpha*(theta*nu) * inner(grad(div(utent)), grad(q))*dx # (iii) for u degree 2, grad(div(utent)) is nonzero, can we drop it then?
            - (1.0/k)*div(utent)*q*dx
            #+ (1.0/k)*inner(utent, grad(q))*dx (iv)
            #- (1.0/k)*inner(dot(utent, n), q)*ds (iv)
            # Boundary terms from partial integration:
            - (1.0/rho) * inner(dot(grad(p_ext), n), q)*ds # (ii)
            #- alpha*(theta*nu) * inner(dot(grad(div(utent)), n), q)*ds # (iii) for u degree 2, grad(div(utent)) is nonzero, can we drop it then?
            )

        # Velocity correction equation
        a_u_corr = [(
            inner(ut, v)*dx
            + alpha*(k*theta*nu) * inner(grad(ut), grad(v))*dx # (i)
            # Boundary terms from partial integration:
            - alpha*(k*theta*nu) * inner(dot(grad(ut), n), v)*ds # (i)
            ) for i in dims]
        b_u_corr = [(
            inner(u[0][i], v)*dx
            + alpha*(k*theta*nu) * inner(grad(utent[i]), grad(v))*dx # (ii)
            + k/rho * inner(p[+1].dx(i), v)*dx
            - k/rho * inner(p_ext.dx(i), v)*dx
            # Boundary terms from partial integration:
            + alpha*(k*theta*nu) * inner(dot(grad(utent[i]), n), v)*ds # (ii)
            ) for i in dims]

        # FIXME: Setup solvers here

        # Preassemble
        lhs_p = assemble(a_p)
        rhs_p = assemble(b_p)
        for bc in bcs_p:
            bc.apply(lhs_p)

        lhs_u_corr = [assemble(a_u_corr[i]) for i in dims]
        rhs_u_corr = [assemble(b_u_corr[i]) for i in dims]
        for i in dims:
            for bc in bcs_u_corr:
                bc.apply(lhs_u_corr[i])

        # Yield initial data for postprocessing
        state = (u[1], p[1])
        yield Namespace(timestep=start_timestep, t=float(t), spaces=spaces, state=state,
                        boundary_conditions=bcs,
                        controls=controls,
                        observations=observations,
                        cost_functionals=cost_functionals)

        # Time loop
        for timestep in xrange(start_timestep+1,len(timesteps)):
            # Set t0,t to the beginning and end of the
            # interval we're computing over now
            t0 = float(t)
            t.assign(timesteps[timestep])

            # Advance boundary conditions and body force from time t0 to time t
            problem.advance(t0, t, timestep, spaces, state,
                            bcs, body_force, controls) # FIXME: bcs, body_force?

            # Set u^-1 to be u at t0-dt
            for i in dims:
                u[-1][i].assign(u[0][i])

            # Set u^0 to be u at t0
            for i in dims:
                u[0][i].assign(u[1][i])

            # Find tentative u
            lhs_u_tent = [assemble(a_u_tent[i]) for i in dims]
            rhs_u_tent = [assemble(b_u_tent[i]) for i in dims]
            for i in dims:
                for bc in bcs_u_tent:
                    bc.apply(rhs_u_tent[i])
            for i in dims:
                solve(lhs_u_tent[i], rhs_u_tent[i], utent[i].vector())
            # TODO: Setup solvers outside of loop
            #for i in dims:
            #    solver_u_tent[i].solve()

            # Find p^1 (pressure at time t)
            for bc in bcs_p:
                bc.apply(rhs_p)
            solve(lhs_p, rhs_p, p[1].vector())
            # TODO: Setup solvers outside of loop
            #solver_p.solve()

            # Find u^1 (final velocity at time t)
            for i in dims:
                for bc in bcs_u_corr:
                    bc.apply(rhs_u_corr[i])
            for i in dims:
                solve(lhs_u_corr[i], rhs_u_corr[i], u[1][i].vector())
            # TODO: Setup solvers outside of loop
            #for i in dims:
            #    solver_u_corr[i].solve()

            # Yield computed data for postprocessing
            yield Namespace(timestep=timestep, t=float(t), spaces=spaces, state=state,
                            boundary_conditions=bcs,
                            controls=controls,
                            observations=observations,
                            cost_functionals=cost_functionals)
