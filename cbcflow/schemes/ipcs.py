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


class IPCS(NSScheme):
    "Incremental pressure-correction scheme."

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
        t = Time(t0=timesteps[start_timestep])
        
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

        # Get initial conditions
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_split(u0, p0, spaces, ics)
        u1.assign(u0)
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        # Tentative velocity step
        u_mean = 0.5 * (u + u0)
        u_diff = (u - u0)
        F_u_tent = ((1/k) * inner(v, u_diff) * dx()
                    + inner(v, grad(u0)*u0) * dx()
                    + inner(epsilon(v), sigma(u_mean, p0, nu)) * dx()
                    - beta * nu * inner(grad(u_mean).T*n, v) * ds()
                    + inner(v, p0*n) * ds()
                    - inner(v, f) * dx())

        a_u_tent = lhs(F_u_tent)
        L_u_tent = rhs(F_u_tent)

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx()
        L_p_corr = inner(grad(q), grad(p0))*dx() - (1/k)*q*div(u1)*dx()

        # Velocity correction
        a_u_corr = inner(v, u)*dx()
        L_u_corr = inner(v, u1)*dx() - k*inner(v, grad(p1-p0))*dx()

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
            b = assemble(L_u_tent)
            for bc in bcu: bc.apply(A_u_tent, b)
            timer.completed("u1 construct rhs")

            iter = solve(A_u_tent, u1.vector(), b, *self.params.solver_u_tent)
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
            b = assemble(L_u_corr)
            for bc in bcu: bc.apply(A_u_corr, b)
            timer.completed("u2 construct rhs")

            solver_params = self.params.solver_u_corr
            iter = solve(A_u_corr, u1.vector(), b, *solver_params)
            timer.completed("u2 solve (%s, %d, %d)"%(', '.join(solver_params), b.size(),iter))

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