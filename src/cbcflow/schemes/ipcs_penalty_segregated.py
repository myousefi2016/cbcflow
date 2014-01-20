from __future__ import division

__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2012-02-15"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Martin Alnaes, 2013. (added penalty bcs)

from ..core.nsscheme import *
from ..core.utils import Timer
from ..core.timesteps import compute_regular_timesteps
from ..core.schemeutils import (assign_ics_segregated,
                                make_segregated_velocity_bcs,
                                make_pressure_bcs,
                                make_penalty_pressure_bcs)
from ..core.spaces import NSSpacePoolSegregated


class SegregatedPenaltyIPCS(NSScheme):
    "Segregated incremental pressure-correction scheme with penalty terms for pressure BCs."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,

            # Penalty pressure bcs is necessary for optimization
            use_penalty_pressure_bcs=True,
            penalty_gamma=100.0,
            )
        return params

    def solve(self, problem, update, restart=None):
        # Spatial parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        dims = range(mesh.topology().dim())

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

        # Function spaces
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)
        U = spaces.U
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
        if self.params.use_penalty_pressure_bcs:
            bcp = []
            a_pbc, L_pbc = make_penalty_pressure_bcs(problem, spaces, bcs,
                                                     self.params.penalty_gamma, q, p)
        else:
            bcp = make_pressure_bcs(problem, spaces, bcs)
            a_pbc, L_pbc = 0, 0

        # Problem parameters
        nu = Constant(problem.params.mu/problem.params.rho, name="nu")
        k  = Constant(dt, name="dt")
        f  = as_vector(problem.body_force(spaces, t))

        # Tentative velocity forms
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

        # Pressure correction forms
        a_p_corr = inner(grad(q), grad(p))*dx() + a_pbc
        L_p_corr = inner(grad(q), grad(p0))*dx() - (1/k)*q*div(u1)*dx() + L_pbc

        # Velocity correction forms
        a_u_corr = [inner(v, u)*dx() for r in dims]
        L_u_corr = [v*u1[r]*dx() - k*inner(v, grad(p1-p0)[r])*dx() for r in dims]

        # Assemble matrices
        A_u_tent = [assemble(a) for a in a_u_tent]
        A_p_corr = assemble(a_p_corr)
        A_u_corr = [assemble(a) for a in a_u_corr]

        # TODO: There's still a bug with LinearSolver in dolfin-adjoint
        # TODO: There's still a bug with LUSolver in dolfin-adjoint
        # Create solvers
        if self.params.solver_u_tent[0] == "lu":
            solver_u_tent = LUSolver()
        else:
            solver_u_tent = KrylovSolver(*self.params.solver_u_tent)

        if self.params.solver_u_corr[0] == "lu":
            solver_u_corr = LUSolver()
        else:
            solver_u_corr = KrylovSolver(*self.params.solver_u_corr)

        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0:
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        p_corr_problem = LinearVariationalProblem(a_p_corr, L_p_corr, p1, bcp)
        solver_p_corr = LinearVariationalSolver(p_corr_problem)
        (solver_p_corr.parameters["linear_solver"],
         solver_p_corr.parameters["preconditioner"]) = solver_p_params
        solver_p_corr.parameters["symmetric"] = True
        #solver_p_corr.solve() # Adding this line (which is a bug) causes the reported dolfin-adjoint warning

        update(u0, p0, float(t), start_timestep, spaces)

        # Profiling object
        timer = Timer(self.params.enable_timer)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Compute tentative velocity step
            for d in dims:
                b = assemble(L_u_tent[d])
                for bc in bcu: bc[d].apply(A_u_tent[d], b)
                timer.completed("u_tent construct rhs")

                iter = solver_u_tent.solve(A_u_tent[d], u1[d].vector(), b)
                timer.completed("u_tent solve (%s, %d dofs, %d iter)" % (
                    ', '.join(self.params.solver_u_tent), b.size(), iter))

            # Pressure correction
            solver_p_corr.solve()
            iter = -1
            timer.completed("p_corr solve (%s, %d dofs)" % (
                ', '.join(solver_p_params), b.size()))

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
        scope = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return scope
