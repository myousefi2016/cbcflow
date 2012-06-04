#from __future__ import division

__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-02-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.

from solverbase import *

class Solver(SolverBase):
    "Incremental pressure-correction scheme."

    def __init__(self, options):
        assert not options['segregated']
        SolverBase.__init__(self, options)

    def solve(self, problem):

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = self.select_timestep(problem)

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)

        # Remove boundary stress term if problem is periodic
        if is_periodic(bcp):
            beta = Constant(0)
        else:
            beta = Constant(1)

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u0 = interpolate(u0, V)
        u1 = Function(V)

        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f[0]
        n  = FacetNormal(mesh)

        # Tentative velocity step
        u_mean = 0.5 * (u + u0)
        u_diff = (u - u0)
        F_u_tent = ((1/k) * inner(v, u_diff) * dx
                    + inner(v, grad(u0)*u0) * dx
                    + inner(epsilon(v), sigma(u_mean, p0, nu)) * dx
                    - beta * nu * inner(grad(u_mean).T*n, v) * ds
                    + inner(v, p0*n) * ds
                    - inner(v, f) * dx)

        a_u_tent = lhs(F_u_tent)
        L_u_tent = rhs(F_u_tent)

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx
        L_p_corr = inner(grad(q), grad(p0))*dx - (1/k)*q*div(u1)*dx

        # Velocity correction
        a_u_corr = inner(v, u)*dx
        L_u_corr = inner(v, u1)*dx - k*inner(v, grad(p1-p0))*dx

        # Assemble matrices
        A_u_tent = assemble(a_u_tent)
        A_p_corr = assemble(a_p_corr)
        A_u_corr = assemble(a_u_corr)

        # Time loop
        self.start_timing()
        for t in t_range:
            # Get boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)
            self.timer("update & fetch bc")

            # Compute tentative velocity step
            b = assemble(L_u_tent)
            for bc in bcu: bc.apply(A_u_tent, b)
            self.timer("u1 construct rhs")
            solver = self.options['solver.u_tent']
            iter = solve(A_u_tent, u1.vector(), b, *solver)
            self.timer("u1 solve (%s, %d, %d)"%(', '.join(solver), b.size(), iter))

            # Pressure correction
            b = assemble(L_p_corr)
            if len(bcp) == 0 or is_periodic(bcp):
                solver = self.options['solver.p'] or self.options['solver.p_neumann']
                normalize(b)
            else:
                solver = self.options['solver.p'] or self.options['solver.p_dirichlet']
            for bc in bcp: bc.apply(A_p_corr, b)
            self.timer("p construct rhs")
            iter = solve(A_p_corr, p1.vector(), b, *solver)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            self.timer("p solve (%s, %d, %d)"%(', '.join(solver), b.size(), iter))

            # Velocity correction
            b = assemble(L_u_corr)
            for bc in bcu: bc.apply(A_u_corr, b)
            self.timer("u2 construct rhs")
            solver = self.options['solver.u_corr']
            iter = solve(A_u_corr, u1.vector(), b, *solver)
            self.timer("u2 solve (%s, %d, %d)"%(', '.join(solver), b.size(),iter))

            # Update
            self.update(problem, t, u1, p1)
            u0.assign(u1)
            p0.assign(p1)

        return u1, p1

    def __str__(self):
        return "IPCS"
