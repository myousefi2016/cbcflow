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
        SolverBase.__init__(self, options)

    def solve(self, problem):

        solver_u_tent      = "gmres", "jacobi"
        solver_p_periodic  = "cg", "ilu"
        solver_p_dirichlet = "gmres", "ml_amg"
        solver_u_corr      = "bicgstab", "ilu"

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        if self.options['segregated']:
            V = FunctionSpace(mesh, "CG", self.options['u_order'])
        else:
            V = VectorFunctionSpace(mesh, "CG", self.options['u_order'])
        Q = FunctionSpace(mesh, "CG", 1)

        # Get initial and boundary conditions
        ics = problem.initial_conditions(V, Q)
        bcs = problem.boundary_conditions(V, Q, t)
        u0, p0 = ics[:-1], ics[-1]
        bcu, bcp = bcs[:-1], bcs[-1]

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
        dims = range(len(u0))
        u0 = [interpolate(_u0, V) for _u0 in u0]
        u1 = [Function(V) for _u0 in u0]

        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

        # To avoid indexing in non-segregated forms
        if not self.options['segregated']:
            u0_ = u0[0]
            u1_ = u1[0]
            f_  = f[0]

        # Tentative velocity step
        if self.options['segregated']:
            F_u_tent = []
            for d in dims:
                u_mean = 0.5 * (u + u0[d])
                u_diff = (u - u0[d])
                F_u_tent += [(1/k) * inner(v, u_diff) * dx
                             + v * sum(u0[r]*u0[d].dx(r) for r in dims) * dx
                             + nu * inner(grad(v), grad(u_mean)) * dx
                             + inner(grad(v), 2*nu*grad(u_mean)) * dx - v.dx(d) * p0 * dx
                             + inner(v, p0*n[d]) * ds
                             - v * f[d] * dx]
        else:
            u_mean = 0.5 * (u + u0_)
            u_diff = (u - u0_)
            F_u_tent = [(1/k) * inner(v, u_diff) * dx
                        + inner(v, grad(u0_)*u0_) * dx
                        + inner(epsilon(v), sigma(u_mean, p0, nu)) * dx
                        - beta * nu * inner(grad(u_mean).T*n, v) * ds
                        + inner(v, p0*n) * ds
                        - inner(v, f_) * dx]

        a_u_tent = [lhs(F) for F in F_u_tent]
        L_u_tent = [rhs(F) for F in F_u_tent]

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx
        if self.options['segregated']:
            L_p_corr = inner(grad(q), grad(p0))*dx - (1/k)*q*sum(u1[r].dx(r) for r in dims)*dx
        else:
            L_p_corr = inner(grad(q), grad(p0))*dx - (1/k)*q*div(u1_)*dx

        # Velocity correction
        a_u_corr = [inner(v, u)*dx for r in dims]
        if self.options['segregated']:
            L_u_corr = [v*u1[r]*dx - k*inner(v, grad(p1-p0)[r])*dx for r in dims]
        else:
            L_u_corr = [inner(v, u1_)*dx - k*inner(v, grad(p1-p0))*dx]

        # Assemble matrices
        A_u_tent = [assemble(a) for a in a_u_tent]
        A_p_corr = assemble(a_p_corr)
        A_u_corr = [assemble(a) for a in a_u_corr]

        # Time loop
        self.start_timing()
        for t in t_range:
            # Get boundary conditions
            bcs = problem.boundary_conditions(V, Q, t)
            bcu, bcp = bcs[:-1], bcs[-1]
            self.timer("update & fetch bc")

            # Compute tentative velocity step
            for A, L, u1_comp, bcu_comp in zip(A_u_tent, L_u_tent, u1, bcu):
                b = assemble(L)
                for bc in bcu_comp: bc.apply(A, b)
                self.timer("u1 construct rhs")
                iter = solve(A, u1_comp.vector(), b, *solver_u_tent)
                self.timer("u1 solve (%s, %d, %d)"%(', '.join(solver_u_tent), A.size(0), iter))

            # Pressure correction
            b = assemble(L_p_corr)
            if len(bcp) == 0 or is_periodic(bcp):
                solver_p = solver_p_periodic
                normalize(b)
            else:
                solver_p = solver_p_dirichlet
            for bc in bcp: bc.apply(A_p_corr, b)
            self.timer("p construct rhs")
            iter = solve(A_p_corr, p1.vector(), b, *solver_p)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            self.timer("p solve (%s, %d, %d)"%(', '.join(solver_p), A_p_corr.size(0), iter))

            # Velocity correction
            for A, L, u1_comp, bcu_comp in zip(A_u_corr, L_u_corr, u1, bcu):
                b = assemble(L)
                for bc in bcu_comp: bc.apply(A, b)
                self.timer("u2 construct rhs")
                iter = solve(A, u1_comp.vector(), b, *solver_u_corr)
                self.timer("u2 solve (%s, %d, %d)"%(', '.join(solver_u_corr), A.size(0),iter))

            # Update
            self.update(problem, t, self._desegregate(u1), p1)
            for r in dims: u0[r].assign(u1[r])
            p0.assign(p1)

        return self._desegregate(u1), p1

    def __str__(self):
        return "IPCS"
