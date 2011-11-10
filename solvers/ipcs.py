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

        if str(problem)=="Aneurysm":
            pc = "jacobi"
        else:
            pc = "ilu"

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        if self.options['segregated']:
            V = FunctionSpace(mesh, "CG", 1)
        else:
            V = VectorFunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)
        DG = FunctionSpace(mesh, "DG", 0)

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
        if self.options['segregated']:
            dims = range(len(u0));
            u0 = [interpolate(_u0, V) for _u0 in u0]
            u1 = [Function(V) for d in dims]
            u0_, u1_, bcu_ = u0, u1, bcu
        else:
            dims = [0]
            u0 = interpolate(u0, V)
            u1 = Function(V)
            u0_, u1_, bcu_ = [u0], [u1], [bcu]  # Always lists

        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

        # Tentative velocity step
        F_u_tent = []
        for d in dims:
            u_avg = 0.5 * (u + u0_[d])
            u_diff = (u - u0_[d])
            if self.options['segregated']:
                # FIXME, include also pressure term
                #            + inner(epsilon(v), sigma(Ux, p0, nu))*dx
                F_u_tent.append((1/k) * inner(v, u_diff) * dx
                                + v * sum(u0[r]*u0[d].dx(r) for r in dims) * dx
                                + inner(grad(v), grad(u_avg)) * dx
                                + inner(v, p0*n[d]) * ds
                                - v * f[d] * dx)
            else:
                F_u_tent.append((1/k) * inner(v, u_diff) * dx
                                + inner(v, grad(u0)*u0) * dx
                                + inner(epsilon(v), sigma(u_avg, p0, nu))*dx
                                + inner(v, p0*n)*ds
                                - beta*nu*inner(grad(u_avg).T*n, v)*ds
                                - inner(v, f)*dx
                                )
        a_u_tent = [lhs(F) for F in F_u_tent]
        L_u_tent = [rhs(F) for F in F_u_tent]

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx
        if self.options['segregated']:
            L_p_corr = inner(grad(q), grad(p0))*dx - (1/k)*q*sum(u1[r].dx(r) for r in dims)*dx
        else:
            L_p_corr = inner(grad(q), grad(p0))*dx - (1/k)*q*div(u1)*dx

        # Velocity correction
        a_u_corr = [inner(v, u)*dx for r in dims]
        if self.options['segregated']:
            L_u_corr = [v*u1[r]*dx - k*inner(v, grad(p1-p0)[r])*dx for r in dims]
        else:
            L_u_corr = [inner(v, u1)*dx - k*inner(v, grad(p1 - p0))*dx]

        # Assemble matrices
        A_u_tent = [assemble(a) for a in a_u_tent]
        A_p_corr = assemble(a_p_corr)
        A_u_corr = [assemble(a) for a in a_u_corr]

        # Time loop
        self.start_timing()
        for t in t_range:
            # Get boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)

            # Compute tentative velocity step
            for A, L, _u1, _bcu in zip(A_u_tent, L_u_tent, u1_, bcu_):
                b = assemble(L)
                for bc in _bcu: bc.apply(A, b)
                solve(A, _u1.vector(), b, "gmres", "ilu")
            self.timer("u1")

            # Pressure correction
            b = assemble(L_p_corr)
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            for bc in bcp: bc.apply(A_p_corr, b)
            if is_periodic(bcp):
                solve(A_p_corr, p1.vector(), b)
            else:
                solve(A_p_corr, p1.vector(), b, 'gmres', 'hypre_amg')
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            self.timer("p")

            # Velocity correction
            for A, L, _u1, _bcu in zip(A_u_corr, L_u_corr, u1_, bcu_):
                b = assemble(L)
                for bc in _bcu: bc.apply(A, b)
                solve(A, _u1.vector(), b, "gmres", pc)
            self.timer("u2")

            # Copy component-velocities into a vector field
            if self.options['segregated']:
                VV = VectorFunctionSpace(mesh, "Lagrange", 1)
                u1_v = Function(VV)
                for r in dims:
                    # FIXME: Depends on u1_v layout
                    n = len(u1[0].vector())
                    u1_v.vector()[r*n:r*n+n] = u1[r].vector()
            else:
                u1_v = u1
            self.timer("copy")

            # Update
            self.update(problem, t, u1_v, p1)
            for r in dims: u0_[r].assign(u1_[r])
            p0.assign(p1)
            self.timer("upd")

        return u1_v, p1

    def __str__(self):
        return "IPCS"
