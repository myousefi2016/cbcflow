#from __future__ import division

__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-02-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.

from solverbase import *
from rhsgenerator import *

class Solver(SolverBase):
    "Incremental pressure-correction scheme."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

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
        dims = range(len(u0));
        u0 = [interpolate(_u0, V) for _u0 in u0]
        u1 = [Function(V) for d in dims]
        if not self.options['segregated']:
            # To avoid indexing in non-segregated forms
            u0_ = u0[0]
            u1_ = u1[0]

        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

        # Tentative velocity step
        M  = assemble(inner(v, u) * dx)
        K1 = assemble((1/k) * inner(v, u) * dx)
        if self.options['segregated']:
            K2 = assemble(inner(epsilon(v), nu*epsilon(u)) * dx
                          + 0.5 * nu * inner(grad(v), grad(u)) * dx)
            A_u_tent = []
            rhs_u_tent = []
            for d in dims:
                A = K1.copy()
                A += K2

                rhs = RhsGenerator(V)
                rhs += (-inner(v, p*n[d]) * ds
                         + v.dx(d) * p * dx
                         , p0)
                rhs += K1, u0[d]
                rhs -= K2, u0[d]
                rhs += M, f[d]
                rhs += -v * sum(u0[r]*u0[d].dx(r) for r in dims) * dx

                A_u_tent.append(A)
                rhs_u_tent.append(rhs)
        else:
            K2 = assemble(inner(epsilon(v), nu*epsilon(u)) * dx
                          - 0.5 * beta * nu * inner(v, grad(u).T*n) * ds)
            A = K1.copy()
            A += K2

            rhs = RhsGenerator(V)
            rhs += (-inner(v, p*n) * ds
                     + inner(epsilon(v), p*Identity(u.cell().d)) * dx
                     , p0)
            rhs += K1, u0_
            rhs -= K2, u0_
            rhs += M, f
            rhs += -inner(v, grad(u0_)*u0_) * dx

            A_u_tent, rhs_u_tent = [A], [rhs]

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx)
        if self.options['segregated']:
            rhs_p_corr = RhsGenerator(Q)
            rhs_p_corr += A_p_corr, p0
            for r in dims:
                rhs_p_corr += -(1/k) * q * u.dx(r) * dx, u1[r]
        else:
            rhs_p_corr = RhsGenerator(Q)
            rhs_p_corr += A_p_corr, p0
            rhs_p_corr += -(1/k)*q*div(u)*dx, u1_

        # Velocity correction
        A_u_corr = [M.copy() for r in dims]
        if self.options['segregated']:
            rhs_u_corr = []
            for r in dims:
                Kp = assemble(-k*inner(v, grad(p)[r])*dx)
                rhs = RhsGenerator(V)
                rhs += M, u1[r]
                rhs += Kp, p1
                rhs -= Kp, p0
                rhs_u_corr.append(rhs)
        else:
            Kp = assemble(-k*inner(v, grad(p))*dx)
            rhs_u_corr = RhsGenerator(V)
            rhs_u_corr += M, u1_
            rhs_u_corr += Kp, p1
            rhs_u_corr -= Kp, p0
            rhs_u_corr = [rhs_u_corr]

        # Apply BCs to matrices
        for A, bcs in zip(A_u_tent, bcu) + zip(A_u_corr, bcu):
            for bc in bcs:
                bc.apply(A)
        for bc in bcp:
            bc.apply(A_p_corr)

        # Time loop
        self.start_timing()
        for t in t_range:
            # Get boundary conditions
            bcs = problem.boundary_conditions(V, Q, t)
            bcu, bcp = bcs[:-1], bcs[-1]
            self.timer("update & fetch bc")

            # Compute tentative velocity step
            for A, rhs, u1_comp, bcu_comp in zip(A_u_tent, rhs_u_tent, u1, bcu):
                b = rhs()
                for bc in bcu_comp: bc.apply(b)
                self.timer("u1 assemble & bc")
                iter = solve(A, u1_comp.vector(), b, "gmres", "ilu")
                self.timer("u1 solve (%d, %d)"%(A.size(0), iter))

            # Pressure correction
            b = rhs_p_corr()
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            for bc in bcp: bc.apply(b)
            self.timer("p assemble & bc")
            if is_periodic(bcp):
                iter = solve(A_p_corr, p1.vector(), b, "cg", "ilu")
            else:
                iter = solve(A_p_corr, p1.vector(), b, 'gmres', 'ml_amg')
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            self.timer("p solve (%d, %d)"%(A_p_corr.size(0), iter))

            # Velocity correction
            for A, rhs, u1_comp, bcu_comp in zip(A_u_corr, rhs_u_corr, u1, bcu):
                b = rhs()
                for bc in bcu_comp: bc.apply(b)
                self.timer("u2 assemble & bc")
                iter = solve(A, u1_comp.vector(), b, "bicgstab", "jacobi")
                self.timer("u2 solve (%d, %d)"%(A.size(0),iter))

            # Update
            self.update(problem, t, self._desegregate(u1), p1)
            for r in dims: u0[r].assign(u1[r])
            p0.assign(p1)

        return self._desegregate(u1), p1

    def __str__(self):
        return "IPCS"
