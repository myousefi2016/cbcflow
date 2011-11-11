#from __future__ import division

__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-02-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.

from solverbase import *

class RHS(object):
    """The "instructions" to create b_u_tent: A list of (matrix,function)s and maybe a form."""
    def __init__(self):
        self.matvecs = []
        self.forms = []
        self.vecs = []
    def __iadd__(self, ins):
        if isinstance(ins, tuple):
            form, func = ins
            self.matvecs.append((assemble(form), func))
        elif isinstance(ins, GenericVector):
            self.vecs.append(ins)
        elif isinstance(ins, Form):
            self.forms.append(ins)
        else:
            raise RuntimeError, "Unknown RHS generator"
        return self
    def __call__(self):
        mat, func = self.matvecs[0]
        b = mat * func.vector()
        for mat, func in self.matvecs[1:]:
            b += mat * func.vector()
        for vec in self.vecs:
            b += vec
        for form in self.forms:
            assemble(form, tensor=b, add_values=True, reset_sparsity=False)
        return b

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
            V = FunctionSpace(mesh, "CG", 1)
        else:
            V = VectorFunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)
        #DG = FunctionSpace(mesh, "DG", 0)

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
        if self.options['segregated']:
            a_u_tent = []
            rhs_u_tent = []
            for d in dims:
                # FIXME, include also pressure term epsilon:sigma
                a_u_tent += [(1/k) * inner(v, u) * dx
                             # + inner(epsilon(v), nu*epsilon(u)) * dx
                             + nu * 0.5 * inner(grad(v), grad(u)) * dx]
                rhs = RHS()
                rhs += (-inner(v, p*n[d]) * ds
                         #- inner(epsilon(v), p*Indentity(u.cell().d)) * dx
                         , p0)
                rhs += (-nu * 0.5 * inner(grad(v), grad(u)) * dx
                         #- inner(epsilon(v), nu*epsilon(u)) * dx
                         + (1/k) * inner(v, u) * dx
                         , u0[d])
                rhs += v * f[d] * dx - v * sum(u0[r]*u0[d].dx(r) for r in dims) * dx
                rhs_u_tent.append(rhs)
        else:
            a_u_tent = ((1/k) * inner(v, u) * dx
                        + inner(epsilon(v), nu*epsilon(u)) * dx
                        - beta * nu * 0.5 * inner(grad(u).T*n, v) * ds)
            rhs_u_tent = RHS()
            rhs_u_tent += (-inner(v, p*n) * ds
                           + inner(epsilon(v), p*Identity(u.cell().d)) * dx
                            , p0)
            rhs_u_tent += (beta * nu * 0.5 * inner(grad(u).T*n, v) * ds
                           - inner(epsilon(v), nu*epsilon(u)) * dx
                           + (1/k) * inner(v, u) * dx
                           , u0)
            rhs_u_tent += inner(v, f) * dx + inner(v, grad(u0)*u0) * dx

            a_u_tent, rhs_u_tent = [a_u_tent], [rhs_u_tent]

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx
        if self.options['segregated']:
            rhs_p_corr = RHS()
            rhs_p_corr += inner(grad(q), grad(p)) * dx, p0
            for r in dims:
                rhs_p_corr += -(1/k) * q * u.dx(r) * dx, u1[r]
        else:
            rhs_p_corr = RHS()
            rhs_p_corr += inner(grad(q), grad(p))*dx, p0
            rhs_p_corr += -(1/k)*q*div(u)*dx, u1

        # Velocity correction
        a_u_corr = [inner(v, u)*dx for r in dims]
        if self.options['segregated']:
            rhs_u_corr = []
            for r in dims:
                rhs = RHS()
                rhs += v*u*dx, u1[r]
                rhs += -k*inner(v, grad(p)[r])*dx, p1
                rhs +=  k*inner(v, grad(p)[r])*dx, p0
                rhs_u_corr.append(rhs)
        else:
            rhs_u_corr = RHS()
            rhs_u_corr += inner(v, u)*dx, u1
            rhs_u_corr += -k*inner(v, grad(p))*dx, p1
            rhs_u_corr +=  k*inner(v, grad(p))*dx, p0
            rhs_u_corr = [rhs_u_corr]

        # Assemble matrices
        A_u_tent = [assemble(a) for a in a_u_tent]
        A_p_corr = assemble(a_p_corr)
        A_u_corr = [assemble(a) for a in a_u_corr]


        # Time loop
        self.start_timing()
        for t in t_range:
            # Get boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)
            self.timer("fetch bc")

            # Compute tentative velocity step
            for A, rhs, u1_comp, bcu_comp in zip(A_u_tent, rhs_u_tent, u1_, bcu_):
                b = rhs()
                self.timer("u1 assemble")
                for bc in bcu_comp: bc.apply(A, b)
                self.timer("u1 bc")
                iter = solve(A, u1_comp.vector(), b, "gmres", "jacobi")
                self.timer("u1 solve (%d, %d)"%(A.size(0), iter))

            # Pressure correction
            b = rhs_p_corr()
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            for bc in bcp: bc.apply(A_p_corr, b)
            self.timer("p assemble & bc")
            if is_periodic(bcp):
                iter = solve(A_p_corr, p1.vector(), b)
            else:
                iter = solve(A_p_corr, p1.vector(), b, 'gmres', 'ml_amg')
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            self.timer("p solve (%d, %d)"%(A_p_corr.size(0), iter))

            # Velocity correction
            for A, rhs, u1_comp, bcu_comp in zip(A_u_corr, rhs_u_corr, u1_, bcu_):
                b = rhs()
                for bc in bcu_comp: bc.apply(A, b)
                self.timer("u2 assemble & bc")
                iter = solve(A, u1_comp.vector(), b, "bicgstab", "ilu")
                self.timer("u2 solve (%d, %d)"%(A.size(0),iter))

            # Update
            self.update(problem, t, u1, p1)
            for r in dims: u0_[r].assign(u1_[r])
            p0.assign(p1)
            self.timer("update")

        return u1, p1

    def __str__(self):
        return "IPCS"
