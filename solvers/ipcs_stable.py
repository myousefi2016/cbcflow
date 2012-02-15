from __future__ import division

__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2012-02-15"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *
from rhsgenerator import *

class Solver(SolverBase):
    "Incremental pressure-correction scheme."

    def __init__(self, options):
        assert options['segregated']
        SolverBase.__init__(self, options)

    def solve(self, problem):

        solver_u_tent_params      = "gmres", "hypre_euclid"
        solver_p_periodic_params  = "gmres", "hypre_euclid"
        solver_p_dirichlet_params = "gmres", "ml_amg"
        solver_u_corr_params      = "bicgstab", "hypre_euclid"

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        V = FunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)

        # Get initial and boundary conditions
        ics = problem.initial_conditions(V, Q)
        bcs = problem.boundary_conditions(V, Q, t)
        u0, p0 = ics[:-1], ics[-1]
        bcu, bcp = bcs[:-1], bcs[-1]

        # Interpolate initial conditions to functions
        u_curr = [interpolate(_, V) for _ in u0]
        p_curr = interpolate(p0, Q)

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Problem parameters etc
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

        dim  = len(u1)
        dims = range(dim)

        # Functions

        u_prev = [_.copy() for _ in u0] # u^{n-1}
        u_curr = [_.copy() for _ in u0] # u^{n}
        u_next = [_.copy() for _ in u0] # u^{n+1}
        p_curr = p0.copy()              # p^{n}
        p_next = p0.copy()              # p^{n+1}

        # Tentative velocity step
        a1 = (1/k) * inner(v, u) * dx
        a2 = 0.5 * inner(grad(v), nu*grad(u)) * dx

        a_conv = v * sum((0.75*u_curr[r] - 0.25*u_prev[r]) * u.dx(r) for r in dims) * dx
        Kconv = Matrix() # assembled from a_conv in the time loop

        # Create the static part of the coefficient matrix for the tentative
        # velocity step. The convection is added in the time loop. We use a
        # single matrix for all dimensions, which means that the BCs must match
        # (they must apply to the full vector space at the vertex, not a
        # subspace.
        A_u_tent = assemble(a1+a2)

        # Create matrices for generating the RHS
        B = assemble(a1-a2)
        M = assemble(v*u*dx)

        # Define how to create the RHS for the tentative velocity. The RHS is
        # (of course) different for each dimension.
        rhs_u_tent = [None]*dim
        for d in dims:
            C = assemble(-v*p*n[d]*ds + v.dx(d)*p*dx)
            rhs_u_tent[d] = RhsGenerator(V)
            rhs_u_tent[d] += B, u_curr[d]
            rhs_u_tent[d] += C, p_curr
            rhs_u_tent[d] += M, f[d]
            rhs_u_tent[d] -= Kconv, u_curr[d]

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx)
        rhs_p_corr = RhsGenerator(Q)
        rhs_p_corr += A_p_corr, p_curr
        for d in dims:
            Ku = assemble(-(1/k)*q*u.dx(d)*dx)
            rhs_p_corr += Ku, u_next[d]

        # Velocity correction. Like for the tentative velocity, a single LHS is used.
        A_u_corr = M
        rhs_u_corr = [None]*dim
        for d in dims:
            Kp = assemble(-k*inner(v, grad(p)[d])*dx)
            rhs_u_corr[d] = RhsGenerator(V)
            rhs_u_corr[d] += M, u_next[d]
            rhs_u_corr[d] += Kp, p_next
            rhs_u_corr[d] -= Kp, p_curr

        # Apply BCs to LHS
        for bc in bcu[0]:
            bc.apply(A_u_tent)
            bc.apply(A_u_corr)
        for bc in bcp:
            bc.apply(A_p_corr)

        # Create solvers
        if is_periodic(bcp): solver_p_params = solver_p_periodic_params
        else:                solver_p_params = solver_p_dirichlet_params
        solver_u_tent = LinearSolver(*solver_u_tent_params)
        solver_p_corr = LinearSolver(*solver_p_params)
        solver_u_corr = LinearSolver(*solver_u_corr_params)

        for A,S in [(A_u_tent, solver_u_tent), (A_p_corr, solver_p_corr), (A_u_corr, solver_u_corr)]:
            S.set_operator(A)
            S.parameters['preconditioner']['reuse'] = True

        # Time loop
        self.start_timing()
        for t in t_range:
            # Get boundary conditions
            bcs = problem.boundary_conditions(V, Q, t)
            bcu, bcp = bcs[:-1], bcs[-1]
            self.timer("update & fetch bc")

            # Assemble the u-dependent convection matrix. It is important that
            # it is assembled into the same tensor, because the tensor is
            # stored in rhs. (And it's faster).
            if Kconv.size(0) > 0:
                A_u_tent.axpy(-1.0, Kconv, True)
                assemble(a_conv, tensor=Kconv, reset_sparsity=False)
            else:
                assemble(a_conv, tensor=Kconv, reset_sparsity=True)
            A_u_tent.axpy(1.0, Kconv, True)
            self.timer("u_tent assemble convection & construct lhs")

            # Compute tentative velocity step
            for d in dims:
                b = rhs_u_tent[d](bcs=bcu[d])
                self.timer("u_tent construct rhs")
                iter = solver_u_tent.solve(u_next[d].vector(), b)
                self.timer("u_tent solve (%s, %d dofs, %d iter)"%(', '.join(solver_u_tent_params), b.size(), iter))

            # Pressure correction
            b = rhs_p_corr(bcs=bcp)
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            self.timer("p_corr construct rhs")
            iter = solver_p_corr.solve(p_next.vector(), b)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p_next.vector())
            self.timer("p_corr solve (%s, %d dofs, %d iter)"%(', '.join(solver_p_params), b.size(), iter))

            # Velocity correction
            for d in dims:
                b = rhs_u_corr[d](bcs=bcu[d])
                self.timer("u_corr construct rhs")
                iter = solver_u_corr.solve(u_next[d].vector(), b)
                self.timer("u_corr solve (%s, %d dofs, %d iter)"%(', '.join(solver_u_corr_params), b.size(), iter))

            # Update
            self.update(problem, t, u_next, p_next)
            for d in dims: u_prev[d].assign(u_curr[d])
            for d in dims: u_curr[d].assign(u_next[d])
            p_curr.assign(p_next)

        return self._list_or_function(u_next), p_curr

    def __str__(self):
        name = "IPCS_stable"
        return name
