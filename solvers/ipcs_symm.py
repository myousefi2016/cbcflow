from __future__ import division

__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2012-02-01"
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

        solver_u_tent_params      = "cg", "hypre_euclid"
        solver_p_periodic_params  = "cg", "hypre_euclid"
        solver_p_dirichlet_params = "cg", "ml_amg"
        solver_u_corr_params      = "cg", "hypre_euclid"

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = self.select_timestep(problem)

        # Define function spaces
        V = FunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)

        # Get initial and boundary conditions
        ics = problem.initial_conditions(V, Q)
        bcs = problem.boundary_conditions(V, Q, t)
        u0, p0 = ics[:-1], ics[-1]
        bcu, bcp = bcs[:-1], bcs[-1]

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        dim  = len(u0)
        dims = range(dim)
        u0 = [interpolate(_u0, V) for _u0 in u0]
        u1 = [Function(V) for d in dims]

        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

        # Tentative velocity step
        a1 = (1/k) * inner(v, u) * dx
        a2 = 0.5 * inner(grad(v), nu*grad(u)) * dx
        a_conv = -v * sum(u0[r]*u.dx(r) for r in dims) * dx
        Kconv = Matrix() # assembled from a_conv in the time loop

        # Create the final coefficient matrix for the tentative velocity
        # step. We use a single matrix for all dimensions, which means that the
        # BCs must match (they must apply to the full vector space at the
        # vertex, not a subspace.
        A_u_tent, A_u_tent_asymm = symmetric_assemble(a1+a2, bcs=bcu[0])
        A_u_tent_asymm.compress()

        # Create matrices for generating the RHS
        B = assemble(a1-a2)
        M_symm, M_asymm = symmetric_assemble(v*u*dx, bcs=bcu[0])
        M_asymm.compress()

        # Define how to create the RHS for the tentative velocity. The RHS is
        # (of course) different for each dimension.
        rhs_u_tent = [None]*dim
        for d in dims:
            C = assemble(-v*p*n[d]*ds + v.dx(d)*p*dx)
            rhs_u_tent[d] = RhsGenerator(V)
            rhs_u_tent[d] += B, u0[d]
            rhs_u_tent[d] += C, p0
            rhs_u_tent[d] += M_symm, f[d]
            rhs_u_tent[d] += M_asymm, f[d]
            rhs_u_tent[d] += Kconv, u0[d]

        # Pressure correction
        A_p_corr, A_p_corr_asymm = symmetric_assemble(inner(grad(q), grad(p))*dx, bcs=bcp)
        A_p_corr_asymm.compress()
        rhs_p_corr = RhsGenerator(Q)
        rhs_p_corr += A_p_corr, p0
        rhs_p_corr += A_p_corr_asymm, p0
        for d in dims:
            Ku = assemble(-(1/k)*q*u.dx(d)*dx)
            rhs_p_corr += Ku, u1[d]

        # Velocity correction. Like for the tentative velocity, a single LHS is used.
        A_u_corr, A_u_corr_asymm = M_symm, M_asymm
        rhs_u_corr = [None]*dim
        for d in dims:
            Kp = assemble(-k*inner(v, grad(p)[d])*dx)
            rhs_u_corr[d] = RhsGenerator(V)
            rhs_u_corr[d] += M_symm, u1[d]
            rhs_u_corr[d] += M_asymm, u1[d]
            rhs_u_corr[d] += Kp, p1
            rhs_u_corr[d] -= Kp, p0

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

            # Assemble the u0-dependent convection matrix. It is important that
            # it is assembled into the same tensor, because it is stored in
            # rhs. (And it's faster).
            assemble(a_conv, tensor=Kconv, reset_sparsity=(Kconv.size(0)==0))

            # Compute tentative velocity step
            for d in dims:
                b = rhs_u_tent[d](bcs=bcu[d], symmetric_mod=A_u_tent_asymm)
                self.timer("u0 construct rhs")
                iter = solver_u_tent.solve(u1[d].vector(), b)
                self.timer("u0 solve (%s, %d dofs, %d iter)"%(', '.join(solver_u_tent_params), b.size(), iter))

            # Pressure correction
            b = rhs_p_corr(bcs=bcp, symmetric_mod=A_p_corr_asymm)
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            self.timer("p1 construct rhs")
            iter = solver_p_corr.solve(p1.vector(), b)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            self.timer("p1 solve (%s, %d dofs, %d iter)"%(', '.join(solver_p_params), b.size(), iter))

            # Velocity correction
            for d in dims:
                b = rhs_u_corr[d](bcs=bcu[d], symmetric_mod=A_u_corr_asymm)
                self.timer("u1 construct rhs")
                iter = solver_u_corr.solve(u1[d].vector(), b)
                self.timer("u1 solve (%s, %d dofs, %d iter)"%(', '.join(solver_u_corr_params), b.size(), iter))

            # Update
            self.update(problem, t, u1, p1)
            for d in dims: u0[d].assign(u1[d])
            p0.assign(p1)

        return self._list_or_function(u1), p1

    def __str__(self):
        name = "IPCS_symm"
        return name
