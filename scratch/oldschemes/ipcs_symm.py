from __future__ import division


from cbcflow.core.nsscheme import *
from cbcflow.core.rhsgenerator import *
from cbcflow.core.timesteps import compute_regular_timesteps
from cbcflow.core.utils import Timer

class IPCS_symm(NSScheme):
    "Incremental pressure-correction scheme, using symmetric assemble."

    def __init__(self, params):
        NSScheme.__init__(self, params, segregated=False)

    def solve(self, problem, update, restart=None):
        # Get problem parameters
        mesh = problem.mesh
        dt, timesteps = compute_regular_timesteps(problem)
        t = timesteps[0]
        if restart: # FIXME: Broken
            dt, t, timesteps = restart.select_timestep(dt, problem.T)
        dx = problem.dx
        ds = problem.ds
        dims = range(mesh.topology().dim())

        # Define function spaces
        V = FunctionSpace(mesh, "CG", self.params.u_degree)
        Q = FunctionSpace(mesh, "CG", self.params.p_degree)

        # Get initial conditions
        if restart: # FIXME: Broken
            u0 = restart.u(t, V)
            p0 = restart.p(t, Q)
        else:
            u0, p0 = problem.initial_conditions(V, Q)
            u0 = as_vector([project(_, V) for _ in u0])
            p0 = interpolate(p0, Q)

        # Get boundary conditions
        bcu, bcp = self.fetch_bcs(problem, u0, p0, t)

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions and parameters
        u1 = as_vector([Function(V) for d in dims])
        p1 = Function(Q)

        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(V))
        n  = FacetNormal(mesh)

        # Tentative velocity step
        a1 = (1/k) * inner(v, u) * dx()
        a2 = 0.5 * inner(grad(v), nu*grad(u)) * dx()
        a_conv = -v * dot(grad(u),u0) * dx()
        Kconv = Matrix() # assembled from a_conv in the time loop

        # Create the final coefficient matrix for the tentative velocity
        # step. We use a single matrix for all dimensions, which means that the
        # BCs must match (they must apply to the full vector space at the
        # vertex, not a subspace.
        A_u_tent, A_u_tent_asymm = symmetric_assemble(a1+a2, bcs=bcu[0]) # FIXME: Does this assumption hold now?
        A_u_tent_asymm.compress()

        # Create matrices for generating the RHS
        B = assemble(a1-a2)
        M_symm, M_asymm = symmetric_assemble(v*u*dx(), bcs=bcu[0]) # FIXME: Does this assumption hold now?
        M_asymm.compress()

        # Define how to create the RHS for the tentative velocity. The RHS is
        # (of course) different for each dimension.
        rhs_u_tent = [None]*len(dims)
        for d in dims:
            C = assemble(-v*p*n[d]*ds() + v.dx(d)*p*dx())
            rhs_u_tent[d] = RhsGenerator(V)
            rhs_u_tent[d] += B, u0[d]
            rhs_u_tent[d] += C, p0
            rhs_u_tent[d] += M_symm, f[d]
            rhs_u_tent[d] += M_asymm, f[d]
            rhs_u_tent[d] += Kconv, u0[d]

        # Pressure correction
        A_p_corr, A_p_corr_asymm = symmetric_assemble(inner(grad(q), grad(p))*dx(), bcs=bcp)
        A_p_corr_asymm.compress()
        rhs_p_corr = RhsGenerator(Q)
        rhs_p_corr += A_p_corr, p0
        rhs_p_corr += A_p_corr_asymm, p0
        for d in dims:
            Ku = assemble(-(1/k)*q*u.dx(d)*dx())
            rhs_p_corr += Ku, u1[d]

        # Velocity correction. Like for the tentative velocity, a single LHS is used.
        A_u_corr, A_u_corr_asymm = M_symm, M_asymm
        rhs_u_corr = [None]*len(dims)
        for d in dims:
            Kp = assemble(-k*inner(v, grad(p)[d])*dx())
            rhs_u_corr[d] = RhsGenerator(V)
            rhs_u_corr[d] += M_symm, u1[d]
            rhs_u_corr[d] += M_asymm, u1[d]
            rhs_u_corr[d] += Kp, p1
            rhs_u_corr[d] -= Kp, p0

        # Create solvers
        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0 or is_periodic(bcp):
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        solver_u_tent = LinearSolver(*self.params.solver_u_tent)
        solver_p_corr = LinearSolver(*solver_p_params)
        solver_u_corr = LinearSolver(*self.params.solver_u_corr)

        for A,S in [(A_u_tent, solver_u_tent), (A_p_corr, solver_p_corr), (A_u_corr, solver_u_corr)]:
            S.set_operator(A)
            if 'preconditioner' in S.parameters:
                S.parameters['preconditioner']['reuse'] = True

        # Time loop
        timer = Timer(self.params.enable_timer)
        for timestep in xrange(1,len(timesteps)):
            t = timesteps[timestep]

            # Get boundary conditions
            bcu, bcp = self.fetch_bcs(problem, u0, p0, t)
            timer.completed("update & fetch bc")

            # Assemble the u0-dependent convection matrix. It is important that
            # it is assembled into the same tensor, because it is stored in
            # rhs. (And it's faster).
            assemble(a_conv, tensor=Kconv, reset_sparsity=(Kconv.size(0)==0))
            timer.completed("assemble a_conv")

            # Compute tentative velocity step
            for d in dims:
                b = rhs_u_tent[d](bcs=bcu[d], symmetric_mod=A_u_tent_asymm)
                timer.completed("u0 construct rhs")

                iter = solver_u_tent.solve(u1[d].vector(), b)
                timer.completed("u0 solve (%s, %d dofs, %d iter)"%(', '.join(self.params.solver_u_tent), b.size(), iter))

            # Pressure correction
            b = rhs_p_corr(bcs=bcp, symmetric_mod=A_p_corr_asymm)
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            timer.completed("p1 construct rhs")

            iter = solver_p_corr.solve(p1.vector(), b)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            timer.completed("p1 solve (%s, %d dofs, %d iter)"%(', '.join(solver_p_params), b.size(), iter))

            # Velocity correction
            for d in dims:
                b = rhs_u_corr[d](bcs=bcu[d], symmetric_mod=A_u_corr_asymm)
                timer.completed("u1 construct rhs")

                iter = solver_u_corr.solve(u1[d].vector(), b)
                timer.completed("u1 solve (%s, %d dofs, %d iter)"%(', '.join(self.params.solver_u_corr), b.size(), iter))

            # Update postprocessing
            update(u1, p1, t, timestep)

            # Rotate functions for next timestep
            for d in dims: u0[d].assign(u1[d])
            p0.assign(p1)
