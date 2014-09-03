from __future__ import division



from cbcflow.core.nsscheme import *
from cbcflow.core.timesteps import compute_regular_timesteps
from cbcflow.core.utils import Timer

class IPCS_p1p1_seg(NSScheme):
    "Incremental pressure-correction scheme, using P1-P1 elements, segregated version."

    def __init__(self, params):
        NSScheme.__init__(self, params, segregated=True)

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            u_degree=1,
            p_degree=1,
            )
        return params

    def solve(self, update, problem):
        # Get problem parameters
        mesh = problem.mesh
        dt, timesteps = compute_regular_timesteps(problem)
        t = timesteps[0]
        dx = problem.dx
        ds = problem.ds
        dims = range(mesh.topology().dim())

        # Define function spaces
        V = FunctionSpace(mesh, "CG", self.params.u_degree)
        Q = FunctionSpace(mesh, "CG", self.params.p_degree)

        # Get initial conditions
        u0, p0 = problem.initial_conditions(V, Q)
        u0 = as_vector([project(_, V) for _ in u0])
        p0 = project(p0, Q)

        # Get boundary conditions
        bcu, bcp = self.fetch_bcs(problem, V, Q, t)

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u1 = as_vector([Function(V) for _u0 in u0])
        p1 = interpolate(p0, Q)
        controls = []
        nu = problem.kinematic_viscosity(controls)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(V))
        n  = FacetNormal(mesh)

        # Tentative velocity step
        F_u_tent = []
        for d in dims:
            u_mean = 0.5 * (u + u0[d])
            u_diff = (u - u0[d])
            F_u_tent += [(1/k) * inner(v, u_diff) * dx()
                         + v * sum(u0[r]*u0[d].dx(r) for r in dims) * dx()
                         + inner(grad(v), nu*grad(u_mean)) * dx()
                         - v.dx(d) * p0 * dx()
                         + v * p0 * n[d] * ds()
                         - v * f[d] * dx()]

        a_u_tent = [lhs(F) for F in F_u_tent]
        L_u_tent = [rhs(F) for F in F_u_tent]

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx()
        L_p_corr = inner(grad(q), grad(p0))*dx() - (1/k)*q*sum(u1[r].dx(r) for r in dims)*dx()

        # Velocity correction
        a_u_corr = [inner(v, u)*dx() for r in dims]
        L_u_corr = [v*u1[r]*dx() - k*inner(v, grad(p1-p0)[r])*dx() for r in dims]

        # Assemble matrices
        A_u_tent = [assemble(a) for a in a_u_tent]
        A_p_corr = assemble(a_p_corr)
        A_u_corr = [assemble(a) for a in a_u_corr]

        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0:
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        # Time loop
        timer = Timer(self.params.enable_timer)
        for timestep in xrange(1,len(timesteps)):
            t = timesteps[timestep]

            # Get boundary conditions
            bcs = problem.boundary_conditions(V, Q, t)
            bcu, bcp = bcs[:-1], bcs[-1]
            timer.completed("update & fetch bc")

            # Compute tentative velocity step
            for A, L, u1_comp, bcu_comp in zip(A_u_tent, L_u_tent, u1, bcu):
                b = assemble(L)
                for bc in bcu_comp: bc.apply(A, b)
                timer.completed("u1 construct rhs")

                solver_params = self.params.solver_u_tent
                iter = solve(A, u1_comp.vector(), b, *solver_params)
                timer.completed("u1 solve (%s, %d, %d)"%(', '.join(solver_params), A.size(0), iter))

            # Pressure correction
            b = assemble(L_p_corr)
            if len(bcp) == 0 or is_periodic(bcp):
                normalize(b)
            for bc in bcp: bc.apply(A_p_corr, b)
            timer.completed("p construct rhs")

            iter = solve(A_p_corr, p1.vector(), b, *solver_p_params)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            timer.completed("p solve (%s, %d, %d)"%(', '.join(solver_p_params), A_p_corr.size(0), iter))

            # Velocity correction
            for A, L, u1_comp, bcu_comp in zip(A_u_corr, L_u_corr, u1, bcu):
                b = assemble(L)
                for bc in bcu_comp: bc.apply(A, b)
                timer.completed("u2 construct rhs")

                solver_params = self.params.solver_u_corr
                iter = solve(A, u1_comp.vector(), b, *solver_params)
                timer.completed("u2 solve (%s, %d, %d)"%(', '.join(solver_params), A.size(0),iter))

            # Update postprocessing
            update(u1, p1, t, timestep)

            # Rotate functions for next timestep
            for r in dims: u0[r].assign(u1[r])
            p0.assign(p1)
