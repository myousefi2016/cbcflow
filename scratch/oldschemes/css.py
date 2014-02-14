from __future__ import division

_author__ = "Kristian Valen-Sendstad <kvs@simula.no>"


from cbcflow.core.nsscheme import *
from cbcflow.core.timesteps import compute_regular_timesteps
from cbcflow.core.utils import Timer

class CSS(NSScheme):
    "Consistent splitting scheme by Guermond and Shen."

    def __init__(self, params):
        NSScheme.__init__(self, params, segregated=False)

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            order=2,
            u_degree=2,
            p_degree=1,
            )
        return params

    def solve(self, problem, update):
        # Get problem data
        mesh = problem.mesh
        dt, timesteps = compute_regular_timesteps(problem)
        t = timesteps[0]
        dx = problem.dx
        ds = problem.ds

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", self.params.u_degree)
        Q = FunctionSpace(mesh, "CG", self.params.p_degree)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)
        bcpsi = homogenize(bcp)
        pbar = problem.pressure_bc(Q)

        # Remove boundary stress term is problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u0  = interpolate(u0, V)
        u1  = interpolate(u0, V)
        p0  = interpolate(p0, Q)
        p1  = interpolate(p0, Q)
        p2  = interpolate(p0, Q)
        nu  = Constant(problem.params.mu/problem.params.rho)
        k   = Constant(dt)
        n   = FacetNormal(mesh)
        #f   = as_vector(problem.body_force())[0]
        f   = as_vector(problem.body_force())
        psi = Function(Q)

        # Tentative pressure
        if self.params.scheme.order == 1:
            ps = p1
        else:
            ps = 2*p1 - p0

        # Tentative velocity step
        F1 = (1/k)*inner(v, u - u0)*dx() + inner(v, grad(u0)*u0)*dx() \
            + inner(epsilon(v), sigma(u, ps, nu))*dx() \
            - beta*nu*inner(v, grad(u).T*n)*ds() + inner(v, pbar*n)*ds() \
            - inner(v, f)*dx()
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Pressure correction
        a2 = inner(grad(q), grad(p))*dx()
        L2 = (1/k)*inner(grad(q), u1 - u0)*dx() - (1/k)*inner(q*n, u1 - u0)*ds()

        # Pressure update
        a3 = q*p*dx()
        L3 = q*(ps + psi - nu*div(u1))*dx()

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Time loop
        timer = Timer(self.params.enable_timer)
        for timestep in xrange(1,len(timesteps)):
            t = timesteps[timestep]

            # Get boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)

            # Compute tentative velocity step
            b = assemble(L1)
            [bc.apply(A1, b) for bc in bcu]
            solve(A1, u1.vector(), b, "gmres", "ilu")

            # Compute pressure correction
            b = assemble(L2)
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            [bc.apply(A2, b) for bc in bcpsi]
            if is_periodic(bcp):
                solve(A2, psi.vector(), b)
            else:
                solve(A2, psi.vector(), b, "gmres", "hypre_amg")
            if len(bcp) == 0 or is_periodic(bcp): normalize(psi.vector())

            # Compute updated pressure
            b = assemble(L3)
            if len(bcp) == 0: normalize(b)
            [bc.apply(A3, b) for bc in bcp]
            solve(A3, p2.vector(), b, "gmres", "ilu")

            # Update postprocessing
            update(u1, p1, t, timestep)

            # Rotate functions for next timestep
            u0.assign(u1)
            p0.assign(p1)
            p1.assign(p2)

