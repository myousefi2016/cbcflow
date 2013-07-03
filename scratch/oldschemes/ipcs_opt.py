from __future__ import division

__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2011-11-11"
__copyright__ = "Copyright (C) 2011 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# TODO: Import only what's needed here and elsewhere:
from ..core.nsscheme import *
from ..core.rhsgenerator import *
from ..core.utils import Timer, is_periodic, epsilon
from ..core.timesteps import compute_regular_timesteps

class IPCS_opt(NSScheme):
    "Incremental pressure-correction scheme, optimized version."

    def __init__(self, params):
        NSScheme.__init__(self, params, segregated=False)

    def solve(self, problem, update, restart=None):
        # Get problem parameters
        mesh = problem.mesh
        dt, timesteps = compute_regular_timesteps(problem)
        t = timesteps[0]
        if restart: # FIXME: broken
            dt, t, timesteps = restart.select_timestep(dt, problem.T)
        dx = problem.dx
        ds = problem.ds
        dims = range(mesh.topology().dim())

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", self.params.u_degree)
        Q = FunctionSpace(mesh, "CG", self.params.p_degree)

        # Get initial conditions
        if restart: # FIXME: Broken
            u0 = restart.u(t, V)
            p0 = restart.p(t, Q)
        else:
            u0, p0 = problem.initial_conditions(V, Q)
            u0 = as_vector([interpolate(_, V) for _ in u0])
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

        # Functions
        u1 = as_vector([Function(V) for d in dims])
        p1 = Function(Q)

        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force())
        n  = FacetNormal(mesh)

        # To avoid indexing in non-segregated forms
        u0_ = u0[0]
        u1_ = u1[0]
        f_  = f[0]

        # Tentative velocity step
        M  = assemble(inner(v, u) * dx())
        K1 = assemble((1/k) * inner(v, u) * dx())
        K2 = assemble(inner(epsilon(v), nu*epsilon(u)) * dx()
                      - 0.5 * beta * nu * inner(v, grad(u).T*n) * ds())
        A = K1+K2
        K3 = assemble(-inner(v, p*n)*ds() + div(v)*p*dx())

        rhs = RhsGenerator(V)
        rhs += K1, u0_
        rhs -= K2, u0_
        rhs += K3, p0
        rhs += M, f_
        rhs += -inner(v, grad(u0_)*u0_) * dx()

        A_u_tent, rhs_u_tent = [A], [rhs]

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx())
        Ku = assemble(-(1/k)*q*div(u)*dx())
        rhs_p_corr = RhsGenerator(Q)
        rhs_p_corr += A_p_corr, p0
        rhs_p_corr += Ku, u1_

        # Velocity correction
        A_u_corr = [M.copy() for d in dims]
        Kp = assemble(-k*inner(v, grad(p))*dx())
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

        # Create solvers
        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0:
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        solver_u_tent = [LinearSolver(*self.params.solver_u_tent) for d in dims]
        solver_p_corr = LinearSolver(*solver_p_params)
        solver_u_corr = [LinearSolver(*self.params.solver_u_corr) for d in dims]

        for A,S in zip(A_u_tent, solver_u_tent) \
                + [(A_p_corr, solver_p_corr)] \
                + zip(A_u_corr, solver_u_corr):
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

            # Compute tentative velocity step
            for d, S, rhs, u1_comp, bcu_comp in zip(dims, solver_u_tent, rhs_u_tent, u1, bcu):
                b = rhs()
                for bc in bcu_comp: bc.apply(b)
                timer.completed("u0 construct rhs")

                iter = S.solve(u1_comp.vector(), b)
                timer.completed("u0 solve (%s, %d, %d)"%(', '.join(self.params.solver_u_tent), A.size(0), iter))

            # Pressure correction
            b = rhs_p_corr()
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            for bc in bcp: bc.apply(b)
            timer.completed("p1 construct rhs")

            iter = solver_p_corr.solve(p1.vector(), b)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())
            timer.completed("p1 solve (%s, %d, %d)"%(', '.join(solver_p_params), A_p_corr.size(0), iter))

            # Velocity correction
            for S, rhs, u1_comp, bcu_comp in zip(solver_u_corr, rhs_u_corr, u1, bcu):
                b = rhs()
                for bc in bcu_comp: bc.apply(b)
                timer.completed("u1 construct rhs")

                iter = S.solve(u1_comp.vector(), b)
                timer.completed("u1 solve (%s, %d, %d)"%(', '.join(self.params.solver_u_corr), A.size(0),iter))

            # Update postprocessing
            update(u1, p1, t, timestep)

            # Rotate functions for next timestep
            for r in dims: u0[r].assign(u1[r])
            p0.assign(p1)
