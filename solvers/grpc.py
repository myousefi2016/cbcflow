__author__ = "Kent-Andre Mardal <kent-and@simula.no>"
__date__ = "2008-08-05"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Harish Narayanan, 2009.
# Modified by Anders Logg, 2009-2010.

from solverbase import *
from numpy import linspace
class Solver(SolverBase):
    "cG(2/1)cG(1) with generalized Richardson iteration on the Schur complement."

    def __init__(self, options):
        assert not options['segregated']
        SolverBase.__init__(self, options)

    def select_timestep(self, problem):
        dt, t, trange = SolverBase.select_timestep(self, problem)
        if str(problem) in ["Channel"]:
            dt /= 3.0
            n = int(trange[-1] / dt + 1.0)
            dt = trange[-1] / n
            t = dt
            trange = linspace(0,trange[-1],n+1)[1:]
            if MPI.process_number() == 0:
                print "Number of time steps increased to", len(trange)
        return dt, t, trange

    def solve(self, problem, restart=None):

        # Get mesh and time step
        mesh = problem.mesh
        dt, t, trange = self.select_timestep(problem)
        if restart:
            dt, t, trange = restart.select_timestep(dt, problem.T)

        # Parameters for Uzawa iteration
        tol = problem.tolerance(problem)

        #FIXME, channel want go further down.
        tol = 1.0e-11 # slightly stricter critera to avoid drivencavity to dance
        maxiter = 100

        # Gives good convergence on drivencavity and channel (tau_y2 is not important since the fluid is not very viscous)
        tau_y1 = 2.0
        tau_y2 = 2.0

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)

        # Remove boundary stress term is problem is periodic
        if is_periodic(bcp):
            beta = Constant(0)
        else:
            beta = Constant(1)

        # FIXME: This should be moved to problem.boundary_conditions()
        pbar = problem.pressure_bc(Q)

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)
        h = CellSize(mesh)

        # Functions
        u0  = interpolate(u0, V)
        u1  = interpolate(u0, V)
        p01 = interpolate(p0, Q)
        k   = Constant(dt)
        n   = FacetNormal(mesh)
        f   = problem.f[0]
        nu  = problem.nu

        # FIXME: Compute residuals more efficiently by matrix-vector multiplication

        # Velocity and pressure residuals
        U = 0.5*(u0 + u1)
        P = p01
        Ru = inner(v, u1 - u0)*dx + k*inner(v, (grad(U)*U))*dx \
           + k*inner(epsilon(v), sigma(U, P, nu))*dx \
           - beta*k*nu*inner(v, grad(U).T*n)*ds + k*inner(v, pbar*n)*ds \
           - k*inner(v, f)*dx
        Rp = k*q*div(U)*dx

        # Add stabilization in case of Pq - Pq
        if V.ufl_element().degree() == Q.ufl_element().degree():
            d2 = Constant(0.1)*h*h
            Rp += d2*k*inner(grad(q), grad(P))*dx

        # Assemble preconditioners
        ax  = inner(v, u)*dx + 0.5*k*inner(v, (grad(u)*u0))*dx \
            + 0.5*k*2*nu*inner(epsilon(v), epsilon(u))*dx \
            - 0.5*k*nu*inner(v, grad(u).T*n)*ds

        ay1 = k**2*(inner(grad(q), grad(p)))*dx
        ay2 = k**2*((1.0/(nu*k)) * q*p)*dx

        Kx        = assemble(ax, bcs=bcu)
        Ky1, Ky1a = symmetric_assemble(ay1, bcs=bcp)
        Ky2, Ky2a = symmetric_assemble(ay2, bcs=bcp)

        # Create solvers
        Kx.solver  = LinearSolver("gmres", "jacobi")
        Ky1.solver = LinearSolver("cg", "jacobi")
        Ky2.solver = LinearSolver("cg", "jacobi")
        for K in [Kx, Ky1, Ky2]:
            K.solver.set_operator(K)
            K.solver.parameters["preconditioner"]["reuse"] = True

        # Get solution vectors
        x = u1.vector()
        y = p01.vector()
        delta_x = Vector(x) # copies parallel layout; zero'd later
        delta_y = Vector(y)

        # Time loop
        self.start_timing()
        for t in trange:

            bcu, bcp = problem.boundary_conditions(V, Q, t)
            [bc.apply(x) for bc in bcu]
            [bc.apply(y) for bc in bcp]
            [bc.homogenize() for bc in bcu+bcp]

            # GRPC iteration
            for iter in range(-1, maxiter):

                # Velocity residual
                rx = assemble(Ru, bcs=bcu)

                # Check for convergence
                if iter >= 0:
                    r = sqrt(norm(rx)**2 + norm(ry)**2)
                    if has_converged(r, iter, "GRPC", maxiter, tol): break

                # Velocity update
                delta_x.zero()
                Kx.solver.solve(delta_x, rx)
                x.axpy(-1.0, delta_x) #x -= delta_x

                # Pressure residual
                ry = assemble(Rp, bcs=bcp)

                # Pressure update 1
                ry1 = ry - Ky1a*ry
                delta_y.zero()
                if is_periodic(bcp):
                    solve(Ky1, delta_y, ry1)
                else:
                    Ky1.solver.solve(delta_y, ry1)
                if len(bcp) == 0 or is_periodic(bcp): normalize(delta_y)
                y.axpy(-tau_y1, delta_y) #y -= tau_y1*delta_y

                # Pressure update 2
                ry2 = ry - Ky2a*ry
                delta_y.zero()
                Ky2.solver.solve(delta_y, ry2)
                y.axpy(-tau_y2, delta_y) #y -= tau_y2*delta_y


            # Update
            # FIXME: Might be inaccurate for functionals of the pressure since p01 is the value at the midpoint
            self.update(problem, t, u1, p01)
            u0.assign(u1)

        return u1, p01

    def __str__(self):
        return "GRPC"
