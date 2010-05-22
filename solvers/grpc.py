__author__ = "Kent-Andre Mardal <kent-and@simula.no>"
__date__ = "2008-08-05"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Harish Narayanan, 2009.
# Modified by Anders Logg, 2009-2010.

from solverbase import *

class Solver(SolverBase):
    "cG(2/1)cG(1) with generalized Richardson iteration on the Schur complement."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        # Get mesh and time step
        mesh = problem.mesh
        dt, t, trange = problem.timestep(problem)

        # Parameters for Uzawa iteration
        tol = problem.tolerance(problem)

        if str(problem)=="Aneurysm":                                                                                                                                             
            pc = "jacobi"                                                                                                                                                                     
        else:                                                                                                                                                            
            pc = "ilu"  


        #FIXME, channel want go further down.
#        tol = 1.0e-6 # slightly stricter critera to avoid drivencavity to dance
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

        # FIXME: Move this somewhere else
        if mesh.topology().dim() == 2:
            shape = triangle
        else:
            shape = tetrahedron

        # Functions
        u0  = interpolate(u0, V)
        u1  = interpolate(u0, V)
        p01 = interpolate(p0, Q)
        k   = Constant(dt)
        n   = shape.n
        f   = problem.f
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
            beta  = Constant(0.1)
            d2 = beta*h*h
            Rp = Rp + d2*k*inner(grad(q), grad(P))*dx

        # Assemble preconditioners
        ax  = inner(v, u)*dx + 0.5*k*inner(v, (grad(u)*u0))*dx \
            + 0.5*k*2*nu*inner(epsilon(v), epsilon(u))*dx \
            - 0.5*k*nu*inner(v, grad(u).T*n)*ds

        ay1 = k**2*(inner(grad(q), grad(p)))*dx
        ay2 = k**2*((1.0/(nu*k)) * q*p)*dx

        Kx  = assemble(ax)
        Ky1 = assemble(ay1)
        Ky2 = assemble(ay2)
        [bc.apply(Kx) for bc in bcu]
        [bc.apply(Ky1) for bc in bcp]

        # Get solution vectors
        x = u1.vector()
        y = p01.vector()
        delta_x = Vector(x.size())
        delta_y = Vector(y.size())

        # Time loop
        self.start_timing()
        for t in trange:
            bcu, bcp = problem.boundary_conditions(V, Q, t)
            # GRPC iteration
            for iter in range(maxiter):

                # Velocity update
                rx = assemble(Ru)
                [bc.apply(rx, x) for bc in bcu]
                delta_x.zero()
                solve(Kx, delta_x, rx, 'gmres', pc)
                x.axpy(-1.0, delta_x)

                # Pressure update
                ry = assemble(Rp)
                delta_y.zero()
                if is_periodic(bcp):
                    solve(Ky1, delta_y, ry)
                else:
                    solve(Ky1, delta_y, ry, 'cg', 'amg_hypre')
                if len(bcp) == 0 or is_periodic(bcp): normalize(delta_y)
                y.axpy(-tau_y1, delta_y)

                delta_y.zero()
                solve(Ky2, delta_y, ry, 'cg', 'jacobi')
                y.axpy(-tau_y2, delta_y)

                # FIXME: A stricter convergence test should reassmble the residuals here

                # Check for convergence
                r = sqrt(norm(rx)**2 + norm(ry)**2)
                if has_converged(r, iter, "GRPC", maxiter, tol): break

            # Update
            # FIXME: Might be inaccurate for functionals of the pressure since p01 is the value at the midpoint
            self.update(problem, t, u1, p01)
            u0.assign(u1)

        return u1, p01

    def __str__(self):
        return "GRPC"
